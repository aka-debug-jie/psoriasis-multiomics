"""Recover missing Jupyter frontend notices from integrity-checked npm packages."""
from pathlib import Path, PurePosixPath
import base64
import concurrent.futures
import hashlib
import json
import tarfile
import urllib.parse
import urllib.request

from collect_component_metadata import is_notice_name, safe_part, sha256_file

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]


def collect_frontend(notices, sources, acquire):
    catalog = next(notices.glob('jupyter/**/third-party-licenses.json'))
    packages = json.loads(catalog.read_text())['packages']
    selected = [p for p in packages if not p.get('extractedText') or p.get('licenseId') == 'EPL-2.0']
    pins = json.loads((REPOSITORY_ROOT / 'config/frontend_license_sources.json').read_text())['packages']
    pinned = {(item['name'], item['version']): item for item in pins}

    def process(package):
        name, version = package['name'], package['versionInfo']
        row = {'name': name, 'version': version, 'license': package['licenseId']}
        try:
            url = 'https://registry.npmjs.org/' + urllib.parse.quote(name, safe='') + '/' + version
            with urllib.request.urlopen(url, timeout=60) as response:
                metadata = json.load(response)
            assert metadata['name'] == name and metadata['version'] == version
            folder = safe_part(name) + '-' + safe_part(version)
            archive = sources / 'npm' / (folder + '.tgz')
            source = acquire(metadata['dist']['tarball'], archive)
            integrity = metadata['dist']['integrity'].split()[0]
            algorithm, encoded = integrity.split('-', 1)
            assert algorithm in ('sha512', 'sha256'), 'Unsupported npm integrity algorithm'
            with archive.open('rb') as handle:
                actual = hashlib.file_digest(handle, algorithm).digest()
            assert base64.b64encode(actual).decode() == encoded, 'npm integrity mismatch'
            captured = []
            with tarfile.open(archive) as tar:
                for member in tar.getmembers():
                    path = PurePosixPath(member.name)
                    if not member.isfile() or not is_notice_name(path.name):
                        continue
                    assert path.parts[0] == 'package' and '..' not in path.parts
                    assert '\\' not in member.name and ':' not in member.name
                    target = notices / 'npm' / folder / path
                    target.parent.mkdir(parents=True, exist_ok=True)
                    target.write_bytes(tar.extractfile(member).read())
                    captured.append({'path': target.relative_to(notices).as_posix(),
                                     'sha256': sha256_file(target)})
            pin = pinned.get((name, version))
            fallback_used = False
            if not captured and pin:
                for item in pin['notices']:
                    expected = item['sha256']
                    notice_name = PurePosixPath(item['name'])
                    assert notice_name.name == item['name'] and notice_name.name not in ('.', '..')
                    if item.get('url'):
                        cache = sources / 'npm-licences' / (expected + '.txt')
                        acquired = acquire(item['url'], cache, expected)
                        source_record = {'url': item['url'], 'cache_path': acquired['path']}
                        content = cache.read_bytes()
                    else:
                        local_file = PurePosixPath(item['local_file'])
                        assert not local_file.is_absolute() and '..' not in local_file.parts
                        source_path = REPOSITORY_ROOT.joinpath(*local_file.parts)
                        content = source_path.read_bytes()
                        actual = hashlib.sha256(content).hexdigest()
                        assert actual == expected, 'Pinned local notice checksum mismatch'
                        source_record = {'local_file': local_file.as_posix()}
                    actual = hashlib.sha256(content).hexdigest()
                    assert actual == expected, 'Pinned upstream notice checksum mismatch'
                    target = notices / 'npm' / folder / 'upstream' / notice_name.name
                    target.parent.mkdir(parents=True, exist_ok=True)
                    target.write_bytes(content)
                    captured.append({'path': target.relative_to(notices).as_posix(),
                                     'sha256': actual, 'source': source_record})
                fallback_used = bool(captured)
            row.update(registry=url, integrity=integrity, source=source,
                       repository=(pin.get('repository') if fallback_used else metadata.get('repository')),
                       git_head=metadata.get('gitHead'),
                       notices=captured, status='captured' if captured else 'notice_not_found')
        except Exception as error:
            row.update(status='collection_incomplete', error=str(error))
        return row

    with concurrent.futures.ThreadPoolExecutor(max_workers=6) as pool:
        rows = list(pool.map(process, selected))
    result = {'source_catalog': catalog.relative_to(notices).as_posix(),
              'catalog_packages': len(packages), 'review_status': 'pending_manual_review',
              'packages': rows}
    (notices / 'npm-components.json').write_text(json.dumps(result, indent=2) + '\n')
    return result
