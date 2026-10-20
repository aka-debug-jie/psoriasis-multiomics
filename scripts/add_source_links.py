#!/usr/bin/env python3
"""Append source-download documentation to a verified release archive."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import tempfile
import zipfile


RELEASE = "https://github.com/aka-debug-jie/psoriasis-multiomics/releases/tag/v4.0.0"
RELEASE_DOWNLOAD = "https://github.com/aka-debug-jie/psoriasis-multiomics/releases/download/v4.0.0"
REPOSITORY = "https://github.com/aka-debug-jie/psoriasis-multiomics"
CHUNK_SIZE = 1024 * 1024


def _sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for chunk in iter(lambda: stream.read(CHUNK_SIZE), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _member_hash(archive, name):
    digest = hashlib.sha256()
    with archive.open(name, "r") as stream:
        for chunk in iter(lambda: stream.read(CHUNK_SIZE), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _zip_info_signature(info):
    """Metadata that should survive a central-directory rewrite."""
    return (
        info.filename,
        info.orig_filename,
        info.date_time,
        info.compress_type,
        info.create_system,
        info.create_version,
        info.extract_version,
        info.reserved,
        info.external_attr,
        info.internal_attr,
        info.flag_bits,
        info.volume,
        info.CRC,
        info.compress_size,
        info.file_size,
        info.header_offset,
        info.extra,
        info.comment,
    )


def _source_text(platform, source_commit):
    index_name = f"sources-{platform}-index.json"
    return (
        "Source code downloads\n"
        "======================\n\n"
        f"Release: {RELEASE}\n"
        f"Source archive index: {RELEASE_DOWNLOAD}/{index_name}\n"
        f"Source tree at the validated commit: {REPOSITORY}/tree/{source_commit}\n\n"
        "Download the individual source ZIP files listed in the index and verify "
        "each file against its SHA256 value. Each source ZIP is a separate file; "
        "do not concatenate the ZIP files. Source downloads are optional and are "
        "not required before students run the course.\n"
    ).encode("utf-8")


def add_source_links(archive_path, manifest_path):
    archive_path = Path(archive_path)
    manifest_path = Path(manifest_path)
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if not isinstance(manifest, dict):
        raise ValueError("manifest must be a JSON object")
    validation = manifest.get("validation")
    if not isinstance(validation, dict):
        raise ValueError("manifest validation must be an object")
    platform = validation.get("platform")
    source_commit = validation.get("source_commit")
    if not isinstance(platform, str) or not platform or not isinstance(source_commit, str) or not source_commit:
        raise ValueError("manifest validation must contain platform and source_commit")
    if manifest.get("archive") != archive_path.name:
        raise ValueError("manifest archive does not match the input archive filename")
    if manifest.get("bytes") != archive_path.stat().st_size:
        raise ValueError("archive size does not match manifest")
    original_sha256 = _sha256_file(archive_path)
    if manifest.get("sha256") != original_sha256:
        raise ValueError("archive SHA256 does not match manifest")

    added_name = f"{archive_path.stem}/resources/licenses/SOURCE_DOWNLOAD.txt"
    added_bytes = _source_text(platform, source_commit)
    with zipfile.ZipFile(archive_path, "r") as original:
        infos = original.infolist()
        names = [info.filename for info in infos]
        if len(names) != len(set(names)):
            raise ValueError("archive contains duplicate member names")
        if added_name in names:
            raise ValueError(f"archive already contains {added_name}")
        original_members = {
            info.filename: (_member_hash(original, info.filename), _zip_info_signature(info))
            for info in infos
        }

    temp_dir = Path(tempfile.mkdtemp(prefix=".source-links-", dir=archive_path.parent))
    temp_archive = temp_dir / archive_path.name
    temp_manifest = temp_dir / manifest_path.name
    backup_archive = temp_dir / "original.zip"
    backup_manifest = temp_dir / "original-manifest.json"
    archive_replaced = False
    manifest_replaced = False
    try:
        shutil.copy2(archive_path, temp_archive)
        with zipfile.ZipFile(temp_archive, "a") as output:
            info = zipfile.ZipInfo(added_name)
            info.compress_type = zipfile.ZIP_DEFLATED
            output.writestr(info, added_bytes)

        with zipfile.ZipFile(temp_archive, "r") as updated:
            updated_infos = updated.infolist()
            updated_names = [info.filename for info in updated_infos]
            if len(updated_names) != len(set(updated_names)):
                raise ValueError("updated archive contains duplicate member names")
            if updated.testzip() is not None:
                raise ValueError("updated archive failed ZIP integrity testing")
            if set(updated_names) != set(names) | {added_name}:
                raise ValueError("updated archive member names differ from expected")
            for info in updated_infos:
                if info.filename == added_name:
                    continue
                if (_member_hash(updated, info.filename), _zip_info_signature(info)) != original_members[info.filename]:
                    raise ValueError(f"existing ZIP member changed: {info.filename}")
            if _member_hash(updated, added_name) != hashlib.sha256(added_bytes).hexdigest():
                raise ValueError("added documentation content failed SHA256 verification")

        updated_sha256 = _sha256_file(temp_archive)
        updated_bytes = temp_archive.stat().st_size
        updated_manifest = dict(manifest)
        updated_manifest["bytes"] = updated_bytes
        updated_manifest["sha256"] = updated_sha256
        updated_manifest["documentation_addition"] = {
            "original_bytes": archive_path.stat().st_size,
            "original_sha256": original_sha256,
            "file": added_name,
            "file_sha256": hashlib.sha256(added_bytes).hexdigest(),
            "original_member_count": len(infos),
            "unchanged": True,
            "unchanged_evidence": "All original member names, streamed SHA256 values, and ZIP metadata were verified after append.",
        }
        temp_manifest.write_text(json.dumps(updated_manifest, indent=2) + "\n", encoding="utf-8")

        shutil.copy2(archive_path, backup_archive)
        shutil.copy2(manifest_path, backup_manifest)
        os.replace(temp_archive, archive_path)
        archive_replaced = True
        os.replace(temp_manifest, manifest_path)
        manifest_replaced = True
    except Exception:
        if archive_replaced:
            os.replace(backup_archive, archive_path)
        if manifest_replaced:
            os.replace(backup_manifest, manifest_path)
        raise
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", help="release ZIP archive to update")
    parser.add_argument("manifest", help="companion JSON manifest")
    args = parser.parse_args()
    try:
        add_source_links(args.archive, args.manifest)
    except (OSError, ValueError, KeyError, zipfile.BadZipFile, json.JSONDecodeError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
