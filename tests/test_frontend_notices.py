import base64
import hashlib
import io
import json
from pathlib import Path
import sys
import tarfile
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from collect_frontend_notices import collect_frontend


class FrontendNoticeTests(unittest.TestCase):
    def run_collection(self, name, version, acquire):
        temporary = tempfile.TemporaryDirectory()
        root = Path(temporary.name)
        notices = root / "notices"
        sources = root / "sources"
        catalog = notices / "jupyter" / "third-party-licenses.json"
        catalog.parent.mkdir(parents=True)
        catalog.write_text(
            json.dumps({"packages": [{"name": name, "versionInfo": version, "licenseId": "MIT"}]}),
            encoding="utf-8",
        )
        archive = io.BytesIO()
        with tarfile.open(fileobj=archive, mode="w:gz") as tar:
            payload = b"package without a notice"
            member = tarfile.TarInfo("package/index.js")
            member.size = len(payload)
            tar.addfile(member, io.BytesIO(payload))
        archive_bytes = archive.getvalue()
        integrity = "sha512-" + base64.b64encode(hashlib.sha512(archive_bytes).digest()).decode()
        tarball_url = "https://registry.npmjs.org/package.tgz"
        metadata = {
            "name": name,
            "version": version,
            "dist": {"tarball": tarball_url, "integrity": integrity},
            "repository": "https://registry.npmjs.org/package",
        }

        def acquire_with_tarball(url, destination, expected=None):
            if url == tarball_url:
                destination.parent.mkdir(parents=True, exist_ok=True)
                destination.write_bytes(archive_bytes)
                return {"url": url, "path": destination.relative_to(sources).as_posix()}
            return acquire(url, destination, expected)

        with patch(
            "collect_frontend_notices.urllib.request.urlopen",
            return_value=io.BytesIO(json.dumps(metadata).encode()),
        ):
            result = collect_frontend(notices, sources, acquire_with_tarball)
        return temporary, notices, sources, result

    def test_uses_exact_local_pin_and_records_relative_source(self):
        temporary, notices, _, result = self.run_collection(
            "json-schema-merge-allof", "0.8.1", lambda *_: self.fail("local pin must not use network")
        )
        self.addCleanup(temporary.cleanup)
        row = result["packages"][0]
        self.assertEqual(row["status"], "captured")
        pin = next(
            item for item in json.loads((ROOT / "config/frontend_license_sources.json").read_text())["packages"]
            if item["name"] == "json-schema-merge-allof" and item["version"] == "0.8.1"
        )
        self.assertEqual(row["repository"], pin["repository"])
        notice = row["notices"][0]
        self.assertEqual(notice["path"], "npm/json-schema-merge-allof-0.8.1/upstream/LICENSE-with-attribution.txt")
        self.assertEqual(notice["sha256"], pin["notices"][0]["sha256"])
        self.assertEqual(notice["source"], {"local_file": pin["notices"][0]["local_file"]})
        self.assertTrue((notices / notice["path"]).is_file())
        self.assertNotIn(str(notices.parent), json.dumps(result))

    def test_url_pin_hash_mismatch_is_not_captured(self):
        calls = []

        def wrong_hash(url, destination, expected):
            calls.append((url, destination, expected))
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes(b"wrong content")
            self.assertEqual(destination.name, expected + ".txt")
            self.assertNotEqual(hashlib.sha256(destination.read_bytes()).hexdigest(), expected)
            raise ValueError("Source checksum mismatch")

        temporary, _, sources, result = self.run_collection("@jupyter/react-components", "0.16.6", wrong_hash)
        self.addCleanup(temporary.cleanup)
        row = result["packages"][0]
        self.assertEqual(row["status"], "collection_incomplete")
        self.assertNotIn("notices", row)
        self.assertEqual(len(calls), 1)
        url, destination, expected = calls[0]
        self.assertTrue(url.startswith("https://raw.githubusercontent.com/"))
        self.assertEqual(destination, sources / "npm-licences" / f"{expected}.txt")


if __name__ == "__main__":
    unittest.main()
