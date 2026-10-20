import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import zipfile

from scripts.add_source_links import add_source_links


class SourceLinksTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.archive_path = self.root / "pmo-en-macos.zip"
        self.manifest_path = self.root / "manifest.json"
        info = zipfile.ZipInfo("launch.command", (2024, 3, 4, 5, 6, 8))
        info.create_system = 3
        info.external_attr = (0o100755 << 16) | 0x10
        info.compress_type = zipfile.ZIP_DEFLATED
        with zipfile.ZipFile(self.archive_path, "w") as archive:
            archive.writestr(info, b"#!/bin/sh\necho ready\n")
            archive.writestr("course/readme.txt", b"course")
        self.original_bytes = self.archive_path.read_bytes()
        self.original_sha256 = hashlib.sha256(self.original_bytes).hexdigest()
        self.manifest = {
            "archive": self.archive_path.name,
            "bytes": len(self.original_bytes),
            "sha256": self.original_sha256,
            "validation": {"platform": "macos-arm64", "source_commit": "a" * 40},
        }
        self._write_manifest()

    def tearDown(self):
        self.temp.cleanup()

    def _write_manifest(self):
        self.manifest_path.write_text(json.dumps(self.manifest), encoding="utf-8")

    def test_append_preserves_members_metadata_and_updates_manifest(self):
        with zipfile.ZipFile(self.archive_path) as before:
            old = {
                info.filename: (hashlib.sha256(before.read(info.filename)).hexdigest(), info.date_time,
                                info.compress_type, info.create_system, info.external_attr)
                for info in before.infolist()
            }
        add_source_links(self.archive_path, self.manifest_path)
        with zipfile.ZipFile(self.archive_path) as after:
            self.assertIsNone(after.testzip())
            self.assertEqual(set(after.namelist()), set(old) | {
                "pmo-en-macos/resources/licenses/SOURCE_DOWNLOAD.txt"
            })
            for info in after.infolist():
                if info.filename == "pmo-en-macos/resources/licenses/SOURCE_DOWNLOAD.txt":
                    continue
                self.assertEqual(
                    (hashlib.sha256(after.read(info.filename)).hexdigest(), info.date_time,
                     info.compress_type, info.create_system, info.external_attr),
                    old[info.filename],
                )
            command = after.getinfo("launch.command")
            self.assertEqual((command.external_attr >> 16) & 0o777, 0o755)
            text = after.read("pmo-en-macos/resources/licenses/SOURCE_DOWNLOAD.txt").decode()
            self.assertIn("sources-macos-arm64-index.json", text)
            self.assertIn("/tree/" + "a" * 40, text)
            self.assertIn("not required before students run", text)
        updated = json.loads(self.manifest_path.read_text(encoding="utf-8"))
        self.assertEqual(updated["validation"], self.manifest["validation"])
        self.assertEqual(updated["bytes"], self.archive_path.stat().st_size)
        self.assertEqual(updated["sha256"], hashlib.sha256(self.archive_path.read_bytes()).hexdigest())
        addition = updated["documentation_addition"]
        self.assertEqual(addition["original_sha256"], self.original_sha256)
        self.assertEqual(addition["original_bytes"], len(self.original_bytes))
        self.assertEqual(addition["original_member_count"], 2)
        self.assertTrue(addition["unchanged"])

    def test_bad_sha256_does_not_change_archive_or_manifest(self):
        self.manifest["sha256"] = "0" * 64
        self._write_manifest()
        manifest_bytes = self.manifest_path.read_bytes()
        with self.assertRaisesRegex(ValueError, "SHA256"):
            add_source_links(self.archive_path, self.manifest_path)
        self.assertEqual(self.archive_path.read_bytes(), self.original_bytes)
        self.assertEqual(self.manifest_path.read_bytes(), manifest_bytes)

    def test_second_append_is_rejected(self):
        add_source_links(self.archive_path, self.manifest_path)
        first_archive = self.archive_path.read_bytes()
        first_manifest = self.manifest_path.read_bytes()
        with self.assertRaisesRegex(ValueError, "already contains"):
            add_source_links(self.archive_path, self.manifest_path)
        self.assertEqual(self.archive_path.read_bytes(), first_archive)
        self.assertEqual(self.manifest_path.read_bytes(), first_manifest)

    def test_cli_supplement_preserves_source_provenance_and_rejects_repeat(self):
        add_source_links(self.archive_path, self.manifest_path)
        first_manifest = json.loads(self.manifest_path.read_text(encoding="utf-8"))
        first_sha256 = hashlib.sha256(self.archive_path.read_bytes()).hexdigest()
        first_bytes = self.archive_path.stat().st_size
        supplement_dir = self.root / "license-documents"
        supplement_dir.mkdir()
        documents = {
            "LICENSE.txt": b"DLL licence text\n",
            "NOTICE.txt": b"DLL attribution\n",
        }
        for name, content in documents.items():
            (supplement_dir / name).write_bytes(content)

        script = Path(__file__).resolve().parents[1] / "scripts" / "add_source_links.py"
        command = [sys.executable, str(script), str(self.archive_path), str(self.manifest_path),
                   "--supplement-dir", str(supplement_dir)]
        completed = subprocess.run(command, capture_output=True, text=True)
        self.assertEqual(completed.returncode, 0, completed.stderr)

        updated = json.loads(self.manifest_path.read_text(encoding="utf-8"))
        self.assertEqual(updated["documentation_addition"], first_manifest["documentation_addition"])
        self.assertEqual(updated["validation"], first_manifest["validation"])
        supplement = updated["licence_supplement"]
        self.assertEqual(supplement["original_sha256"], first_sha256)
        self.assertEqual(supplement["original_bytes"], first_bytes)
        self.assertEqual(supplement["original_member_count"], 3)
        self.assertTrue(supplement["unchanged"])
        expected_files = {
            f"pmo-en-macos/resources/licenses/{name}": hashlib.sha256(content).hexdigest()
            for name, content in documents.items()
        }
        self.assertEqual(supplement["files"], expected_files)
        with zipfile.ZipFile(self.archive_path) as archive:
            self.assertIsNone(archive.testzip())
            self.assertEqual((archive.getinfo("launch.command").external_attr >> 16) & 0o777, 0o755)
            self.assertIn("pmo-en-macos/resources/licenses/SOURCE_DOWNLOAD.txt", archive.namelist())
            for name, content in documents.items():
                self.assertEqual(
                    archive.read(f"pmo-en-macos/resources/licenses/{name}"), content
                )

        supplemented_archive = self.archive_path.read_bytes()
        supplemented_manifest = self.manifest_path.read_bytes()
        repeated = subprocess.run(command, capture_output=True, text=True)
        self.assertNotEqual(repeated.returncode, 0)
        self.assertEqual(self.archive_path.read_bytes(), supplemented_archive)
        self.assertEqual(self.manifest_path.read_bytes(), supplemented_manifest)


if __name__ == "__main__":
    unittest.main()
