import hashlib
import json
from pathlib import Path
import tempfile
import unittest

from scripts.collect_component_metadata import collect


class ComponentMetadataTests(unittest.TestCase):
    def test_records_runtime_dll_hashes_only_for_selected_conda_packages(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            runtime = root / "runtime"
            output = root / "metadata-output"
            metadata = runtime / "conda-meta"
            metadata.mkdir(parents=True)
            selected = ("vc14_runtime", "vcomp14", "ucrt")
            expected = {}
            for name in (*selected, "demo"):
                relative = f"Library/bin/{name}.dll"
                binary = runtime / relative
                binary.parent.mkdir(parents=True, exist_ok=True)
                content = f"fixture binary for {name}".encode()
                binary.write_bytes(content)
                if name in selected:
                    expected[name] = {
                        "path": relative,
                        "bytes": len(content),
                        "sha256": hashlib.sha256(content).hexdigest(),
                    }
                (metadata / f"{name}.json").write_text(
                    json.dumps(
                        {
                            "name": name,
                            "version": "1.0",
                            "build": "0",
                            "files": [relative, "Library/bin/readme.txt", "missing.dll"],
                        }
                    ),
                    encoding="utf-8",
                )

            manifest = collect(runtime, output)
            components = {item["name"]: item for item in manifest["components"] if item["ecosystem"] == "conda"}
            self.assertEqual(set(components), {*selected, "demo"})
            for name in selected:
                self.assertEqual(components[name]["binary_files"], [expected[name]])
            self.assertNotIn("binary_files", components["demo"])
            self.assertNotIn(str(runtime), json.dumps(manifest))

    def test_collects_three_ecosystems_notices_and_jupyter_file_without_leaking_paths(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            runtime = root / "runtime"
            output = root / "metadata-output"
            package_source = runtime / "pkgs" / "demo-1.0-0"
            (runtime / "conda-meta").mkdir(parents=True)
            (package_source / "info" / "licenses").mkdir(parents=True)
            (package_source / "info" / "recipe").mkdir(parents=True)
            conda_license = package_source / "info" / "licenses" / "LICENSE.txt"
            recipe = package_source / "info" / "recipe" / "meta.yaml"
            conda_license.write_text("conda license\n", encoding="utf-8")
            recipe.write_text("about:\n  license: BSD-3-Clause\n", encoding="utf-8")
            (runtime / "conda-meta" / "demo.json").write_text(
                json.dumps(
                    {
                        "name": "demo",
                        "version": "1.0",
                        "build": "py_0",
                        "license": "BSD-3-Clause",
                        "url": f"file://{package_source}/demo-1.0.tar.bz2",
                        "sha256": "a" * 64,
                        "link": {"source": str(package_source)},
                    }
                ),
                encoding="utf-8",
            )

            r_package = runtime / "course-library" / "rDemo"
            r_package.mkdir(parents=True)
            (r_package / "DESCRIPTION").write_text(
                "Package: rDemo\nVersion: 0.2.1\nLicense: MIT + file LICENSE\n"
                "Repository: CRAN\nRemoteType: cran\nRemoteRepo: rDemo\n"
                "RemoteSha: abcdef123456\nBiocViews: Example, Testing\n",
                encoding="utf-8",
            )
            (r_package / "LICENSE").write_text("R package license\n", encoding="utf-8")

            dist_info = runtime / "lib" / "python3.12" / "site-packages" / "pyDemo-2.3.dist-info"
            (dist_info / "licenses").mkdir(parents=True)
            (dist_info / "METADATA").write_text(
                "Metadata-Version: 2.3\nName: pyDemo\nVersion: 2.3\n"
                "License-Expression: MIT\nHome-page: https://example.invalid/pyDemo\n\n",
                encoding="utf-8",
            )
            py_license = dist_info / "licenses" / "NOTICE"
            py_license.write_text("Python package notice\n", encoding="utf-8")

            jupyter_license = runtime / "share" / "jupyter" / "lab" / "static" / "third-party-licenses.json"
            jupyter_license.parent.mkdir(parents=True)
            jupyter_license.write_text('[{"name":"ui-dependency","license":"MIT"}]\n', encoding="utf-8")

            outside_notice = root / "outside-secret.txt"
            outside_notice.write_text("must not be copied\n", encoding="utf-8")
            try:
                (package_source / "info" / "licenses" / "EXTERNAL_NOTICE").symlink_to(outside_notice)
            except (OSError, NotImplementedError):
                pass

            runtime_snapshot = {
                path.relative_to(runtime).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
                for path in runtime.rglob("*")
                if path.is_file() and not path.is_symlink()
            }
            manifest = collect(runtime, output)
            self.assertEqual(
                runtime_snapshot,
                {
                    path.relative_to(runtime).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
                    for path in runtime.rglob("*")
                    if path.is_file() and not path.is_symlink()
                },
            )
            saved = json.loads((output / "components.json").read_text(encoding="utf-8"))
            self.assertEqual(manifest, saved)
            by_ecosystem = {item["ecosystem"]: item for item in saved["components"]}
            self.assertTrue({"conda", "r", "python", "jupyter"}.issubset(by_ecosystem))
            self.assertEqual(by_ecosystem["conda"]["build"], "py_0")
            self.assertEqual(by_ecosystem["conda"]["artifact_sha256"], "a" * 64)
            self.assertIsNone(by_ecosystem["conda"]["artifact_url"])
            self.assertEqual(by_ecosystem["r"]["source_metadata"]["RemoteSha"], "abcdef123456")
            self.assertEqual(by_ecosystem["python"]["source_metadata"]["Home-page"], "https://example.invalid/pyDemo")
            for component in saved["components"]:
                self.assertEqual(component["source_status"], "unreviewed")
                for notice in component["notices"]:
                    notice_path = output / notice["path"]
                    self.assertTrue(notice_path.is_file())
                    self.assertEqual(
                        hashlib.sha256(notice_path.read_bytes()).hexdigest(), notice["sha256"]
                    )
                    self.assertFalse(Path(notice["path"]).is_absolute())
            self.assertTrue(
                (output / "jupyter/share/jupyter/lab/static/third-party-licenses.json").is_file()
            )
            self.assertFalse(
                (output / "conda/demo-1.0-py_0/info/licenses/EXTERNAL_NOTICE").exists()
            )
            serialized = json.dumps(saved, ensure_ascii=False)
            self.assertNotIn(str(runtime), serialized)
            self.assertNotIn(str(package_source), serialized)
            self.assertNotIn(str(outside_notice), serialized)

    def test_rejects_output_inside_runtime(self):
        with tempfile.TemporaryDirectory() as temporary:
            runtime = Path(temporary) / "runtime"
            runtime.mkdir()
            with self.assertRaises(ValueError):
                collect(runtime, runtime / "output")


if __name__ == "__main__":
    unittest.main()
