import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


class BootstrapEnvTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.project = Path(self.temp.name) / "project"
        research = self.project / "research"
        (research / "scripts").mkdir(parents=True)
        (research / "records").mkdir()
        (research / "requirements-data.txt").write_text("demo-package==1.2.3\n", encoding="utf-8")
        (research / "records/python_version.txt").write_text("Python 3.12.7\n", encoding="utf-8")
        (research / "records/conda_explicit.txt").write_text("@EXPLICIT\nhttps://example.test/r.tar.bz2\n", encoding="utf-8")
        (research / "renv.lock").write_text('{"R":{"Version":"4.4.2"},"Packages":{}}\n', encoding="utf-8")
        shutil.copy2(ROOT / "research/scripts/bootstrap_env.sh", research / "scripts/bootstrap_env.sh")
        shutil.copy2(ROOT / "research/scripts/bootstrap_data_env.sh", research / "scripts/bootstrap_data_env.sh")
        self.stub_bin = Path(self.temp.name) / "bin"
        self.stub_bin.mkdir()
        self.log = Path(self.temp.name) / "calls.log"
        self._write_command_stubs()

    def tearDown(self):
        self.temp.cleanup()

    def _write_executable(self, path, content):
        path.write_text(content, encoding="utf-8")
        path.chmod(0o755)

    def _write_rscript(self, path):
        self._write_executable(path, """#!/usr/bin/env bash
printf 'Rscript %s RENV_PATHS_ROOT=%s RENV_PATHS_CACHE=%s\\n' "$*" "${RENV_PATHS_ROOT:-}" "${RENV_PATHS_CACHE:-}" >> "$CALL_LOG"
if [[ "${R_MODE:-ok}" == drift && "$*" == *" -"* ]]; then
  echo 'locked package version mismatch: demo expected 1.0 found 2.0' >&2
  exit 19
fi
""")

    def _write_python(self, path):
        self._write_executable(path, """#!/usr/bin/env bash
printf 'env-python %s\\n' "$*" >> "$CALL_LOG"
if [[ "${1:-}" == --version ]]; then echo 'Python 3.12.7'; exit 0; fi
if [[ "${1:-}" == - ]]; then
  cat >/dev/null
  if [[ "${PY_MODE:-ok}" == version-drift ]]; then echo 'Python version mismatch: expected Python 3.12.7, found 3.12.8' >&2; exit 17; fi
  if [[ "${PY_MODE:-ok}" == drift ]]; then echo 'demo-package expected 1.2.3, found 9.9.9' >&2; exit 17; fi
  exit 0
fi
if [[ "${1:-}" == -m && "${2:-}" == pip ]]; then exit 0; fi
exit 0
""")

    def _write_command_stubs(self):
        self._write_executable(self.stub_bin / "conda", """#!/usr/bin/env bash
printf 'conda %s\\n' "$*" >> "$CALL_LOG"
if [[ "${1:-}" == create ]]; then
  prefix=''
  while (($#)); do
    if [[ "$1" == --prefix ]]; then prefix="$2"; shift 2; else shift; fi
  done
  mkdir -p "$prefix/bin"
  cat > "$prefix/bin/Rscript" <<'RS'
#!/usr/bin/env bash
printf 'Rscript %s RENV_PATHS_ROOT=%s RENV_PATHS_CACHE=%s\\n' "$*" "${RENV_PATHS_ROOT:-}" "${RENV_PATHS_CACHE:-}" >> "$CALL_LOG"
RS
  chmod +x "$prefix/bin/Rscript"
fi
""")
        self._write_executable(self.stub_bin / "python3.12", """#!/usr/bin/env bash
printf 'bootstrap-python %s\\n' "$*" >> "$CALL_LOG"
if [[ "${1:-}" == --version ]]; then echo 'Python 3.12.7'; exit 0; fi
if [[ "${1:-}" == -m && "${2:-}" == venv ]]; then
  mkdir -p "$3/bin"
  cat > "$3/bin/python" <<'PY'
#!/usr/bin/env bash
printf 'env-python %s\\n' "$*" >> "$CALL_LOG"
if [[ "${1:-}" == --version ]]; then echo 'Python 3.12.7'; exit 0; fi
if [[ "${1:-}" == - ]]; then cat >/dev/null; exit 0; fi
if [[ "${1:-}" == -m && "${2:-}" == pip ]]; then exit 0; fi
exit 0
PY
  chmod +x "$3/bin/python"
fi
""")

    def _environment(self, **overrides):
        env = os.environ.copy()
        env.update({"PATH": f"{self.stub_bin}:{env['PATH']}", "CALL_LOG": str(self.log)})
        env.update(overrides)
        return env

    def _run(self, script, **overrides):
        return subprocess.run(
            ["bash", str(self.project / "research/scripts" / script)],
            cwd=self.project,
            env=self._environment(**overrides),
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )

    def test_existing_r_environment_only_validates_and_preserves_locks(self):
        rscript = self.project / "research/envs/r442/bin/Rscript"
        rscript.parent.mkdir(parents=True)
        self._write_rscript(rscript)
        lock = self.project / "research/records/conda_explicit.txt"
        lock_before = lock.read_bytes()
        renv_before = (self.project / "research/renv.lock").read_bytes()

        result = self._run("bootstrap_env.sh")

        self.assertEqual(result.returncode, 0, result.stderr)
        calls = self.log.read_text(encoding="utf-8")
        self.assertIn("Rscript --vanilla -", calls)
        self.assertNotIn("conda ", calls)
        self.assertEqual(lock.read_bytes(), lock_before)
        self.assertEqual((self.project / "research/renv.lock").read_bytes(), renv_before)
        self.assertFalse((self.project / "research/records/r_packages.csv").exists())

    def test_r_version_drift_fails_without_restore(self):
        rscript = self.project / "research/envs/r442/bin/Rscript"
        rscript.parent.mkdir(parents=True)
        self._write_rscript(rscript)

        result = self._run("bootstrap_env.sh", R_MODE="drift")

        self.assertNotEqual(result.returncode, 0)
        calls = self.log.read_text(encoding="utf-8")
        self.assertNotIn("conda ", calls)
        self.assertNotIn("renv::restore", calls)

    def test_new_r_environment_uses_both_existing_locks(self):
        lock = self.project / "research/records/conda_explicit.txt"
        lock_before = lock.read_bytes()
        renv_before = (self.project / "research/renv.lock").read_bytes()

        result = self._run("bootstrap_env.sh")

        self.assertEqual(result.returncode, 0, result.stderr)
        calls = self.log.read_text(encoding="utf-8")
        self.assertIn(f"conda create --yes --prefix {self.project}/research/envs/r442 --file {lock}", calls)
        self.assertIn('renv::restore(project="research", library=.Library, lockfile="research/renv.lock", prompt=FALSE)', calls)
        self.assertIn(f"RENV_PATHS_ROOT={self.project}/research/cache/renv-root", calls)
        self.assertIn(f"RENV_PATHS_CACHE={self.project}/research/cache/renv", calls)
        self.assertEqual(lock.read_bytes(), lock_before)
        self.assertEqual((self.project / "research/renv.lock").read_bytes(), renv_before)

    def test_existing_python_environment_validates_without_reinstalling(self):
        python = self.project / "research/envs/data/bin/python"
        python.parent.mkdir(parents=True)
        self._write_python(python)

        result = self._run("bootstrap_data_env.sh")

        self.assertEqual(result.returncode, 0, result.stderr)
        calls = self.log.read_text(encoding="utf-8")
        self.assertIn("env-python -", calls)
        self.assertNotIn("pip install", calls)
        self.assertNotIn("bootstrap-python", calls)
        self.assertFalse((self.project / "research/records/python_freeze.txt").exists())

    def test_python_package_drift_fails_without_reinstall(self):
        python = self.project / "research/envs/data/bin/python"
        python.parent.mkdir(parents=True)
        self._write_python(python)

        result = self._run("bootstrap_data_env.sh", PY_MODE="drift")

        self.assertNotEqual(result.returncode, 0)
        self.assertIn("expected 1.2.3", result.stderr)
        calls = self.log.read_text(encoding="utf-8")
        self.assertNotIn("pip install", calls)

    def test_python_version_drift_fails_without_reinstall(self):
        python = self.project / "research/envs/data/bin/python"
        python.parent.mkdir(parents=True)
        self._write_python(python)

        result = self._run("bootstrap_data_env.sh", PY_MODE="version-drift")

        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Python version mismatch", result.stderr)
        calls = self.log.read_text(encoding="utf-8")
        self.assertNotIn("pip install", calls)

    def test_new_python_environment_uses_version_and_pinned_requirements(self):
        result = self._run("bootstrap_data_env.sh")

        self.assertEqual(result.returncode, 0, result.stderr)
        calls = self.log.read_text(encoding="utf-8")
        self.assertIn("bootstrap-python --version", calls)
        self.assertIn("bootstrap-python -m venv", calls)
        self.assertIn("env-python -m pip install -r", calls)
        self.assertIn("env-python -", calls)


if __name__ == "__main__":
    unittest.main()
