"""Install the pinned official PLINK binary locally; no analysis is run."""
from pathlib import Path
import hashlib
import io
import json
import subprocess
import zipfile
import requests

ROOT = Path(__file__).resolve().parents[2]
URL = "https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20260927.zip"


def main():
    target = ROOT / "research/tools/plink-1.9.0"
    target.mkdir(parents=True, exist_ok=True)
    archive = target / "plink_linux_x86_64_20260927.zip"
    if not archive.exists():
        response = requests.get(URL, timeout=60)
        response.raise_for_status()
        archive.write_bytes(response.content)
    with zipfile.ZipFile(archive) as z:
        for name in z.namelist():
            if "/" not in name and name not in (".", ".."):
                (target / name).write_bytes(z.read(name))
    binary = target / "plink"
    binary.chmod(0o755)
    version = subprocess.check_output([str(binary), "--version"], text=True).strip()
    record = {"url": URL, "archive_sha256": hashlib.sha256(archive.read_bytes()).hexdigest(),
              "binary_sha256": hashlib.sha256(binary.read_bytes()).hexdigest(), "version": version,
              "paper_version": "not specified; installed version is not claimed to be the author's version",
              "smr": "not installed: official Yang Lab host was unreachable during acquisition"}
    (ROOT / "research/records/external_tools.json").write_text(json.dumps(record, indent=2) + "\n")
    print(version)


if __name__ == "__main__":
    main()
