#!/usr/bin/env python3
"""Collect installed component metadata and license notices without copying payloads."""
from __future__ import annotations

import argparse
import email.parser
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import sys
from typing import Any, Iterable
from urllib.parse import urlsplit

NOTICE_NAMES = ("license", "licence", "copying", "notice")
R_FIELDS = (
    "Package",
    "Version",
    "License",
    "Repository",
    "RemoteType",
    "RemoteRepo",
    "RemoteSha",
    "BiocViews",
    "RemoteUsername",
    "RemoteHost",
    "RemoteRef",
    "Priority",
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def safe_part(value: Any) -> str:
    """Make metadata suitable for a relative output directory name."""
    text = str(value or "unknown")
    return re.sub(r"[^A-Za-z0-9._+-]+", "_", text).strip("._") or "unknown"


def local_url(value: Any) -> str | None:
    """Keep only literal HTTP(S) artifact URLs; never expose local paths."""
    if not isinstance(value, str):
        return None
    value = value.strip()
    parsed = urlsplit(value)
    if parsed.scheme in ("https", "http") and parsed.hostname and not parsed.username and not parsed.password:
        return value
    return None


def safe_metadata_value(value: str | None) -> str | None:
    """Retain raw metadata except values that would expose a local machine path."""
    if value is None:
        return None
    text = value.strip()
    if text.lower().startswith("file://") or re.match(r"^(?:/|\\\\|[A-Za-z]:[\\/])", text):
        return None
    return value


def ensure_output_directory(path: Path, output: Path) -> None:
    """Create output parents without traversing pre-existing symlink directories."""
    relative = path.relative_to(output)
    current = output
    for part in relative.parts:
        current = current / part
        if current.is_symlink():
            raise ValueError("output contains a symlink directory")
        current.mkdir(exist_ok=True)


def has_symlink_component(path: Path) -> bool:
    """Return whether any existing component in a source path is a symlink."""
    current = Path(path.anchor)
    for part in path.parts[1:]:
        current = current / part
        if current.is_symlink():
            return True
    return False


def parse_dcf(path: Path) -> dict[str, str]:
    """Parse the simple field/continuation format used by R DESCRIPTION files."""
    fields: dict[str, str] = {}
    current: str | None = None
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line[:1] in (" ", "\t") and current:
            fields[current] += "\n" + line.strip()
        elif ":" in line:
            key, value = line.split(":", 1)
            current = key.strip()
            fields[current] = value.strip()
        else:
            current = None
    return fields


def iter_regular_files(root: Path) -> Iterable[Path]:
    """Walk without following symlinks and yield regular files only."""
    if root.is_symlink() or not root.is_dir():
        return
    for current, dirs, files in os.walk(root, followlinks=False):
        current_path = Path(current)
        dirs[:] = [name for name in dirs if not (current_path / name).is_symlink()]
        for name in files:
            path = current_path / name
            if not path.is_symlink() and path.is_file():
                yield path


def is_notice_name(name: str) -> bool:
    lowered = name.lower()
    return any(
        lowered == stem or lowered.startswith(stem + ".") or lowered.startswith(stem + "-")
        for stem in NOTICE_NAMES
    )


def copy_notice_files(
    sources: Iterable[Path], output: Path, destination: Path
) -> list[dict[str, str]]:
    """Copy selected regular files and return output-relative path/hash records."""
    notices: list[dict[str, str]] = []
    seen: set[str] = set()
    for source_root in sources:
        for source_file in iter_regular_files(source_root):
            if not is_notice_name(source_file.name):
                continue
            relative = source_file.relative_to(source_root)
            # Keep files with identical relative names from separate roots distinct.
            relative_key = relative.as_posix()
            if relative_key in seen:
                continue
            seen.add(relative_key)
            target = destination / "notices" / relative
            ensure_output_directory(target.parent, output)
            if target.is_symlink():
                continue
            shutil.copyfile(source_file, target)
            notices.append(
                {
                    "path": target.relative_to(output).as_posix(),
                    "sha256": sha256_file(target),
                }
            )
    return notices


def copy_tree_files(source: Path, output: Path, destination: Path) -> list[dict[str, str]]:
    """Copy all regular files from a metadata-only source tree, skipping symlinks."""
    notices: list[dict[str, str]] = []
    for source_file in iter_regular_files(source):
        relative = source_file.relative_to(source)
        target = destination / relative
        ensure_output_directory(target.parent, output)
        if target.is_symlink():
            continue
        shutil.copyfile(source_file, target)
        notices.append(
            {"path": target.relative_to(output).as_posix(), "sha256": sha256_file(target)}
        )
    return notices


def conda_components(runtime: Path, output: Path) -> list[dict[str, Any]]:
    components: list[dict[str, Any]] = []
    for record_path in sorted((runtime / "conda-meta").glob("*.json")):
        try:
            record = json.loads(record_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        name = record.get("name")
        version = record.get("version")
        build = record.get("build")
        if not name or not version:
            continue
        directory = output / "conda" / "-".join(
            safe_part(value) for value in (name, version, build or "unknown")
        )
        link = record.get("link")
        source_value = link.get("source") if isinstance(link, dict) else None
        notices: list[dict[str, str]] = []
        if isinstance(source_value, str) and source_value:
            package_source = Path(source_value).expanduser()
            # link.source is used only for local lookup; it is never serialized.
            if not has_symlink_component(package_source):
                for metadata_subdir in ("info/licenses", "info/recipe"):
                    metadata_source = package_source / metadata_subdir
                    if not has_symlink_component(metadata_source):
                        notices.extend(
                            copy_tree_files(metadata_source, output, directory / metadata_subdir)
                        )
        components.append(
            {
                "ecosystem": "conda",
                "name": name,
                "version": version,
                "build": build,
                "license": safe_metadata_value(record.get("license")),
                "artifact_url": local_url(record.get("url")),
                "artifact_sha256": record.get("sha256"),
                "notices": notices,
                "source_status": "unreviewed",
            }
        )
    return components


def r_library_roots(runtime: Path) -> list[Path]:
    return [runtime / "course-library", runtime / "R" / "library", runtime / "lib" / "R" / "library"]


def r_components(runtime: Path, output: Path) -> list[dict[str, Any]]:
    components: list[dict[str, Any]] = []
    seen: set[Path] = set()
    for library in r_library_roots(runtime):
        for description in sorted(library.glob("*/DESCRIPTION")):
            if description.is_symlink() or not description.is_file():
                continue
            package_dir = description.parent
            resolved = package_dir.resolve()
            if resolved in seen:
                continue
            seen.add(resolved)
            fields = parse_dcf(description)
            name = fields.get("Package")
            version = fields.get("Version")
            if not name or not version:
                continue
            destination = output / "r" / f"{safe_part(name)}-{safe_part(version)}"
            notices = copy_notice_files([package_dir], output, destination)
            components.append(
                {
                    "ecosystem": "r",
                    "name": name,
                    "version": version,
                    "build": None,
                    "license": safe_metadata_value(fields.get("License")),
                    "artifact_url": None,
                    "artifact_sha256": None,
                    "source_metadata": {key: safe_metadata_value(fields.get(key)) for key in R_FIELDS},
                    "notices": notices,
                    "source_status": "unreviewed",
                }
            )
    return components


def python_site_packages(runtime: Path) -> list[Path]:
    roots = [runtime / "Lib" / "site-packages"]
    roots.extend(sorted((runtime / "lib").glob("python*/site-packages")))
    return roots


def python_components(runtime: Path, output: Path) -> list[dict[str, Any]]:
    components: list[dict[str, Any]] = []
    seen: set[Path] = set()
    for site_packages in python_site_packages(runtime):
        for metadata in sorted(site_packages.glob("*.dist-info/METADATA")):
            if metadata.is_symlink() or not metadata.is_file():
                continue
            dist_info = metadata.parent
            resolved = dist_info.resolve()
            if resolved in seen:
                continue
            seen.add(resolved)
            try:
                message = email.parser.Parser().parsestr(metadata.read_text(encoding="utf-8", errors="replace"))
            except OSError:
                continue
            name = message.get("Name")
            version = message.get("Version")
            if not name or not version:
                continue
            destination = output / "python" / f"{safe_part(name)}-{safe_part(version)}"
            # Preserve legal notices from dist-info only; do not copy package payloads.
            notice_roots = [dist_info]
            notices = copy_notice_files(notice_roots, output, destination)
            license_value = safe_metadata_value(message.get("License-Expression") or message.get("License"))
            components.append(
                {
                    "ecosystem": "python",
                    "name": name,
                    "version": version,
                    "build": None,
                    "license": license_value,
                    "artifact_url": None,
                    "artifact_sha256": None,
                    "source_metadata": {
                        "Name": name,
                        "Version": version,
                        "License-Expression": safe_metadata_value(message.get("License-Expression")),
                        "License": safe_metadata_value(message.get("License")),
                        "Home-page": local_url(message.get("Home-page")),
                    },
                    "notices": notices,
                    "source_status": "unreviewed",
                }
            )
    return components


def jupyter_license_component(runtime: Path, output: Path) -> list[dict[str, Any]]:
    candidates = (
        runtime / "share/jupyter/lab/static/third-party-licenses.json",
        runtime / "Lib/share/jupyter/lab/static/third-party-licenses.json",
        runtime / "lib/share/jupyter/lab/static/third-party-licenses.json",
    )
    source = next((path for path in candidates if path.is_file() and not path.is_symlink()), None)
    if source is None:
        return []
    destination = output / "jupyter/share/jupyter/lab/static/third-party-licenses.json"
    ensure_output_directory(destination.parent, output)
    if destination.is_symlink():
        return []
    shutil.copyfile(source, destination)
    notice = {"path": destination.relative_to(output).as_posix(), "sha256": sha256_file(destination)}
    return [
        {
            "ecosystem": "jupyter",
            "name": "third-party-licenses",
            "version": None,
            "build": None,
            "license": None,
            "artifact_url": None,
            "artifact_sha256": None,
            "notices": [notice],
            "source_status": "unreviewed",
        }
    ]


def collect(runtime: Path, output: Path) -> dict[str, Any]:
    runtime = runtime.expanduser().resolve()
    output = output.expanduser().resolve()
    try:
        output.relative_to(runtime)
    except ValueError:
        pass
    else:
        raise ValueError("output must be outside runtime so runtime remains read-only")
    if not runtime.is_dir():
        raise ValueError("runtime directory does not exist")
    output.mkdir(parents=True, exist_ok=True)
    components = []
    components.extend(conda_components(runtime, output))
    components.extend(r_components(runtime, output))
    components.extend(python_components(runtime, output))
    components.extend(jupyter_license_component(runtime, output))
    manifest = {"components": components}
    (output / "components.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runtime", required=True, type=Path, help="read-only runtime root")
    parser.add_argument("--output", required=True, type=Path, help="metadata output directory outside runtime")
    args = parser.parse_args(argv)
    try:
        manifest = collect(args.runtime, args.output)
    except ValueError as exc:
        parser.error(str(exc))
    except OSError:
        parser.error("filesystem operation failed")
    print(f"Wrote metadata for {len(manifest['components'])} components.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
