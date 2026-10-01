#!/usr/bin/env python3
"""Refresh SHA-256 hashes in MANIFEST.sha256 and add final repository helper scripts."""

from __future__ import annotations

import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parent
MANIFEST = ROOT / "MANIFEST.sha256"
EXTRA_TRACKED = ("make_submission_archive.py", "refresh_manifest.py")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    entries: set[str] = set()
    for raw in MANIFEST.read_text(encoding="utf-8").splitlines():
        if not raw.strip():
            continue
        try:
            _, rel = raw.split("  ", 1)
        except ValueError as exc:
            raise SystemExit(f"Malformed manifest line: {raw!r}") from exc
        entries.add(rel.removeprefix("./"))

    for rel in EXTRA_TRACKED:
        if (ROOT / rel).is_file():
            entries.add(rel)

    rows = []
    for rel in sorted(entries):
        path = ROOT / rel
        if not path.is_file():
            raise SystemExit(f"Manifest path is missing: {rel}")
        rows.append(f"{sha256(path)}  ./{rel}")

    MANIFEST.write_text("\n".join(rows) + "\n", encoding="utf-8")
    print(f"Refreshed {len(rows)} manifest entries.")


if __name__ == "__main__":
    main()
