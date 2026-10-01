#!/usr/bin/env python3
"""Build and validate the exact IEEE TQE submission bundle."""

from __future__ import annotations

import hashlib
import re
import subprocess
import sys
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent
MAIN = "Qubits_Are_Not_Enough_TQE"
BUNDLE = ROOT / "tqe_submission_bundle.zip"


def run(*args: str) -> None:
    subprocess.run(args, cwd=ROOT, check=True)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def verify_manifest() -> None:
    for raw in (ROOT / "MANIFEST.sha256").read_text(encoding="utf-8").splitlines():
        if not raw.strip():
            continue
        expected, rel = raw.split("  ", 1)
        rel = rel.removeprefix("./")
        path = ROOT / rel
        if not path.is_file():
            raise SystemExit(f"Manifest path is missing: {rel}")
        actual = sha256(path)
        if actual != expected:
            raise SystemExit(
                f"Manifest mismatch for {rel}: expected {expected}, obtained {actual}"
            )


def enforce_preflight() -> None:
    text = (ROOT / "PDF_PREFLIGHT.txt").read_text(encoding="utf-8")
    if text.startswith("STATUS: STALE"):
        raise SystemExit("PDF preflight is still marked stale after rebuild.")

    type3 = [int(v) for v in re.findall(r"Type 3 fonts:\s*(\d+)", text)]
    unembedded = [int(v) for v in re.findall(r"Non-embedded fonts:\s*(\d+)", text)]
    if len(type3) < 2 or len(unembedded) < 2:
        raise SystemExit("Preflight did not report both submission PDFs.")
    if any(type3):
        raise SystemExit(f"Type-3 font gate failed: {type3}")
    if any(unembedded):
        raise SystemExit(f"Font-embedding gate failed: {unembedded}")


def referenced_figures() -> list[Path]:
    pattern = re.compile(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}")
    out: set[Path] = set()

    for source_name in (f"{MAIN}.tex", f"{MAIN}_supplementary.tex"):
        source = ROOT / source_name
        text = source.read_text(encoding="utf-8")
        for name in pattern.findall(text):
            candidate = ROOT / name
            if not candidate.is_file():
                candidate = ROOT / "figures" / name
            if not candidate.is_file():
                raise SystemExit(f"Referenced figure is missing: {name}")
            out.add(candidate)

    return sorted(out)


def main() -> None:
    run(sys.executable, "validate_outputs.py", "--full")
    run("make", "manuscript")
    run("make", "supplementary")
    run(sys.executable, "make_preflight.py")
    enforce_preflight()

    run(sys.executable, "refresh_manifest.py")
    verify_manifest()

    files = [
        ROOT / f"{MAIN}.pdf",
        ROOT / f"{MAIN}.tex",
        ROOT / f"{MAIN}.bbl",
        ROOT / f"{MAIN}_supplementary.pdf",
        ROOT / f"{MAIN}_supplementary.tex",
        ROOT / f"{MAIN}_supplementary.bbl",
        ROOT / "references.bib",
        ROOT / "IEEEtran.cls",
        ROOT / "PDF_PREFLIGHT.txt",
    ]
    files.extend(referenced_figures())

    forbidden_names = {"ieeeaccess.cls", "logo.png", "notaglinelogo.png", "bullet.png"}
    for path in files:
        if path.name in forbidden_names or "IEEEAccess" in path.name:
            raise SystemExit(f"Access-only artefact leaked into TQE package: {path.name}")
        if not path.is_file():
            raise SystemExit(f"Submission file is missing: {path.relative_to(ROOT)}")

    if BUNDLE.exists():
        BUNDLE.unlink()

    with zipfile.ZipFile(BUNDLE, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for path in files:
            zf.write(path, path.relative_to(ROOT))

    print(f"Wrote {BUNDLE.name} with {len(files)} files.")
    for path in files:
        print(f"  {path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
