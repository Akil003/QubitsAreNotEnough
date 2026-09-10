"""Regenerate PDF_PREFLIGHT.txt from the built PDFs.

PDF_PREFLIGHT.txt previously had no generator and drifted out of date (it reported 19
pages against a 23-page manuscript). This script rebuilds it from whatever PDFs are
present, so the preflight record is reproducible like every other artefact in the package.

Type 3 fonts are reported explicitly because IEEE rejects them; matplotlib emits Type 3
by default, which is why the plotting scripts set ``pdf.fonttype = 42``.

Usage:  python make_preflight.py
"""
from __future__ import annotations

import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent
TARGETS = ("Qubits_Are_Not_Enough_TQE.pdf", "Qubits_Are_Not_Enough_TQE_IEEEtran.pdf")


def _pdffonts(pdf: Path) -> tuple[int, int, int]:
    """Return (total fonts, Type 3 count, non-embedded count)."""
    out = subprocess.run(["pdffonts", str(pdf)], capture_output=True, text=True).stdout
    rows = [r for r in out.splitlines()[2:] if r.strip()]
    type3 = sum(1 for r in rows if "Type 3" in r)
    not_emb = sum(1 for r in rows if r.split() and r.split()[-1] == "no")
    return len(rows), type3, not_emb


def _pdfinfo(pdf: Path) -> dict[str, str]:
    out = subprocess.run(["pdfinfo", str(pdf)], capture_output=True, text=True).stdout
    info = {}
    for line in out.splitlines():
        if ":" in line:
            k, v = line.split(":", 1)
            info[k.strip()] = v.strip()
    return info


def main() -> None:
    lines: list[str] = []
    for name in TARGETS:
        pdf = ROOT / name
        if not pdf.exists():
            continue
        info = _pdfinfo(pdf)
        total, type3, not_emb = _pdffonts(pdf)
        lines += [
            f"PDF: {name}",
            f"Pages: {info.get('Pages', '?')}",
            f"Page size: {info.get('Page size', '?')}",
            f"PDF version: {info.get('PDF version', '?')}",
            f"Encrypted: {info.get('Encrypted', '?')}",
            f"Fonts: {total}",
            f"Type 3 fonts: {type3}   (IEEE requires 0)",
            f"Non-embedded fonts: {not_emb}   (IEEE requires 0)",
            "",
        ]
    (ROOT / "PDF_PREFLIGHT.txt").write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
