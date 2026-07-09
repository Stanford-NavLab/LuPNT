#!/usr/bin/env python
"""Shrink the embedded PNG outputs in the tutorial notebooks in place.

The docs render the notebooks with their stored outputs (nbsphinx_execute = "never"),
so matplotlib figures live as base64 PNGs inside each ``ex*.ipynb``. Full-res, full-colour
PNGs make the committed notebooks large; this pass re-encodes each image output with an
adaptive palette (and an optional max-width cap), which is lossless enough for documentation
figures while cutting the byte size several-fold. Interactive Plotly / Cesium outputs are
left untouched (they are hosted externally, not embedded).

    pixi run python scripts/optimize_notebook_images.py            # optimize all ex*.ipynb
    pixi run python scripts/optimize_notebook_images.py --dry-run  # report only
    pixi run python scripts/optimize_notebook_images.py --max-width 1200 --colors 256
"""

from __future__ import annotations

import argparse
import base64
import io
import sys
from pathlib import Path

from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = ROOT / "python" / "examples"


def optimize_png(b64: str, max_width: int, colors: int) -> str | None:
    """Return a smaller base64 PNG, or None if the re-encode did not help."""
    raw = base64.b64decode(b64)
    img = Image.open(io.BytesIO(raw))
    orig_mode = img.mode
    # Flatten transparency onto white so palette conversion stays clean.
    if img.mode in ("RGBA", "LA") or (img.mode == "P" and "transparency" in img.info):
        bg = Image.new("RGB", img.size, (255, 255, 255))
        rgba = img.convert("RGBA")
        bg.paste(rgba, mask=rgba.split()[-1])
        img = bg
    else:
        img = img.convert("RGB")
    if max_width and img.width > max_width:
        h = round(img.height * max_width / img.width)
        img = img.resize((max_width, h), Image.LANCZOS)
    pal = img.convert("P", palette=Image.ADAPTIVE, colors=colors)
    out = io.BytesIO()
    pal.save(out, format="PNG", optimize=True)
    new = out.getvalue()
    if len(new) >= len(raw):  # keep the original if we did not actually shrink it
        return None
    return base64.b64encode(new).decode("ascii"), len(raw), len(new), orig_mode, img.size


def process(nb_path: Path, max_width: int, colors: int, dry_run: bool):
    import nbformat

    nb = nbformat.read(nb_path, as_version=4)
    before = after = 0
    n_imgs = 0
    for cell in nb.cells:
        if cell.get("cell_type") != "code":
            continue
        for o in cell.get("outputs", []):
            data = o.get("data")
            if not data or "image/png" not in data:
                continue
            v = data["image/png"]
            b64 = "".join(v) if isinstance(v, list) else v
            res = optimize_png(b64, max_width, colors)
            n_imgs += 1
            before += len(b64)
            if res is None:
                after += len(b64)
                continue
            new_b64, raw_len, new_len, _, _ = res
            after += len(new_b64)
            if not dry_run:
                data["image/png"] = new_b64
    if not dry_run and after < before:
        nbformat.write(nb, nb_path)
    return n_imgs, before, after


def main() -> int:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "--max-width", type=int, default=1200, help="cap image width in px (0 = no resize)"
    )
    p.add_argument("--colors", type=int, default=256, help="adaptive palette size")
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--only", nargs="+", help="substring filter on notebook stem")
    args = p.parse_args()

    nbs = sorted(EXAMPLES.glob("ex[0-9]*.ipynb"))
    if args.only:
        nbs = [n for n in nbs if any(s in n.stem for s in args.only)]

    tot_b = tot_a = 0
    print(f"{'notebook':34s} {'imgs':>4s} {'before':>9s} {'after':>9s}  saved")
    for nb in nbs:
        n, b, a = process(nb, args.max_width, args.colors, args.dry_run)
        tot_b += b
        tot_a += a
        if n:
            print(f"{nb.name:34s} {n:4d} {b/1024:8.0f}K {a/1024:8.0f}K  {100*(b-a)/b:4.0f}%")
    print("-" * 68)
    print(
        f"{'TOTAL image bytes':34s} {'':4s} {tot_b/1048576:7.1f}M {tot_a/1048576:7.1f}M  "
        f"{100*(tot_b-tot_a)/max(1,tot_b):4.0f}%{'  (dry-run)' if args.dry_run else ''}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
