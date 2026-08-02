"""Convert a GRAIL/PDS spherical-harmonic gravity file (.sha / *_sha.tab) to LuPNT's
`.cof` format, optionally truncating the degree/order.

LuPNT's `ReadHarmonicGravityField` (cpp/lupnt/environment/body.cc) parses a fixed-column
`.cof` layout:

    POTFIELD<n:3><m:3>  0 <GM> <R> <norm>
    RECOEF  <n:3><m:3>   <Cnm> <Snm>

with `n`/`m` read from byte offsets 8 and 11, so the column widths matter. Coefficients are
fully normalised in both formats and transfer across unchanged (the parser de-normalises
when called with `normalized=true`, which `Body::Moon` does).

Truncation matters: a 1200x1200 field is ~721k coefficients (~83 MB); the shipped
`grgm900c.cof` is itself truncated to 360x360 (~65k, ~4 MB), which is the convention
followed here.

Usage:
    pixi run python scripts/convert_sha_gravity.py IN.sha OUT.cof --degree 360 \
        [--model-name GRGM1200B]
"""

import argparse
from pathlib import Path


def convert(src: Path, dst: Path, degree: int, model_name: str) -> None:
    with src.open() as f:
        gm_str, r_str, n_max_str, m_max_str = f.readline().split()
        gm, radius = float(gm_str), float(r_str)
        n_max_file, m_max_file = int(n_max_str), int(m_max_str)

        if degree > min(n_max_file, m_max_file):
            raise SystemExit(f"requested degree {degree} exceeds file's {n_max_file}x{m_max_file}")

        rows, skipped = [], 0
        for line in f:
            parts = line.split()
            if len(parts) < 4:
                continue
            n, m = int(parts[0]), int(parts[1])
            if n > degree:
                skipped += 1
                continue  # not `break`: do not assume the file is sorted by degree
            if m > degree:
                skipped += 1
                continue
            rows.append((n, m, float(parts[2]), float(parts[3])))

    rows.sort(key=lambda r: (r[0], r[1]))  # LuPNT's parser stops at the first n > requested

    header = [
        "COMMENT   5",
        f"C {model_name} lunar gravity model, converted from {src.name}.",
        "C Converted by scripts/convert_sha_gravity.py for LuPNT's .cof reader.",
        f"C Source model is {n_max_file}x{m_max_file}; truncated here to {degree}x{degree}",
        "C to bound memory -- matching the shipped grgm900c.cof, also truncated to 360x360.",
        f"C GM and reference radius are taken from the source file header.",
    ]

    with dst.open("w") as f:
        for line in header:
            f.write(line + "\n")
        # Fixed columns: 'POTFIELD' occupies 0-7, degree 8-10, order 11-13, then the values.
        f.write(f"POTFIELD{degree:3d}{degree:3d}  0 {gm:.14e} {radius:.14e} {1.0:.14e}\n")
        for n, m, c, s in rows:
            # 'RECOEF' + 2 spaces occupies 0-7, degree 8-10, order 11-13, then the values.
            f.write(f"RECOEF  {n:3d}{m:3d}   {c:.14e} {s:.14e}\n")

    print(f"{src.name} ({n_max_file}x{m_max_file}) -> {dst.name} ({degree}x{degree})")
    print(f"  GM     : {gm:.14e} m^3/s^2")
    print(f"  R      : {radius:.14e} m")
    print(f"  kept   : {len(rows)} coefficients (dropped {skipped})")
    print(f"  size   : {src.stat().st_size / 1e6:.1f} MB -> {dst.stat().st_size / 1e6:.1f} MB")


if __name__ == "__main__":
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("src", type=Path)
    p.add_argument("dst", type=Path)
    p.add_argument(
        "--degree", type=int, default=360, help="max degree/order to keep (default: 360)"
    )
    p.add_argument("--model-name", default="", help="model name for the comment header")
    a = p.parse_args()
    convert(a.src, a.dst, a.degree, a.model_name or a.src.stem)
