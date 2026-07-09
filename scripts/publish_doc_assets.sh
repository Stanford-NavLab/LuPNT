#!/usr/bin/env bash
#
# Publish the interactive documentation assets (Plotly HTML + Cesium scenes) staged under
# build/doc_assets/ to the lupnt-doc-assets GitHub Pages repo, so the <iframe> URLs embedded
# in the tutorial notebooks resolve. These files are intentionally NOT committed to this repo
# (they are large, regenerated interactive HTML) — they live in a separate Pages repo.
#
# Prerequisites:
#   1. `pixi run render-notebooks`   (writes build/doc_assets/{plots,cesium}/*.html)
#   2. push access to the assets repo
#   3. one-time: enable Pages on the assets repo (Settings -> Pages -> Deploy from branch:
#      main / root). The base URL must match _doc_assets.py's ASSETS_BASE
#      (https://stanford-navlab.github.io/lupnt-doc-assets by default).
#
# Usage:
#   scripts/publish_doc_assets.sh [assets_repo_url]
#     default url: git@github.com:Stanford-NavLab/lupnt-doc-assets.git
set -euo pipefail

ASSETS_REPO="${1:-git@github.com:Stanford-NavLab/lupnt-doc-assets.git}"
SRC="$(git rev-parse --show-toplevel)/build/doc_assets"

[ -d "$SRC" ] || { echo "ERROR: $SRC not found — run 'pixi run render-notebooks' first." >&2; exit 1; }

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
echo "Cloning $ASSETS_REPO ..."
git clone --depth 1 "$ASSETS_REPO" "$WORK"

mkdir -p "$WORK/plots" "$WORK/cesium"
cp -f "$SRC"/plots/*.html   "$WORK/plots/"   2>/dev/null || true
cp -f "$SRC"/cesium/*.html  "$WORK/cesium/"  2>/dev/null || true
touch "$WORK/.nojekyll"     # let Pages serve the files verbatim

cd "$WORK"
git add -A
if git diff --cached --quiet; then
  echo "Assets already up to date — nothing to push."
  exit 0
fi
git commit -m "Update interactive documentation assets"
git push
echo
echo "Done. If this is the first push, enable Pages on the assets repo:"
echo "  Settings -> Pages -> Build and deployment -> Deploy from a branch -> main / (root)"
