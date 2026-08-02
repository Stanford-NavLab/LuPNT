# The LuPNT Manual (LaTeX)

A single self-contained PDF: the design and mathematical specifications of
LuPNT, with full derivations, unit and time-scale contracts, and the
implementing C++/Python code inline.

This tree lives at the repository root, deliberately **outside** `docs/`: the
Sphinx site and this manual serve different audiences, and the manual is
intended to stand on its own. Nothing here depends on `docs/` building, or on
`docs/` existing at all.

## Build

No system LaTeX installation is required — `pixi exec` fetches [tectonic] into
a throwaway environment, so nothing is added to `pixi.toml` or the lockfile.

```bash
cd latex_docs
./build.sh              # -> build/lupnt.pdf
./build.sh --keep       # keep .aux/.bbl for debugging
```

`build.sh` uses `tectonic` directly if it is on `PATH`, otherwise `pixi exec`.
Tectonic runs BibTeX itself, so `refs.bib` needs no separate pass — but it does
need to sit beside `lupnt.tex`, which is why the bibliography is committed here
rather than generated.

With a traditional TeX installation the equivalent is:

```bash
pdflatex lupnt && bibtex lupnt && pdflatex lupnt && pdflatex lupnt
```

Two passes after BibTeX are required: `cleveref` resolves cross-references on
the second.

## Layout

| Path | Contents |
|---|---|
| `lupnt.tex` | main document — title page, front matter, `\part`/`\include` structure |
| `preamble.tex` | every package, listing style, TikZ style, and macro the chapters may use |
| `chapters/*.tex` | one file per chapter, `\include`d in the order set by `lupnt.tex` |
| `refs.bib` | the bibliography; keys are stable and cited by name |
| `STYLE.md` | the contract every chapter is written against |
| `build.sh` | the build entry point |

**This directory is text only** — `.tex`, `.bib`, `.sh`, `.md`. No images, no
built PDFs; `.gitignore` enforces it. Keeping binaries out is what stops the
library from growing every time a figure is regenerated. See *Figures* below
for where they live instead.

Chapters are *fragments*: no `\documentclass`, no `\usepackage`, no
`\begin{document}`. Everything a chapter needs is in `preamble.tex`, which is
what makes `\includeonly{chapters/filters}` work for a fast single-chapter
rebuild while you edit.

## Adding a chapter

1. Read `STYLE.md`. It fixes the macros (`\code`, `\file`, `\cpp`, `\lupnt`,
   the time-scale and vector macros), the listing languages, the admonition
   environments, and the label prefixes.
2. Write `chapters/<name>.tex`, starting with `\chapter{...}` and
   `\label{chap:<slug>}`.
3. Add `\include{chapters/<name>}` to `lupnt.tex` under the right `\part`.
4. Add any new bibliography entries to `refs.bib` — do not invent a key in the
   prose and leave the entry unwritten; an undefined citation compiles to a
   silent `[?]`.

## Figures

Nothing binary is committed here. Figures come from two places, both resolved
by `\graphicspath` in `preamble.tex`:

| Source | For | Example |
|---|---|---|
| `../assets/`, `../docs/_static/` | brand assets already in git | `lupnt_logo_horizontal.svg` |
| `$LUPNT_DATA_PATH/docs/` | generated or large artefacts, fetched not committed | `tutorials/*.png` |

The title-page logo is the repository logo, `assets/lupnt_logo_horizontal.svg`.
LaTeX cannot read SVG, so `build.sh` renders it once to
`$LUPNT_DATA_PATH/docs/lupnt_logo_horizontal.pdf` (via `rsvg-convert`, fetched
through `pixi exec` if it is not installed). If that fails the build still
succeeds and the title page shows a placeholder.

`build.sh` resolves `LUPNT_DATA_PATH` — from the environment (pixi activation
sets it), then `.env`, then the default `data/LuPNT_data` — and writes the
absolute path into `figpath.tex`, which `preamble.tex` reads. `figpath.tex` is
generated and gitignored.

Chapters pull data-directory figures with `\datafig{name.pdf}` rather than
`\includegraphics`. **A missing artefact does not break the build**: `\datafig`
falls back to a visible placeholder box, so the manual still compiles from a
bare clone with no data directory fetched.

The tutorial figures are extracted from the example notebooks' stored outputs.
To refresh them:

```bash
pixi run python latex_docs/scripts/extract_notebook_figures.py
```

Note that `data/LuPNT_data` is itself fetched from a public GitHub Release
(`cmake/FetchLuPNTData.cmake`), so putting a figure there makes it available to
*your* builds; distributing it to other clones means republishing that release asset.

Everything else is drawn in TikZ in the chapter that uses it, so it stays in
sync with the prose, adds no bytes, and needs no build step. Prefer TikZ for
anything new.

## Relationship to `docs/`

The chapters were originally derived from `docs/**/*.rst` and then expanded
well past it. This manual is now the primary written reference and is intended
to stand alone. Nothing here reads from
`docs/` except the optional `docs/_static` entry on the `\graphicspath`, which
is a fallback and not required.

What would be lost with `docs/` are the auto-generated API references (Doxygen
+ Breathe for C++, autodoc for Python) and the executed notebooks. Those are
generated artefacts, not prose, and are not duplicated here.

## Provenance

LuPNT is developed at the **Stanford NAV Lab**
(<https://github.com/Stanford-NavLab/LuPNT>).

[tectonic]: https://tectonic-typesetting.github.io/
