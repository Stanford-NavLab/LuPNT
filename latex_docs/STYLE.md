# Chapter style contract

Every file under `chapters/` is `\include`d by `lupnt.tex` after `preamble.tex`
has been read. A chapter is therefore a *fragment*: it must not contain
`\documentclass`, `\usepackage`, `\begin{document}`, or a bibliography.

## Skeleton

```latex
\chapter{Time Systems and Conversions}
\label{chap:time-conversions}

Opening paragraph: what this chapter specifies, which source files implement
it, and which references it follows.

\section{...}
\label{sec:tc-...}
```

## Macros available (defined in `preamble.tex`)

| Macro | Use for |
|---|---|
| `\code{x}` | identifier, enum value, CLI flag, short expression |
| `\cpp{X}` / `\py{x}` | a C++ / Python symbol |
| `\file{cpp/lupnt/x.cc}` | a repository path (breaks at `/` and `_`) |
| `\lupnt{}` / `\pylupnt{}` | the library name / the Python package |
| `\vec{r}` `\mat{P}` `\uvec{e}` | 3-vector, matrix, unit vector |
| `\dd` `\transpose` | upright d, transpose superscript |
| `\diag \atantwo \argmin \Tr \erf` | operators |
| `\TT \TAI \TDB \TCB \TCG \TCL \LT \UTC \UTone \GPST` | time scales in math mode |

Units go through `siunitx`: `\SI{32.184}{\second}`, `\si{\meter\per\second}`,
`\num{1.5e-8}`. Extra units declared in the preamble: `\yr`, `\dayunit`,
`\bit`, `\dB`, `\dBHz`, `\dBW`, `\dBi`, `\TECU`, `\au`.

## Environments

- `\begin{lupntnote} ... \end{lupntnote}` — replaces `.. note::`
- `\begin{lupntwarning} ... \end{lupntwarning}` — replaces `.. warning::`
- `\begin{lstlisting}[language=cpp, caption={...}, label={lst:...}]` —
  languages: `cpp`, `py`, `yamlcfg`, `shell`, `bibtexlang`
- Tables: `booktabs` (`\toprule/\midrule/\bottomrule`); wide tables use
  `tabularx` with the ragged-right `L` column type.
- TikZ styles ready to use: `tsnode`, `boxnode`, and the edge classes
  `constant`, `table`, `linear`, `model`, `integral`, `flow`.

## Figures

`latex_docs/` is **text only** — never add an image or a built PDF here.

- **Draw it in TikZ** in the chapter that uses it. This is the default: it adds
  no bytes, stays in sync with the prose, and needs no build step.
- **A raster artefact already committed in the repo** (`assets/`, or
  `docs/_static/`) is on the `\graphicspath`; use plain
  `\includegraphics{name.png}` with no directory prefix. Note LaTeX cannot read
  SVG: an SVG must be rendered to PDF into the data directory first, as
  `build.sh` does for the title-page logo.
- **A generated or large artefact** lives in `$LUPNT_DATA_PATH/docs/` and is
  included with `\datafig{name.pdf}`, which renders a visible placeholder
  instead of failing when the data directory has not been fetched. Never
  `\includegraphics` such a file directly — that turns a missing artefact into
  a broken build.

## Labels and cross references

Prefix by kind and scope: `chap:<slug>`, `sec:<slug>-<name>`,
`eq:<slug>-<name>`, `tab:<slug>-<name>`, `fig:<slug>-<name>`,
`lst:<slug>-<name>`. Reference with `\cref` / `\Cref`, never with a hard-coded
number.

## Citations

`\citep{key}` / `\cite{key}` with keys from `refs.bib` only. Add a new entry to
`refs.bib` rather than inventing a key. The PhD thesis is `iiyama2026thesis`;
chapter/equation numbers in it are given as
`\citep[Ch.~6.2]{iiyama2026thesis}`.

## Code listings

Every listing names the file it came from on its first line:

```latex
\begin{lstlisting}[language=cpp]
// cpp/lupnt/conversions/time_conversions.cc :: TaiToTt
Real TaiToTt(Real t_tai) { return t_tai + TT_TAI_OFFSET; }
\end{lstlisting}
```

Listings must be **ASCII only** — no `→`, `↔`, `μ`, `≈`, en dashes. Use the
math equivalents in surrounding prose instead. The same applies to the rest of
the chapter: non-ASCII characters belong in math mode (`$\mu$`, `$\approx$`),
never bare in the text.

## Depth

The chapters are specifications, not summaries. Each model gets:

1. the state/unit/time-scale contract it assumes,
2. the governing equations, numbered, with every symbol defined,
3. the code that implements them, with the source path,
4. the configuration knobs that change the behaviour, and
5. an explicit *model boundaries* section: what is approximated, what is
   omitted, and what accuracy is expected.
