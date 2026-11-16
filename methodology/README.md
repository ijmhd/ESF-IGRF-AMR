> On the performance branch, these PDFs/TeX and the recorded document build
> belong to frozen v0.1.0. For exact source-line matching or rebuilding those
> historical documents, check out tag `v0.1.0`. Current execution changes and
> code links are in `../docs/PERFORMANCE_UPDATE.md` and `../docs/FUNCTIONS.md`.

# Building the methods documents

`main.tex` is the primary methods manuscript. `source-atlas.tex` reproduces the
complete annotated runtime sources, headers and C tests with line numbers.
The included PDFs correspond to the frozen source release.

From this directory, with TeX Live/MacTeX and latexmk:

```sh
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
latexmk -pdf -interaction=nonstopmode -halt-on-error source-atlas.tex
```

Use the whole repository: the atlas includes `../src`, `../include`, `../tests`
and `../CMakeLists.txt`. `references.bib` records literature and formula sources.
The parameter and function appendices are checked against the actual CLI/source
inventory. Run evidence is summarized in `../validation/qualification.json`.

This is a methods manuscript and teaching reference, not a statement of peer
review. Authorship/event-calibration and publication venue decisions remain with
the scientific contributors. Do not describe implementation-parity tests as
independent observational validation.
