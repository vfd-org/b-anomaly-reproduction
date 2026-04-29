# Paper sources

The compiled PDF is shipped at [`main.pdf`](main.pdf) — open it directly to read.

## Layout

```
paper/
├── main.pdf            # camera-ready PDF (25 pages, ~315 KB)
├── main.tex            # LaTeX entry point
├── references.bib      # bibliography
├── sections/           # 10 section files inputted by main.tex
└── figures/            # F1, F2, F3 (PDF + PNG)
```

## Recompile from source

The paper uses `natbib` + `enumitem` + `graphicx` + `booktabs` + `hyperref` —
all standard packages on any modern LaTeX install.

### Option A — Tectonic (recommended; no sudo, single static binary)

```bash
# install once (~50 MB)
curl -L https://github.com/tectonic-typesetting/tectonic/releases/download/tectonic%400.15.0/tectonic-0.15.0-x86_64-unknown-linux-musl.tar.gz \
  | tar -xz -C ~/.local/bin/

# compile (downloads CTAN packages on first run, caches afterwards)
~/.local/bin/tectonic -X compile main.tex
```

Tectonic handles bibtex automatically; one invocation produces `main.pdf`.

### Option B — TeX Live (`pdflatex` + `bibtex`)

```bash
# Debian/Ubuntu:
sudo apt-get install texlive-latex-recommended texlive-latex-extra \
                     texlive-bibtex-extra texlive-fonts-recommended

# compile (the bibtex pass and one re-run are needed for cross-references)
pdflatex main.tex
bibtex main
pdflatex main.tex
pdflatex main.tex
```

### Option C — Overleaf

Upload the entire `paper/` directory as a project. Overleaf's default
`pdflatex` toolchain will work; set the main document to `main.tex`.

## Regenerating figures

The three figures (F1 kernel shape, F2 bin pulls, F3 cross-dataset
amplitudes) are produced by `../scripts/wo017_paper_figures.py`. They
are regenerated whenever you run `bash ../repro/run_all.sh`.

To regenerate just the figures:

```bash
cd ..
PYTHONPATH=src python3 scripts/wo017_paper_figures.py
```

Output goes into `paper/figures/` as both `.pdf` (vector, used by LaTeX)
and `.png` (preview, used by README.md at the repository root).

## Review history

The paper went through three rounds of internal hostile review. The
Round 2 review found that the linearised Mode-B fit's
$\Delta\mathrm{AIC}=-1.67$ on LHCb 2025 did not survive a non-linear
`flavio.np_prediction` refit (it flipped to $+1.09$). The paper was
rewritten around the negative finding and accepted as preprint-ready
in Round 3.
