# HAlign-4

[![Downloads](https://anaconda.org/malab/halign4/badges/downloads.svg)](https://anaconda.org/malab/halign4)
[![License](https://anaconda.org/malab/halign4/badges/license.svg)](https://anaconda.org/malab/halign4)
[![Platforms](https://anaconda.org/malab/halign4/badges/platforms.svg)](https://anaconda.org/malab/halign4)

[HAlign 4: A New Strategy for Rapidly Aligning Millions of Sequences.](https://doi.org/10.1093/bioinformatics/btae718)

Documentation:

- Detailed usage & examples: [`docs/usage.md`](docs/usage.md)
- Source install & dependencies: [`docs/install.md`](docs/install.md)
- Tests: [`docs/test.md`](docs/test.md)

---

## Install (Conda)

Conda is the recommended installation method for end users.

```bash
conda install -c malab halign4
```

Verify:

```bash
halign4 --version
halign4 -h
```

Source installation: see [`docs/install.md`](docs/install.md).

---

## Quick start

The repository includes small datasets under `test/data/` which are perfect for a first run.

Minimal example:

```bash
halign4 \
  -i test/data/mt1x.fasta.gz \
  -o mt1x.out.fasta
```

---

## Parameters (overview)

The most important parameters are:

- `-i/--input`: input FASTA (required)
- `-o/--output`: output aligned FASTA (required)
- `-w/--workdir`: working directory (optional; default: `./tmp-<random>`)
- `-p/--msa-cmd`: MSA method (keyword: `minipoa`/`mafft`/`clustalo`, or a custom template)
- `-r/--ref`: provide a reference/center FASTA (optional)
- `--align-ref`: provide the pre-aligned MSA for the `-r/--ref` reference set
- `-s/--score`: provide a DNA5 scoring matrix file (see `score.example.tsv`)
- `--gap-open`, `--gap-extend`: control gap open/extension penalties
- `--profile-ref-min`, `--profile-ref-max`, `--profile-ref-min-similarity`: mash-based reference profile selection controls
- `--keep-length`: keep reference length coordinate rules

For the full parameter list and detailed examples, see [`docs/usage.md`](docs/usage.md).

---

## Tests

See [`docs/test.md`](docs/test.md) for how to run tests under the `test/` directory.

---

## Citation

If you use HAlign-4 in academic work, please cite:

HAlign 4: a new strategy for rapidly aligning millions of sequences. Bioinformatics, 2024, 40(12): btae718. https://doi.org/10.1093/bioinformatics/btae718
