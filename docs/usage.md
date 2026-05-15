# HAlign-4 Usage (CLI)

This document explains the command-line arguments of `halign4` and provides runnable examples using the datasets under `test/data/`.

> Notes
>
> - `halign4` needs **two mandatory arguments**: `-i/--input` and `-o/--output`.
> - `-w/--workdir` is optional. If not provided, the program will create a default workdir: `./tmp-<random>`.
> - Some options are validated by CLI11 at parse time (for example `-i` requires an existing file).
> - `-p/--msa-cmd` supports **keywords** (`minipoa` / `mafft` / `clustalo`) and also supports a **custom command template string**.
>   It is **not** required to be a file path.

---

## Quick start

Minimal run (uses built-in defaults):

```bash
./build/halign4 \
  -i test/data/mt1x.fasta.gz \
  -o out.fasta
```

---

## Arguments reference

### Required

#### `-i, --input <path>`
Input sequence file path.

- Expected format: FASTA (optionally gzipped, depending on your build)
- Validation: must exist (`CLI::ExistingFile`)

#### `-o, --output <path>`
Output file path. The main result is written here.

- Validation: required (but does not need to exist)

---

### Optional

#### `-v, --version`
Print version information and exit.

#### `-w, --workdir <path>`
Working directory used for intermediate files.

- The program will create sub-directories under this path (for example `data/`, `temp/`, `result/`).
- **Important**: some builds may require `workdir` to be empty to avoid overwriting previous outputs (see project README).

#### `-t, --thread <int>`
Number of CPU threads.

- Default: hardware concurrency
- Range check: `1 .. 100000`

#### `--kmer-size <int>`
K-mer size used by minimizer / seeding logic.

- Default: `19`
- Range check: `4 .. 31`

#### `--kmer-window <int>`
Minimizer window size **in number of k-mers**.

- Default: `19`

#### `--cons-n <int>`
Number of sequences selected for consensus step (Top-K by sequence length).

- Default: `1000`

#### `--sketch-size <int>`
Sketch size used by Mash/MinHash-related components.

- Default: `3000`

#### `--profile-ref-kmer-len <int>`
K-mer size used when building mash sketches for profile reference search.

- Default: `10`
- Range check: `4 .. 31`
- This is separate from `--kmer-size`, which still controls minimizer/seeding.

#### `-r, --ref <path>`
Provide an explicit center/reference sequence file (FASTA).

- If provided, the program will use these sequences as the reference/center set instead of auto-selecting.
- Validation: must exist (`CLI::ExistingFile`)

#### `-a, --align-ref <path>`
Provide the **pre-aligned MSA** corresponding to `-r/--ref`.

- This option is useful when your reference set has already been aligned and you do **not** want HAlign-4 to run the external MSA step again.
- Validation: must exist (`CLI::ExistingFile`)
- Relationship: this flag is intended to be used together with `-r/--ref`.

Example:

```bash
./build/halign4 \
  -i test/data/covid-test.fasta.gz \
  -o covid.out.fasta \
  -w covid.work \
  -r test/data/covid-ref.fasta.gz \
  --ref-align test/data/covid-ref.aligned.fasta
```

#### `-s, --score <path>`
Provide a DNA5 scoring matrix file for reference alignment.

- Internal base order is `A C G T N` for both rows and columns.
- The file may contain exactly 25 numeric scores, or a header/row labels like `score.example.tsv`.
- Scores must fit signed int8 range `[-128, 127]`.
- Comments start with `#`; commas are accepted as separators.

Example:

```bash
./build/halign4 \
  -i input.fasta \
  -o aligned.fasta \
  -s score.example.tsv \
  --gap-open 10 \
  --gap-extend 2
```

#### `--gap-open <int>`, `--gap-extend <int>`
Set affine gap penalties used by reference alignment.

- `--gap-open` default: `10`
- `--gap-extend` default: `2`
- Both values must be in `[0, 127]`.

#### Profile reference selection: `--profile-ref-*`
When `-r/--ref` contains multiple reference sequences, `seq2profile` ranks candidate references by HAlign's mash sketch search. Sketch hashes are sorted and deduplicated, and hashes shared by all references may have already been removed. HAlign converts the retained-hash Jaccard score to Mash ANI before applying the sequence-similarity threshold.

- `--profile-ref-min <int>`: keep at least this many references after sorting by Mash ANI. Default: `15`.
- `--profile-ref-max <int>`: keep at most this many references. Default: `40`.
- `--profile-ref-min-similarity <float>`: after `--profile-ref-min`, keep additional references only when Mash ANI is at least this value. Default: `0.7`.

#### `-p, --msa-cmd <string>`
MSA command **keyword** or **command template string**.

Supported keywords:

- `minipoa` (default)
- `mafft`
- `clustalo`

If you provide a custom template string, it may contain placeholders:

- `{input}`: required, replaced with the input FASTA path for the external MSA
- `{output}`: required, replaced with the output FASTA path for the external MSA
- `{thread}`: optional, replaced with the integer from `-t/--thread`

Built-in default template (used when you don’t pass `-p`):

```text
minipoa {input} -S -t {thread} -r1 > {output}
```

Built-in MAFFT template:

```text
mafft --thread {thread} --auto {input} > {output}
```

Built-in Clustal Omega template:

```text
clustalo -i {input} -o {output} --threads {thread}
```

Important note:

- In `halign4`, **minipoa is the default high-quality aligner** (via the built-in template above).
- In the examples below we use **MAFFT** only because it’s a common MSA tool and its CLI is easy to demonstrate.

Security note:

- The command is executed via the system shell. Treat template inputs as trusted data.

#### `--keep-length`
Keep reference sequences in `-r/--ref` ungapped in the final MSA.

- Source-level meaning (matches `RefAligner` implementation):
  - When set, the pipeline removes alignment columns that would introduce gaps into the reference sequences.
- When `-r` contains multiple reference sequences:
  - **all** reference sequences are guaranteed to have no inserted gaps.

> Important clarification
>
> This flag is about **reference sequences from `-r`** (the "center/reference FASTA"), not about general query sequences.


#### `--save-workdir`
Keep the working directory after successful completion.

- Default behavior: remove workdir on success

---

## Examples

### Example 1: Minimal dataset (`mt1x`) + demonstrate `-p/--msa-cmd` (using MAFFT keyword)

Dataset:

- `test/data/mt1x.fasta.gz`

Goal:

- Show the most basic run.
- Show how to use `-p` with the new keyword mapping.

Run:

```bash
./build/halign4 \
  -i test/data/mt1x.fasta.gz \
  -o mt1x.out.fasta \
   -w mt1x.work \
   -t 8 \
  -p mafft
```

If you don’t have `mafft` installed, either install it or switch the template file back to a command that exists in your environment.

---

### Example 2: COVID dataset + demonstrate `-r`, `--keep-length`

Dataset:

- `test/data/covid-ref.fasta.gz`
- `test/data/covid-test.fasta.gz`

Background:

- The first record in `covid-ref.fasta.gz` is the Wuhan reference sequence.
- Other sequences are consensus sequences for variants produced by:
  https://github.com/corneliusroemer/pango-sequences

#### 2.1 Keep all reference sequences ungapped (`--keep-length`)

```bash
./build/halign4 \
  -i test/data/covid-test.fasta.gz \
  -o covid.out.fasta \
  -w covid.work \
  -r test/data/covid-ref.fasta.gz \
  --keep-length
```

---

## How to understand `--keep-length` behavior (toy example)

This section uses a tiny, **made-up** example to make the idea concrete. All alignment rows below have the **same length**.

Assume your `-c/--center-path` has **two** reference sequences:

```text
ref1 (first):  ACGTAC
ref2:          ACGTC
```

And you have one query sequence:

```text
q1:            ACGTTAC
```

A generic MSA tool may produce an alignment like this (gaps inserted into references are allowed):

```text
ref1  ACGT-AC
ref2  ACGT--C
q1    ACGTTAC
```

### Case A: default behavior (no `--keep-...` flags)

- Gaps in reference sequences are allowed.
- This can change the *effective* coordinate system of references.

### Case B: `--keep-length`

- Guarantee: **all references in `-r` will not contain inserted gaps**.
- The pipeline will drop columns that would introduce gaps in any reference.

Using the original MSA snippet, dropping the column where `ref1` has `-` *and* the column where `ref2` has `-` yields:

```text
ref1  ACGTAC
ref2  ACGT-C
q1    ACGTAC
```

What to take away:

- Using `--keep-length` will reduce the number of alignment columns (because it removes "reference-gap columns").
- This ensures all reference sequences maintain their original coordinate system without insertions.

> This toy example is meant to build intuition; exact output depends on the real sequences and the chosen MSA method.

---

## Troubleshooting

- **`-p/--msa-cmd` fails at startup even though the tool exists**: `halign4` runs a tiny self-check during argument validation.
  If it fails, try running the expanded command manually to see stderr, or use a custom template.
  On Windows/WSL setups, make sure `halign4` and the external MSA tool are in the **same environment** (both in WSL or both native).
- **Workdir already exists**: remove it, choose a new `-w`, or build/run in Debug mode if your build allows reusing a non-empty workdir.
