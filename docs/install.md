# Install HAlign-4

This document explains how to install **HAlign-4 / halign4**.

- If you want a quick install, prefer **Conda** (see the project `README.md`).
- If you want to build from source (development / customizing), follow the steps below.

> Environment note
>
> `halign4` may call external MSA tools (`minipoa`, `mafft`, `clustalo`) during argument validation (template self-check).
> Please make sure `halign4` and those external tools are installed in the **same environment** (e.g. all in WSL, or all native Linux).

---

## 1. Build from source (Linux / WSL)

### 1.1 Prerequisites

Required build tools:

- CMake >= 3.28
- A C++20 compiler (GCC/Clang)
- Make or Ninja
- OpenMP support (usually via `libgomp` / `libomp`)

Optional but recommended:

- `zlib` (for reading `.gz` FASTA)
- `libcurl` (for URL inputs)

Example (Ubuntu / Debian):

```bash
sudo apt-get update
sudo apt-get install -y \
  build-essential \
  cmake \
  ninja-build \
  zlib1g-dev \
  libcurl4-openssl-dev
```

> If your distro ships an older CMake, consider using `conda` or installing a newer CMake.

### 1.2 Clone and build

```bash
git clone https://github.com/pinglu-zhang/HAlign-4.git
cd HAlign-4

cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

The binary will be:

- `build/halign4`

### 1.3 Quick sanity check

```bash
./build/halign4 --version
./build/halign4 -h
```

---

## 2. Install / prepare external MSA tools

HAlign-4 uses an external MSA method to align consensus and insertion sequences.

### 2.1 MAFFT

Install (Conda recommended):

```bash
conda install -c conda-forge -c bioconda mafft
```

Verify:

```bash
mafft --version
```

### 2.2 Clustal Omega (clustalo)

Install (Conda recommended):

```bash
conda install -c conda-forge -c bioconda clustalo
```

Verify:

```bash
clustalo --version
```

### 2.3 minipoa

Project default MSA tool is `minipoa`.

- Upstream repo: https://github.com/NCl3-lhd/minipoa

Install options:

#### Option A: Conda (recommended)

If available in your channels:

```bash
conda install -c malab -c conda-forge -c bioconda minipoa
```

Verify:

```bash
minipoa -h
```

#### Option B: Build minipoa from source

Follow the instructions in the minipoa repository. A typical pattern is:

```bash
git clone https://github.com/NCl3-lhd/minipoa.git
cd minipoa
make -j
# ensure the resulting binary is on PATH, e.g.
# export PATH="$PWD:$PATH"
```

Verify:

```bash
minipoa -h
```

---

## 3. Common issues

### 3.1 `msa-tool template test failed` at startup

During argument validation, `halign4` runs a tiny self-check command via the system shell.

Checklist:

1. Verify the tool is available on PATH **in the same environment** where `halign4` runs:

```bash
command -v mafft
command -v clustalo
command -v minipoa
```

2. Verify the tool runs without interactive prompts:

```bash
mafft --version
clustalo --version
minipoa -h
```

3. If you use `--msa-tool` with a custom template string, ensure it contains `{input}` and `{output}`.

---

## 4. Next: usage and examples

- CLI usage & examples: [`docs/usage.md`](usage.md)
- Tests: [`docs/test.md`](test.md)

