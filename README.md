# KoNuSeg
[![GitHub Release](https://img.shields.io/github/v/release/BilkentCompGen/KoNuSeg?color=blue&cacheSeconds=60)](https://github.com/BilkentCompGen/KoNuSeg/releases)
[![Last Commit](https://img.shields.io/github/last-commit/BilkentCompGen/KoNuSeg)](https://github.com/BilkentCompGen/KoNuSeg/commits/main)
[![License](https://img.shields.io/github/license/BilkentCompGen/KoNuSeg?cacheSeconds=60)](https://github.com/BilkentCompGen/KoNuSeg/blob/main/LICENSE)

`KoNuSeg` is a C++ tool for fast, efficient repeat detection using k-mer counting and Bloom filters, backed by `btllib`.

---

## Prerequisites & Requirements

Before building `KoNuSeg`, ensure your system meets the following software and library requirements:

* **Compiler**: GCC / G++ version **11** or higher (must support `C++17` and OpenMP).
* **Dependencies**:
  * [btllib](https://github.com/bcgsc/btllib): Library for bioinformatics data structures (Bloom filters, sequence parsing, etc.). This library is included in this repository (version 1.4.9). 
---

## Installation

### Clone the Repository
```bash
git clone git@github.com:BilkentCompGen/KoNuSeg.git
cd KoNuSeg
```

## Compilation & Verification

The project includes a `Makefile` to simplify compilation and testing. Ensure `g++-11` is available at `/usr/bin/g++-11` on your system, or update the `CC` variable in the `Makefile` to match your local compiler path.

### 1. Build the Executable
To compile the optimized production binary:
```bash
make
```

### 2. Sanity check
Before running `KoNuSeg` on your own dataset, execute the included test suite as a sanity check to verify that the binary was built correctly and all bundled dependencies function properly:
```bash
make run_test
```
### 3. To remove generated binary files:
```bash
make clean
```

## Usage

```bash
./bin/konuseg.out -fq <fastq> -fa <fasta> -o <bed> [options]
```

### Required Parameters

| Parameter | Description |
| :--- | :--- |
| `-fq` | Path to input reads file (FASTQ or FASTA, supports `.gz`/`.bz2`/`.xz`; optional when `-load_cbf` is used) |
| `-fa` | Path to input assembly FASTA file (supports `.gz`/`.bz2`/`.xz`) |
| `-o` | Path to output BED file (merged output written to `<o>_merged.bed`) |

### Optional Parameters

| Option | Default | Description |
| :--- | :--- | :--- |
| `-g` | `3000000000` | Genome size (automatically sets `-anum` and `-knum` when those flags are omitted) |
| `-s` | `10010111001110100101110011101001` | Binary spaced-seed mask (must match `-k`; auto-synced if only `-k` or `-s` is given) |
| `-k` | `32` | K-mer length |
| `-fpr` | `0.01` | Target Bloom filter false-positive rate |
| `-w` | `5` | Sliding window size for median outlier filtering |
| `-minseg` | `500` | Minimum segment length (bp) for Stage 2 merging |
| `-c_thresh` | `1` | Copy-number threshold for Stage 1 segmentation and Stage 2 boundary compatibility |
| `-anum` | `1500000000` | Expected distinct k-mers in assembly/reference |
| `-knum` | `1500000000` | Expected distinct k-mers in reads |
| `-save_bf` | — | Save populated Seed Bloom Filter (and `.meta` counts) to file |
| `-load_bf` | — | Load prebuilt Seed Bloom Filter from file |
| `-save_cbf` | — | Save populated Counting Bloom Filter (and `.meta` counts) to file |
| `-load_cbf` | — | Load prebuilt Counting Bloom Filter from file (skips SBF build and read streaming) |

If `-anum` and `-knum` are not known, passing `-g <genome_size>` automatically uses the genome size to dimension both filters.

> **Note on Counter Saturation**: `KoNuSeg` uses 16-bit saturating counters (`btllib::CountingBloomFilter16`, max raw count `65,535`). At coverage $C$, estimated copy numbers saturate safely at $65,535 / C$ (e.g., ~120 copies at 545× effective Bloom filter coverage), which does not affect repeat classification (`>1`, `>2`, `>5`, `>10`). 

---

### Example Command

```bash
./bin/konuseg.out \
    -fa sample_ref.fa \
    -fq sample_reads.fq \
    -o output_segments.bed \
    -g 100000000 \
    -s 11111111111111111111111111111111 \
    -k 32 \
    -fpr 0.01 \
    -w 5 \
    -minseg 64 \
    -c_thresh 1
```
