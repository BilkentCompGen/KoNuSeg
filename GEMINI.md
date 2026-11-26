# KoNuSeg — Project Knowledgebase

## 1. Overview & Scientific Purpose

**KoNuSeg** is a C++17 / OpenMP bioinformatics tool developed by the Alkan Lab (Bilkent University) and Birol Lab (BC Cancer Research Institute / UBC) for fast, space-efficient, mapping-free detection of repetitive sequences and collapsed segmental duplications in *de novo* genome assemblies.

- **Reference Manuscript**: Zambaku K., Roman-Brenes R., Coombe L., Warren R. L., Wong J., Birol İ., & Alkan C. (2026). *"Space-efficient detection of repeats in genome assemblies using Bloom filters."*
- **Core Idea**: Standard *de novo* repeat annotators inspect only the assembled contigs and therefore cannot distinguish a **collapsed duplication** (assembled as a single copy) from a true single-copy locus. Sequencing reads sample all genomic copies proportionally, so k-mers originating from collapsed repeats appear at multiples of the baseline haploid coverage. Instead of expensive multi-mapping read alignments, KoNuSeg uses a **spaced-seed Bloom filter** (`btllib::SeedBloomFilter`) and a **16-bit Counting Bloom filter** (`btllib::CountingBloomFilter16`) to count error-filtered read k-mers along assembly contigs, smooth hash-collision outliers, and partition contigs into copy-number-annotated BED intervals.
- **Primary Downstream Use Cases**:
  1. **Repeat-aware scaffolding**: Masking segments with estimated copy number exceeding a chosen threshold (`>1`, `>2`, `>5`, or `>10`) via `bedtools maskfasta` before scaffolding draft contigs (e.g., GoldRush contigs scaffolded with `ntLink`) to prevent repeat-induced misjoins.
  2. **Self-alignment repeat profiling**: Supplying an assembly FASTA as both contig and read input (converted to FASTQ) to detect exact/near-exact intragenomic duplications at ~1× baseline depth without whole-genome BLAST alignments.

---

## 2. Repository Structure & External API Policy

> **CRITICAL RULE — DO NOT MODIFY `btllib`**:
> `btllib` (`incl/btllib/`, `incl/cpptoml/`, `incl/tests/`, `lib/libbtllib.a`, and `lib/btllib/`) is a vendored external API (v1.4.9) and **MUST NOT be modified under any circumstances**. All bug fixes, optimizations, and feature changes must be implemented strictly in KoNuSeg's own source and header files (`src/*.cpp`, `incl/*.hpp` outside `incl/btllib/`, and `Makefile`).

```text
KoNuSeg/
├── Makefile                  # Build rules (build, debug, clean, run_test)
├── README.md                 # User installation and CLI reference
├── LICENSE                   # BSD 3-Clause License
├── bin/
│   └── konuseg.out           # Compiled executable target
├── incl/
│   ├── konuseg.hpp           # Umbrella header included by konuseg.cpp
│   ├── data_structs.hpp      # Config, Sequence, SegmentInfo, typedefs, globals, roundCount()
│   ├── loadAssembly.hpp      # FASTA loader declaration
│   ├── populateBF.hpp        # SeedBloomFilter builder declaration
│   ├── loadReads.hpp         # FASTQ streamer & CountingBloomFilter16 builder declarations
│   ├── checkRepeats.hpp      # Contig querying, outlier filtering, Stage 1 segmentation
│   ├── mergeSegments.hpp     # Stage 2 segment pooling & copy-number estimation
│   ├── btllib/               # Bundled btllib v1.4.9 headers (Bloom filters, ntHash, SeqReader, sdsl)
│   ├── cpptoml/              # Header-only TOML parser used by btllib serialization
│   └── tests/                # Bundled upstream btllib test files
├── lib/
│   ├── libbtllib.a           # Prebuilt static btllib v1.4.9 archive
│   └── btllib/               # Upstream btllib Python wrapper artifacts
├── src/
│   ├── konuseg.cpp           # CLI parser, pipeline driver, timing logs, main()
│   ├── loadAssembly.cpp      # Loads contigs from FASTA into in-memory SEQLIST
│   ├── populateBF.cpp        # Populates SeedBloomFilter from assembly k-mers
│   ├── loadReads.cpp         # Streams FASTQ reads and populates CountingBloomFilter16
│   ├── checkRepeats.cpp      # Queries CBF, runs noiseFilter() and createSegments()
│   └── mergeSegments.cpp     # Pools short segments into *_merged.bed via processData()
└── test/
    ├── test.fa               # Small synthetic 5-contig FASTA fixture (seq1: 1x, seq2/seq3: 2x)
    └── test.fq               # Matching 5-read FASTQ fixture for sanity checking
```

---

## 3. Pipeline Architecture & Data Flow

The pipeline executes five sequential phases orchestrated by `main()` in `src/konuseg.cpp`:

### Stage 0: Configuration & Assembly Loading (`src/konuseg.cpp`, `src/loadAssembly.cpp`)
- `parseCommandLine()` populates and validates a `Config` instance (`incl/data_structs.hpp`), automatically synchronizing `-k` and `-s` when only one is explicitly passed.
- `loadAssembly()` streams FASTA contigs (plain or `.gz`/`.bz2`/`.xz` via `btllib::DataSource`) into `SEQLIST` in RAM, extracting clean whitespace-trimmed BED sequence IDs so contigs can be reused during both Bloom filter population and repeat querying.

### Stage 1: Populating the Seed Bloom Filter (`src/populateBF.cpp`)
- Computes optimal filter parameters from expected distinct assembly k-mers (`-anum`, $n$) and false-positive rate (`-fpr`, $p$):
  $$m = -\frac{n \ln p}{(\ln 2)^2}, \qquad \text{bytes} = \left\lceil \frac{m}{8} \right\rceil, \qquad k_{\text{hash}} = \max\!\left(1, \operatorname{round}\!\left(\frac{m}{n} \ln 2\right)\right)$$
- Initializes `btllib::SeedBloomFilter` (`SBF`) with `bfBytes` and the user's spaced seed (`-s`) (or loads a prebuilt filter via `-load_bf` and optionally saves via `-save_bf`).
- Partitions contigs into 64k-start `SeqSlice` chunks (`k - 1` overlap) and runs rolling `btllib::SeedNtHash` (`#pragma omp parallel for schedule(dynamic, 1)`), skipping `N`/non-`ACGT` k-mers, inserting spaced-seed hashes into `seed_bf`, and accumulating valid assembly k-mers in global `k_assembly`.

### Stage 2: Streaming Reads & Populating the Counting Bloom Filter (`src/loadReads.cpp`)
- Computes optimal counter slot count `cbfSlots` ($m$) and `optHashFunctions` from expected distinct read k-mers (`-knum`) and `-fpr`, converting to `cbfBytes = cbfSlots * sizeof(uint16_t)` to initialize `btllib::CountingBloomFilter16` (`CBF`, 16-bit atomic counters saturating at `65,535`, capping estimated copy numbers at $65,535 / \text{coverage}$), or loads a prebuilt filter directly via `-load_cbf` (and optionally saves via `-save_cbf`).
- Streams input reads (FASTQ or FASTA, uncompressed or `.gz`/`.bz2`/`.xz` via `btllib::DataSource`) in bounded multi-read batches (`max(64, 8 * threads)` reads up to 64 MiB).
- Across each batch, OpenMP workers roll `btllib::SeedNtHash` over 64k-start `ReadSlice`s:
  - **Absent in SBF** (or `N`-containing): Discarded with zero heap allocations.
  - **Present in SBF**: Reuses `h.hashes()` directly to insert into `cbf` and increments `k_reads`.
- **Bloom Filter Coverage** is computed and logged as the exact continuous ratio `double(k_reads) / double(k_assembly)`.

### Stage 3: Querying, Outlier Filtering, & Unmerged Segmentation (`src/checkRepeats.cpp`)
- Pre-sizes `segTable` to `assembly->size()` and processes contigs lock-free via `#pragma omp parallel for schedule(dynamic, 1)`:
  - **Querying**: Rolls a single `btllib::SeedNtHash` across each contig, querying `cbf->contains(h.hashes())` and storing non-zero `[pos, pos + k)` intervals in a local `CONTIGLIST`.
  - **Outlier Filtering (`noiseFilter()`)**: Slides a centered window of size `-w` (`windowSize`) out-of-place over the contig's raw k-mer counts (`std::vector<uint>`) to replace hash-collision spikes with the local window median.
  - **Unmerged Segmentation (`createSegments()`)**: Normalizes k-mer counts by exact continuous `coverage_depth = float(k_reads) / float(k_assembly)`, extends segments while consecutive k-mers remain within `±count_threshold` (`-c_thresh`) of the running mean copy number (`wMedian / kCount`), and stores each contig's `SEGMENTINFO*` into `segTable->at(idx)`.
- Writes unmerged 0-based half-open `[start, end)` intervals sequentially in **deterministic FASTA contig order** to `<outfile_bed>`.

### Stage 4: Segment Pooling & Merging (`src/mergeSegments.cpp`)
- `mergeSegments()` iterates sequentially through `segTable` to smooth short segments into longer genomic blocks:
  - Concatenates consecutive unmerged segments on the same contig while accumulated segment length `< min_seg_len` (`-minseg`), inter-segment gap `<= 3 * k_len`, and without absorbing short flanking segments into already full-sized (`>= min_seg_len`) blocks whose copy state differs by `> count_threshold`.
  - `processData()` sorts pooled copy numbers, removes the top quartile (`(size + 3) / 4` when `size >= 2`, ensuring at least 1 collision outlier is dropped even in 2–3 segment pools), takes the median of the remaining values through `roundCount()`, and writes merged intervals to `<outfile_bed>_merged.bed`.

### Copy-Number Rounding Heuristic (`roundCount()` in `incl/data_structs.hpp`)
- `< 0.25` -> `0.0`
- `[0.25, 0.65)` -> `0.5`
- `>= 0.65`: rounds down to `floor(x)` if fractional part `<= 0.64`, otherwise rounds up to `floor(x) + 1.0`.
- Segments with representative haploid copy state `<= 1` are interpreted as single-copy/neutral; segments `> 1` are classified as duplicated.

---

## 4. Command-Line Interface & Parameter Guidelines

```bash
./bin/konuseg.out -fq <reads.fastq> -fa <contigs.fasta> -o <output.bed> [options]
```

| Flag | Required | Default | Internal Field | Description |
| :--- | :---: | :--- | :--- | :--- |
| `-fq` | Yes* | — | `Config::infile_fastq` | Input sequencing reads (FASTQ/FASTA; supports `.gz`/`.bz2`/`.xz`; *optional with `-load_cbf`) |
| `-fa` | Yes | — | `Config::infile_fasta` | Input draft assembly contigs in FASTA format (supports `.gz`/`.bz2`/`.xz`) |
| `-o` | Yes | — | `Config::outfile_bed` | Output path for unmerged BED (merged written to `<o>_merged.bed`) |
| `-g` | No | `3000000000` | `Config::genome_size` | Estimated genome size (sets `-anum` and `-knum` when those flags are omitted) |
| `-s` | No | `10010111001110100101110011101001` | `Config::seed1` | Binary spaced-seed mask (`1` = match, `0` = wildcard; auto-synced if only `-k` or `-s` is passed) |
| `-k` | No | `32` | `Config::k_len` | K-mer length |
| `-fpr` | No | `0.01` | `Config::fpr` | Target Bloom filter false-positive rate |
| `-w` | No | `5` | `Config::window_size` | Sliding window size for median outlier filtering |
| `-minseg` | No | `500` | `Config::seg_len` | Minimum segment length (bp) target during Stage 2 merging |
| `-c_thresh` | No | `1` | `Config::count_threshold` | Copy-number threshold for Stage 1 segmentation & Stage 2 boundary compatibility |
| `-anum` | No | `1500000000` | `Config::a_num` | Expected distinct k-mers in the assembly (`SeedBloomFilter` sizing) |
| `-knum` | No | `1500000000` | `Config::k_num` | Expected distinct k-mers in the reads (`CountingBloomFilter16` sizing) |
| `-save_bf` | No | — | `Config::save_bf_path` | Save populated `SeedBloomFilter` (plus `.meta`) to disk |
| `-load_bf` | No | — | `Config::load_bf_path` | Load prebuilt `SeedBloomFilter` from disk |
| `-save_cbf` | No | — | `Config::save_cbf_path` | Save populated `CountingBloomFilter16` (plus `.meta`) to disk |
| `-load_cbf` | No | — | `Config::load_cbf_path` | Load prebuilt `CountingBloomFilter16` from disk (skips SBF build & read streaming) |

### Practical Parameter Guidance (from Manuscript Benchmarks)
- **Compact / moderately repetitive genomes** (*C. elegans* ~100 Mbp, *D. melanogaster* ~140 Mbp): Permissive-to-moderate masking thresholds (masking unmerged segments with copy estimate `>1`, `>2`, or `>5`) substantially improve scaffolding NGA50 (up to 4.72× in *C. elegans* and 2.04× in *D. melanogaster*) while reducing misassemblies.
- **Large / repeat-rich genomes** (*S. lycopersicum* ~830 Mbp / 60% repeats, *H. sapiens* ~3.0 Gbp / 51% repeats): Permissive masking (`>1`) over-fragments long repeat paths; use strict copy thresholds (`>5` or `>10`) to reduce misassemblies while preserving baseline contiguity.
- **Thread Scaling**: Parallelism is controlled via standard OpenMP environment variables (e.g., `OMP_NUM_THREADS=32 ./bin/konuseg.out ...`).

---

## 5. Build, Dependencies, & Verification

### Dependencies
1. **GCC 11 (`g++-11`) with C++17 and OpenMP**: **Specifically GCC 11 (`g++-11`) is required to link against the vendored `lib/libbtllib.a`**, because `lib/libbtllib.a` was compiled as slim LTO bytecode (`.gnu.lto_*`, LTO version 11.3, empty `.text` sections). Newer GCC versions (e.g., GCC 13 with LTO 13.1 on Ubuntu 24.04) will compile all `src/*.cpp` translation units cleanly (`g++ -c`), but the linker (`lto1`) will reject GCC 11.3 LTO bytecode (`bytecode stream in file './lib/libbtllib.a' generated with LTO version 11.3 instead of the expected 13.1`) unless `g++-11` is installed (`sudo apt install g++-11`).
2. **Boost C++ Headers** (`libboost-dev`, optional): Used if installed (`boost/algorithm/string.hpp`, `boost/container/stable_vector.hpp`); `incl/data_structs.hpp` provides standard C++17 fallbacks (`std::vector` and `boost::trim`/`trim_right`) via `__has_include` when Boost is absent.
3. **btllib (v1.4.9/1.4.10)**: Header files (`incl/btllib/`) and static library (`lib/libbtllib.a`) are vendored directly in the repository (**do not modify**).

### Common Commands
```bash
# Build optimized binary
make

# Build with debug symbols
make debug

# Run built-in sanity check on test/test.fa and test/test.fq
make run_test

# Remove compiled binary
make clean
```

---

## 6. Known Implementation Caveats (See `CODEREVIEW.md`)

Developers modifying or extending KoNuSeg should keep the following implementation details in mind (note: **all critical, high, and medium/performance issues `C1`–`C6`, `H1`–`H5`, `M1`–`M6`, and `P1`–`P3` from `CODEREVIEW.md` have been resolved** without modifying `btllib`):
1. **Do Not Modify `btllib`**: Always respect `btllib` API contracts from the caller side (e.g., passing byte counts `(bfBits + 7) / 8` to `SeedBloomFilter` and `cbfSlots * sizeof(uint16_t)` to `CountingBloomFilter16`, using `btllib::DataSource` for compressed stream I/O, and checking `SeedNtHash::roll()` / `!hit_seeds.empty()` on sequences that may contain `N`s).
2. **BED Coordinate Semantics**: Both `<outfile_bed>` and `<outfile_bed>_merged.bed` emit standard **0-based, half-open `[start, end)`** BED intervals directly compatible with `bedtools`.
3. **Global Mutable State**: `k_assembly`, `k_reads`, and `segTable` are defined as `inline` globals in `incl/data_structs.hpp`.
