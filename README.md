# KoNuSeg

[![GitHub Release](https://img.shields.io/github/v/release/BilkentCompGen/KoNuSeg?color=blue)](https://github.com/BilkentCompGen/KoNuSeg/releases)
[![Last Commit](https://img.shields.io/github/last-commit/BilkentCompGen/KoNuSeg)](https://github.com/BilkentCompGen/KoNuSeg/commits/main)
[![License](https://img.shields.io/github/license/BilkentCompGen/KoNuSeg)](https://github.com/BilkentCompGen/KoNuSeg/blob/main/LICENSE)

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
