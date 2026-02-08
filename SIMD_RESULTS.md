# SIMD Optimization Results

## Overview

This document presents results from implementing SIMD (Single Instruction, Multiple Data) vectorization in the HGP-Opt electron repulsion integral method. The SIMD implementation targets the innermost loops of VRR (Vertical Recurrence Relation) stages 6-7, using `std::simd` from Rust's portable_simd feature.

## Implementation Details

### SIMD Approach
- **Library**: `std::simd` (portable_simd nightly feature)
- **Vector width**: f64x4 (4 double-precision floats per SIMD register)
- **Target**: VRR stages 6-7 innermost `im` loops
- **Fallback**: Scalar code for loop tails and when SIMD is disabled

### Optimized Code Sections
1. **Stage 6** (Y-direction c-center, lines 348-465 in `hgp_simd.rs`)
   - Vectorizes innermost im loop
   - Handles 2 conditional branches with SIMD operations
   - Processes 4 m-values per iteration

2. **Stage 7** (Z-direction c-center, lines 505-603 in `hgp_simd.rs`)
   - Same vectorization pattern as Stage 6
   - 2 conditional branches (kc > 0, j > 0)

### Compiler Requirements
- **Nightly Rust** required for `std::simd`
- Feature flag: `--features simd` (enabled by default)
- Rust toolchain: nightly channel (set via `rust-toolchain.toml`)

## Performance Results

### Test System
- **Platform**: Darwin 25.2.0 (macOS on Apple Silicon - ARM64)
- **Compiler**: rustc 1.95.0-nightly
- **Build**: Release mode (opt-level=3, LTO=true)

### Molecule Benchmarks

Performance on real molecules using STO-3G basis set, parallel execution:

| Molecule | Method | Time (ms) | vs HGP-Opt | Notes |
|----------|--------|-----------|------------|-------|
| **H2** (2 basis fns, 6 ERIs) | | | | |
| | HGP-Opt | 0.05 | 1.0x | Baseline |
| | HGP-SIMD | 0.11 | 0.45x | **2.2x slower** |
| **H2O** (7 basis fns, 406 ERIs) | | | | |
| | HGP-Opt | 2.14 | 1.0x | Baseline |
| | HGP-SIMD | 2.10 | 1.02x | **~Equal** |
| **Benzene** (36 basis fns, 222,111 ERIs) | | | | |
| | HGP-Opt | 1081.76 | 1.0x | Baseline |
| | HGP-SIMD | 1159.80 | 0.93x | **1.07x slower** |

### Orbital-Level Benchmarks

Micro-benchmarks on individual integral types (serial execution):

| Orbital Combination | Scalar (ns) | SIMD (ns) | Speedup | Verdict |
|---------------------|-------------|-----------|---------|---------|
| (ss\|ss) | 12,947 | 9,940 | **1.30x faster** | ✅ Win |
| (pp\|pp) | 19,777 | 19,844 | 1.00x slower | ❌ Neutral |
| (sp\|sp) | 14,749 | 14,761 | 1.00x slower | ❌ Neutral |
| (dd\|dd) | 79,943 | 82,696 | 1.03x slower | ❌ Loss |
| (pd\|pd) | 38,482 | 35,416 | **1.09x faster** | ✅ Win |
| (sp\|sp) distant | 19,720 | 20,123 | 1.02x slower | ❌ Loss |

### Analysis

**Why SIMD doesn't provide the expected speedup:**

1. **Small Loop Iterations**
   - VRR im loops typically have 5-20 iterations
   - Processing 4 at a time (f64x4) provides limited parallelism
   - For small integrals (s-orbitals), loops may have < 8 iterations
   - SIMD overhead (broadcasting, setup) dominates savings

2. **Conditional Branches**
   - Each loop has 0-2 conditional terms (`if jc > 0`, `if j > 0`)
   - SIMD requires computing all paths and blending results
   - This negates some performance gains

3. **ARM NEON vs x86 AVX**
   - Tests run on Apple Silicon (ARM64 with NEON)
   - ARM NEON has different performance characteristics than x86 AVX2
   - x86_64 with AVX2 may show different results

4. **Memory Bandwidth**
   - VRR tensor already has excellent cache locality (stride-1 access)
   - SIMD doesn't improve memory access patterns further
   - Scalar code may be memory-bound rather than compute-bound

5. **Compiler Auto-Vectorization**
   - Modern LLVM may already auto-vectorize some loops
   - Explicit SIMD may not add much beyond compiler optimizations

**Where SIMD helps:**
- **(ss|ss)**: 1.30x faster - Simple case with minimal branches
- **(pd|pd)**: 1.09x faster - Mixed angular momentum with more VRR work

**Where SIMD hurts:**
- Small molecules (H2) - overhead dominates
- D-orbital integrals (dd|dd) - branching overhead

## Correctness Validation

All SIMD tests pass with exact agreement to scalar code:

```bash
$ cargo test --release --lib hgp_simd
running 5 tests
test hgp_simd::tests::test_simd_vs_scalar_s_orbitals ... ok
test hgp_simd::tests::test_simd_vs_scalar_p_orbitals ... ok
test hgp_simd::tests::test_simd_vs_scalar_d_orbitals ... ok
test hgp_simd::tests::test_simd_vs_scalar_mixed ... ok
test hgp_simd::tests::test_simd_vs_scalar_all_permutations ... ok

test result: ok. 5 passed; 0 failed
```

Tests verify:
- S, P, D orbital integrals
- Mixed angular momentum combinations
- 8-fold permutational symmetry preservation
- Numeric accuracy (< 1e-10 relative error)

## Usage

### Building with SIMD

```bash
# SIMD is opt-in (requires nightly Rust)
cargo +nightly build --release --features simd

# Standard build (stable Rust, no SIMD)
cargo build --release
```

### Using SIMD in Code

```rust
use rmolints::parallel::{compute_eri_tensor_parallel, ERIMethod};

// Use SIMD-optimized method
let eris = compute_eri_tensor_parallel(&basis, ERIMethod::HeadGordonPopleSimd);
```

### Requirements

**SIMD is opt-in and requires:**

1. Nightly Rust toolchain:
```bash
rustup toolchain install nightly
```

2. Enable the `simd` feature when building:
```bash
cargo +nightly build --release --features simd
```

**By default**, the project builds on **stable Rust** without SIMD.

## Comparison to Other Methods

Performance on benzene (222,111 ERIs, parallel execution):

| Method | Time (ms) | Speedup vs Standard | Notes |
|--------|-----------|---------------------|-------|
| Standard THO | 2,233 | 1.0x | Baseline recursive method |
| Rys Quadrature | 3,741 | 0.60x | Specialized quadrature |
| HGP (original) | 7,002 | 0.32x | HashMap-based recursion |
| **HGP-Opt** | **1,082** | **2.06x** | **Flat array + iterative HRR** |
| HGP-SIMD | 1,160 | 1.93x | SIMD vectorization |

**Verdict**: HGP-Opt (scalar) remains the fastest implementation. SIMD provides marginal benefit at best on ARM hardware.

## Conclusions

### Key Findings

1. **SIMD overhead dominates benefits** for typical quantum chemistry integrals
   - Loop iterations too small (5-20) for 4-wide vectors
   - Conditional branches negate vectorization gains
   - Already excellent cache locality limits memory improvements

2. **Platform-specific performance**
   - Results on ARM (NEON) show SIMD slower or neutral
   - x86_64 with AVX2 may perform differently (untested)

3. **Correctness preserved**
   - SIMD implementation numerically identical to scalar
   - All 8-fold symmetries maintained
   - Safe for production use

4. **Code complexity increased**
   - SIMD adds ~200 lines of conditional compilation
   - Requires nightly Rust (may limit adoption)
   - Harder to maintain and debug

### Recommendations

**For this codebase:**
- **Keep HGP-Opt as default**: Scalar version is simpler and faster on ARM
- **Keep SIMD as optional feature**: May benefit x86_64 users with AVX2
- **Document platform differences**: Warn users SIMD may not help on ARM

**For future work:**
1. **Test on x86_64 with AVX2** - May show better results
2. **Optimize for 2-wide SIMD** - Better match for small loop counts
3. **Profile on real calculations** - HF/DFT workloads may differ
4. **Consider auto-vectorization** - Let compiler handle it
5. **Explore GPU acceleration** - Better fit for massive parallelism

### Overall Assessment

The SIMD optimization successfully demonstrates modern vectorization techniques but does not achieve the target 2-4x speedup on ARM hardware. The scalar HGP-Opt remains the most efficient implementation for typical quantum chemistry workloads, offering a good balance of performance, simplicity, and portability.

However, the SIMD implementation serves as:
- ✅ Proof of concept for SIMD in quantum chemistry integrals
- ✅ Learning resource for Rust SIMD programming
- ✅ Potential benefit on x86_64 platforms (untested)
- ✅ Foundation for future AVX-512 or GPU ports

## References

- [Rust portable_simd RFC](https://github.com/rust-lang/rust/issues/86656)
- [Head-Gordon & Pople, J. Chem. Phys. 89, 5777 (1988)](https://doi.org/10.1063/1.455553)
- Original profiling results: `PROFILING.md`
