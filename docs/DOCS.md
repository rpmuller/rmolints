# rmolints - Complete Documentation

## Table of Contents

1. [Overview](#overview)
2. [Installation](#installation)
3. [Quick Start](#quick-start)
4. [Capabilities](#capabilities)
5. [Performance Analysis](#performance-analysis)
6. [Algorithm Details](#algorithm-details)
7. [Usage Examples](#usage-examples)
8. [Testing](#testing)
9. [Profiling and Optimization](#profiling-and-optimization)
10. [Next Steps](#next-steps)

---

## Overview

**rmolints** is a high-performance molecular integrals library written in Rust, ported from PyQuante2. It provides production-ready implementations of one-electron and two-electron repulsion integrals for quantum chemistry calculations.

### Key Features

- **Complete integral library**: All fundamental one- and two-electron integrals
- **Multiple ERI algorithms**: Five different methods thoroughly benchmarked
- **Schwarz screening**: Skip near-zero shell quartets for up to 2x speedup on larger systems
- **Production-ready basis sets**: STO-3G, 6-31G, 6-31G(d), 6-31G(d,p) for elements H-Ar
- **High accuracy**: Matches PyQuante2 reference to 1e-5 precision
- **Excellent performance**: Optimized Rys is the fastest method (~30% faster than HGP)
- **Parallel computation**: Multithreaded ERI evaluation using Rayon
- **Real molecule support**: H2, H2O, NH3, benzene with multiple basis sets
- **Well-tested**: 56 comprehensive tests covering all modules

### Project Status

✅ **Complete and production-ready** for s, p, d orbitals with STO-3G and 6-31G family basis sets.

---

## Installation

### Requirements

- **Rust**: 1.70 or later (stable)
- **Optional**: Nightly Rust for SIMD features

### Standard Build

```bash
# Clone repository
git clone <repository-url>
cd rmolints

# Build in release mode (required for performance)
cargo build --release

# Run tests
cargo test --release
```

### With SIMD (Optional, Nightly Required)

SIMD provides marginal benefits on x86_64 but is slower on ARM hardware. Use only if you have confirmed benefits on your platform.

```bash
# Install nightly Rust
rustup toolchain install nightly

# Build with SIMD support
cargo +nightly build --release --features simd

# Run SIMD-specific benchmarks
cargo +nightly run --release --features simd --example simd_benchmark
```

---

## Quick Start

### Basic Usage

```rust
use rmolints::common::*;
use rmolints::{one_electron, hgp};

// Create s-orbital centered at origin
let origin = Vec3::new(0.0, 0.0, 0.0);
let s_orbital = CGBF {
    origin,
    shell: (0, 0, 0),  // s-orbital (l=m=n=0)
    primitives: vec![Primitive {
        exponent: 1.0,
        coefficient: 1.0,
    }],
};

// One-electron integrals
let overlap = one_electron::overlap(&s_orbital, &s_orbital);
let kinetic = one_electron::kinetic(&s_orbital, &s_orbital);
let nuclear = one_electron::nuclear_attraction(&s_orbital, &s_orbital, &origin);

// Two-electron integral (recommended method)
let eri = hgp::electron_repulsion_hgp(&s_orbital, &s_orbital,
                                               &s_orbital, &s_orbital);
```

### Parallel Computation

```rust
use rmolints::parallel::{
    compute_eri_tensor_parallel, compute_eri_tensor_screened_parallel,
    ERIMethod, SCHWARZ_THRESHOLD,
};
use rmolints::molecule::Molecule;
use rmolints::basis::{build_basis, BasisSet};

let molecule = Molecule::h2o();
let basis = build_basis(&molecule, BasisSet::_631GStarStar);

// With Schwarz screening (recommended for larger systems: skips 40-60% of integrals)
let eris = compute_eri_tensor_screened_parallel(
    &basis, ERIMethod::HeadGordonPopleContracted, SCHWARZ_THRESHOLD,
);

// Without screening (simpler, fine for small/trivial systems)
let eris = compute_eri_tensor_parallel(&basis, ERIMethod::HeadGordonPopleContracted);

// Returns Vec<(i, j, k, l, value)> for all (surviving) unique integrals
println!("Computed {} unique ERIs", eris.len());
```

### Running Benchmarks

```bash
# Schwarz screening impact (screened vs unscreened, multiple basis sets)
cargo run --release --example schwarz_benchmark

# Comprehensive molecular benchmark
cargo run --release --example molecule_benchmark

# Orbital-level micro-benchmarks
cargo run --release --example benchmark

# Parallel vs serial comparison
cargo run --release --example parallel_benchmark

# Performance ratio analysis
cargo run --release --example ratio_analysis
```

---

## Capabilities

### One-Electron Integrals

All integrals use the Taketa-Huzinaga-O-ohata (THO) method:

| Type | Function | Description |
|------|----------|-------------|
| **Overlap** | `overlap(bra, ket)` | ⟨bra\|ket⟩ |
| **Kinetic** | `kinetic(bra, ket)` | ⟨bra\|T\|ket⟩ |
| **Nuclear** | `nuclear_attraction(bra, ket, center)` | ⟨bra\|V\|ket⟩ |

### Two-Electron Integrals (ERIs)

Five methods available with different performance characteristics:

| Method | Function | Performance | Status |
|--------|----------|-------------|--------|
| **Rys** | `rys::electron_repulsion_rys()` | **Fastest (~30% faster than HGP)** | 🏆 **Recommended** |
| **HGP-Contracted** | `hgp::electron_repulsion_hgp_contracted()` | Fast (VRR-contracted) | Production ready |
| **HGP** | `hgp::electron_repulsion_hgp()` | Good (baseline) | Debugging/verification |
| **Standard** | `two_electron::electron_repulsion()` | 2.2x slower than Rys | Reference impl |
| **HGP-SIMD** | `hgp_simd::electron_repulsion_hgp_simd()` | ~Equal on ARM | Optional (nightly) |

**Recommendation**: Use **Rys** for production code.

### Supported Molecules

Built-in molecule definitions with STO-3G basis sets:

```rust
use rmolints::molecule::Molecule;

let h2 = Molecule::h2(1.4);           // H2 with custom bond length
let water = Molecule::h2o();          // H2O
let benzene = Molecule::benzene();    // C6H6
```

### Basis Set Support

- **Available basis sets**: STO-3G, 6-31G, 6-31G(d), 6-31G(d,p)
- **Elements**: H through Ar (Z=1-18) for all basis sets
- **Recommended for production**: 6-31G(d,p) - industry standard for chemical accuracy
- **Angular momentum**: s, p, d orbitals tested; arbitrary L supported

---

## Performance Analysis

### STO-3G Performance (parallel execution)

| Molecule | Rys (ms) | HGP Contracted (ms) | Rys vs HGP |
|----------|----------|----------------------|------------|
| H2       | 0.094    | 0.111                | 0.84x (faster) |
| H2O      | 0.676    | 0.947                | 0.71x (faster) |
| NH3      | 0.937    | 1.362                | 0.69x (faster) |
| Benzene  | 231.6    | 333.6                | 0.69x (faster) |

Rys is now **~30% faster** than HGP Contracted across all molecule sizes.

### Basis Set Scaling

ERIs scale as O(N⁴) with basis set size. Moving from STO-3G to 6-31G(d,p):

| Molecule | STO-3G Functions | 6-31G(d,p) Functions | STO-3G ERIs | 6-31G(d,p) ERIs | Ratio |
|----------|------------------|----------------------|-------------|-----------------|-------|
| H2       | 2                | 10                   | ~6          | ~3,000          | 500x  |
| H2O      | 7                | 25                   | ~400        | ~100,000        | 250x  |
| Benzene  | 36               | 120                  | ~220,000    | ~26,000,000     | 120x  |

**Key insight**: 6-31G(d,p) increases computational cost by 100-500x but is necessary for chemical accuracy and is the minimum recommended basis set for published work.

### Key Performance Learnings

#### 1. Memory Layout Matters

**Efficient array indexing vs nested structures**: The HGP method achieves a **2.4x speedup** over nested Vec storage by using:
- Pre-allocated efficient array indexing with computed strides
- Cache-optimal memory layout (stride-1 innermost loops)
- Eliminated allocation overhead in hot loops

#### 2. Algorithm Selection

Different ERI methods have distinct performance profiles:

- **Rys**: Best overall — fastest method after hot-path optimization, good for all angular momenta
- **HGP-Contracted**: Close second, VRR-level contraction reduces redundant HRR work
- **HGP**: Original implementation, good for debugging
- **Standard**: Simplest code, good reference implementation

**Performance ratio (independent of system variance):**
- HGP Contracted is ~1.4x slower than Rys (stable across systems)
- Standard is ~2.2x slower than Rys (stable across systems)

#### 3. SIMD Limitations

SIMD vectorization (f64x4) provides **minimal benefit** for typical quantum chemistry integrals:

**Why SIMD doesn't help:**
- Small loop iterations (5-20) limit parallelism opportunity
- Conditional branches require SIMD blending (overhead)
- Memory already cache-optimal in scalar code
- ARM NEON shows slower performance than scalar

**Verdict**: Scalar HGP is simpler and faster. SIMD is optional for x86_64 experimentation.

#### 4. Profiling Results

CPU profiling identified the computational bottlenecks in HGP:

- **VRR tensor computation**: 75-85% of total time
  - Stage 6 (Y-direction): ~25-30% of VRR time
  - Stage 7 (Z-direction): ~40-50% of VRR time
- **HRR recursion**: ~10-15% of total time
- **Helper functions** (Boys function, etc.): ~5-10% of total time

**Optimization focus**: VRR stages 6-7 are the primary target for any future improvements.

#### 5. Rys Optimizations

Two rounds of optimization brought Rys from 3.2x slower than HGP to **~30% faster** than HGP:

**Round 1 (1.9x speedup):**
1. **Flat GMatrix storage**: Replaced `Vec<Vec<f64>>` with efficient array indexing
2. **Cached product centers**: Eliminated 75% of redundant Gaussian product calculations

**Round 2 (1.8–2.3x additional speedup):**
3. **GMatrix reuse**: Pre-allocate once per contracted integral, reuse across all primitive quartets and quadrature points (eliminated ~243 allocations per integral for STO-3G water)
4. **Stack-allocated Rys roots/weights**: Fixed-size `[f64; 8]` arrays instead of `Vec<f64>`, eliminating 2 heap allocations per primitive quartet
5. **Hoisted normalization**: `gaussian_normalization` computed once per primitive, not once per quartet (O(n) instead of O(n⁴))

Combined effect: Rys is now the fastest ERI method, ~30% faster than HGP Contracted.

### Scaling Characteristics

ERI computation scales as O(N⁴) with basis set size N:

**STO-3G Basis Set:**
| Molecule | Basis Functions | Unique ERIs | Ratio |
|----------|----------------|-------------|-------|
| H2 | 2 | 6 | 3:1 |
| H2O | 7 | 406 | 58:1 |
| Benzene | 36 | 222,111 | 6,170:1 |

**6-31G(d,p) Basis Set:**
| Molecule | Basis Functions | Unique ERIs | Ratio |
|----------|----------------|-------------|-------|
| H2 | 10 | 3,025 | 303:1 |
| H2O | 25 | 52,975 | 2,119:1 |
| Benzene | 120 | ~26,000,000 | ~217,000:1 |

**Parallel speedup**: ~3.4x on multi-core system (measured on benzene).

---

## Algorithm Details

### Head-Gordon-Pople (HGP)

The recommended method uses:

1. **VRR (Vertical Recurrence Relation)**: Build intermediate tensor in 7 stages
   - Stages 1-3: A-center recursion (x, y, z)
   - Stages 4-7: C-center recursion (x, y, z)
   - Efficient array indexing storage with pre-computed strides

2. **HRR (Horizontal Recurrence Relation)**: Convert VRR tensor to final integrals
   - Iterative implementation (not recursive)
   - Minimal computational cost (~10-15% of time)

**Key innovation**: VRR tensor with pre-computed strides for efficient array indexing.

### Rys Quadrature

Uses polynomial quadrature with Rys roots and weights:

- Numerical quadrature of incomplete gamma function
- Good for all angular momenta, now the fastest method overall
- Optimizations: GMatrix reuse, stack-allocated roots/weights, hoisted normalization
- Zero heap allocations in the hot primitive loop

**Reference**: Augspurger, Bernholdt, Dykstra, *J. Comp. Chem.* **11**(8), 972-977 (1990).

### Standard THO Method

Traditional recursive method from Taketa, Huzinaga, O-ohata:

- Clear algorithm, easy to understand
- B-array recursion with incomplete gamma function
- Solid reference implementation

**Reference**: Taketa, Huzinaga, O-ohata, *J. Phys. Soc. Japan* **21**, 2313 (1966).

### Schwarz Inequality Screening

Before computing a shell quartet (ij|kl), apply the Cauchy-Schwarz bound:

```
|(ij|kl)| ≤ Q_ij · Q_kl   where Q_ij = √|(ij|ij)|
```

**Algorithm**:
1. Precompute the N×N Schwarz matrix: `Q_ij = √|(ij|ij)|` for all pairs — O(N²) integrals, parallelized
2. For each candidate quartet, test `Q_ij · Q_kl ≥ threshold` before adding it to the work list
3. Only the surviving quartets are computed

**Cost/benefit**: The diagonal integrals add O(N²) overhead, but screening typically eliminates 40-60% of quartets for medium/large molecules, yielding net ~2x speedup.

| System | Basis | Skipped | Speedup |
|--------|-------|---------|---------|
| H2O | STO-3G | 38% | 1.4x |
| Benzene | STO-3G | 47% | 1.95x |
| H2O | 6-31G(d) | 53% | 1.84x |
| Benzene | 6-31G(d) | 59% | 2.25x |

Default threshold: `SCHWARZ_THRESHOLD = 1e-12` (consistent with 1e-12 integral precision).

**Note**: For very small systems (N < ~5) the O(N²) overhead may exceed the savings — use `compute_eri_tensor_parallel` in that case.

---

## Usage Examples

### Example 1: Water Molecule Energy Components

```rust
use rmolints::molecule::Molecule;
use rmolints::basis::{build_basis, BasisSet};
use rmolints::{one_electron, hgp};

fn main() {
    // Build water molecule with production-quality basis set
    let water = Molecule::h2o();
    let basis = build_basis(&water, BasisSet::_631GStarStar);
    let n = basis.len();

    println!("Water molecule: {} basis functions", n);

    // Compute overlap matrix
    let mut overlap = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            overlap[i][j] = one_electron::overlap(&basis[i], &basis[j]);
        }
    }

    // Compute kinetic energy matrix
    let mut kinetic = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            kinetic[i][j] = one_electron::kinetic(&basis[i], &basis[j]);
        }
    }

    // Compute nuclear attraction matrix
    let mut nuclear = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            let mut v_ij = 0.0;
            for atom in &water.atoms {
                v_ij += atom.Z * one_electron::nuclear_attraction(
                    &basis[i], &basis[j], &atom.position
                );
            }
            nuclear[i][j] = v_ij;
        }
    }

    println!("One-electron integrals computed successfully!");
}
```

### Example 2: Parallel ERI Computation

```rust
use rmolints::parallel::{compute_eri_tensor_parallel, ERIMethod};
use rmolints::molecule::Molecule;
use rmolints::basis::{build_basis, BasisSet};
use std::time::Instant;

fn main() {
    let molecule = Molecule::h2o();
    let basis = build_basis(&molecule, BasisSet::_631GStarStar);

    println!("Computing ERIs for {} basis functions...", basis.len());

    let start = Instant::now();
    let eris = compute_eri_tensor_parallel(&basis, ERIMethod::HeadGordonPople);
    let elapsed = start.elapsed().as_millis();

    println!("Computed {} unique ERIs in {} ms", eris.len(), elapsed);
    println!("Average: {:.2} μs per integral",
             elapsed as f64 * 1000.0 / eris.len() as f64);
}
```

### Example 3: Method Comparison

```rust
use rmolints::parallel::{compute_eris_parallel, ERIMethod};
use rmolints::common::*;
use std::time::Instant;

fn benchmark_method(basis: &[CGBF], indices: &[(usize, usize, usize, usize)],
                    method: ERIMethod) -> f64 {
    let start = Instant::now();
    let _results = compute_eris_parallel(basis, indices, method);
    start.elapsed().as_secs_f64() * 1000.0
}

fn main() {
    // Create test basis and indices...

    let hgp_time = benchmark_method(&basis, &indices, ERIMethod::HeadGordonPople);
    let rys_time = benchmark_method(&basis, &indices, ERIMethod::Rys);
    let std_time = benchmark_method(&basis, &indices, ERIMethod::Standard);

    println!("HGP: {:.2} ms (baseline)", hgp_time);
    println!("Rys:     {:.2} ms ({:.2}x)", rys_time, rys_time / hgp_time);
    println!("Standard: {:.2} ms ({:.2}x)", std_time, std_time / hgp_time);
}
```

---

## Testing

### Test Coverage

```
✅ 56/56 tests passing (100%)

Breakdown:
  - Utility functions:        11 tests
  - One-electron integrals:   10 tests
  - Two-electron (standard):   7 tests
  - Rys quadrature:            3 tests
  - Head-Gordon-Pople:         9 tests
  - Parallel computation:      7 tests
  - Molecule support:          3 tests
  - Basis sets:                5 tests
  - Boys function:             2 tests
```

### Running Tests

```bash
# Run all tests
cargo test --release

# Run specific module tests
cargo test --release one_electron
cargo test --release hgp
cargo test --release parallel

# Run with output
cargo test --release -- --nocapture
```

### Validation Strategy

All implementations are validated against PyQuante2 reference values:

- **Accuracy**: 1e-5 relative precision
- **Coverage**: s, p, d orbital combinations
- **Symmetry**: 8-fold permutational symmetry verified
- **Edge cases**: Zero separation, large exponents, mixed angular momentum

---

## Profiling and Optimization

### CPU Profiling

Generate flamegraphs to identify bottlenecks:

```bash
# Parallel profiling (multi-threaded)
cargo run --release --example profile_hgp

# Serial profiling (single-threaded, clearer results)
cargo run --release --example profile_hgp_serial

# Open generated flamegraphs
open flamegraph.svg
open flamegraph_serial.svg
```

### Optimization History

The project has implemented several key optimizations:

1. **HGP Method**
   - VRR tensor with pre-computed strides for efficient indexing
   - Iterative HRR to avoid repeated VRR calls
   - VRR-level contraction: accumulate VRR tensors before HRR (10-40% faster)
   - Excellent cache locality through contiguous memory layout

2. **Rys Optimization** (cumulative ~4x speedup over initial)
   - Optimized GMatrix storage (contiguous arrays)
   - Cached Gaussian product centers
   - GMatrix reuse across primitive quartets and quadrature points
   - Stack-allocated Rys roots/weights (zero heap allocs in hot loop)
   - Hoisted normalization out of O(n⁴) loop

3. **SIMD Exploration** (no benefit on ARM)
   - Implemented f64x4 vectorization
   - Overhead dominated benefits for small loops
   - Opt-in feature for experimentation

### Performance Measurement Best Practices

When benchmarking, account for system variance:

1. **Use performance ratios** instead of absolute times
2. **Run multiple iterations** and take minimum/average
3. **Include warmup runs** to avoid JIT/cache effects
4. **Compare on same system** to isolate changes

See `examples/ratio_analysis.rs` for reference implementation.

---

## Next Steps

### High Priority (Performance)

1. ~~**Larger Basis Sets**~~ ✅ **COMPLETE**
   - ✅ Implemented 6-31G, 6-31G(d), 6-31G(d,p) for all elements H-Ar
   - ✅ Extended all basis sets to elements H through Ar (Z=1-18)
   - Future: cc-pVDZ, cc-pVTZ, and basis set parsing from standard formats

2. **Integral Screening**
   - Skip near-zero integrals based on Schwarz inequality
   - Critical for large molecules (100+ basis functions)
   - Expected 10-100x speedup for sparse systems

3. **Higher Angular Momentum**
   - Thoroughly test f, g orbitals (L ≥ 4)
   - Add complete Rys polynomial tables (currently simplified)
   - Validate against reference codes

### Medium Priority (Features)

4. **Hartree-Fock Solver**
   - SCF implementation using computed integrals
   - DIIS convergence acceleration
   - Demonstrate end-to-end quantum chemistry calculation

5. **Density Functional Theory (DFT)**
   - Exchange-correlation functionals (LDA, GGA)
   - Numerical quadrature grids
   - Hybrid functionals (B3LYP)

6. **Gradient Integrals**
   - Analytic derivatives of integrals
   - Geometry optimization support
   - Frequency calculations

### Low Priority (Integration)

7. **Python Bindings**
   - PyO3 wrapper for Python interoperability
   - NumPy integration for matrices
   - Compare performance to PySCF, Psi4

8. **GPU Offload**
   - Move VRR computation to GPU (CUDA/ROCm)
   - Batch processing for throughput
   - Target 10-100x speedup for large systems

9. **Benchmark vs Production Codes**
   - Compare to Gaussian, ORCA, Psi4
   - Validate correctness on standard test sets
   - Publish performance comparisons

### Infrastructure

10. **Continuous Integration**
    - Automated testing on multiple platforms
    - Performance regression tracking
    - Documentation generation

11. **Extended Documentation**
    - Tutorial for quantum chemistry beginners
    - API documentation with examples
    - Theory background for each method

12. **Code Organization**
    - Split large modules (rys.rs is 800+ lines)
    - Reduce code duplication across methods
    - Improve error handling and validation

---

## References

### Original Publications

- **THO**: Taketa, Huzinaga, O-ohata, "Gaussian-Expansion Methods for Molecular Integrals", *J. Phys. Soc. Japan* **21**, 2313 (1966).

- **Rys**: Augspurger, Bernholdt, Dykstra, "Concise, open-ended implementation of Rys polynomial evaluation of two-electron integrals", *J. Comp. Chem.* **11**(8), 972-977 (1990).

- **HGP**:
  - Head-Gordon & Pople, "A method for two-electron Gaussian integral and integral derivative evaluation using recurrence relations", *J. Chem. Phys.* **89**(9), 5777 (1988).
  - Gill, "Molecular Integrals Over Gaussian Basis Functions", *Adv. Q. Chem.* **25**, 141 (1994).
  - Gill & Pople, "The Prism Algorithm for Two-Electron Integrals", *IJQC* **40**, 753 (1991).

### Implementation References

- **PyQuante2**: https://github.com/rpmuller/pyquante2
- **Rust SIMD**: https://github.com/rust-lang/rust/issues/86656

---

## License

Same as PyQuante2 (modified BSD license)

---

## Contributing

Contributions welcome! Priority areas:

1. Larger basis sets (6-31G family)
2. Integral screening for large molecules
3. Hartree-Fock solver implementation
4. Validation against production codes

Please include tests and benchmarks with all contributions.
