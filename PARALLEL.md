# Parallel Computation in rmolints

## Overview

The `parallel` module provides multithreaded computation of molecular integrals using [Rayon](https://github.com/rayon-rs/rayon) for data parallelism. This is particularly useful for computing large batches of two-electron repulsion integrals (ERIs).

## Quick Start

```rust
use rmolints::parallel::{compute_eri_tensor_parallel, ERIMethod};
use rmolints::common::*;

// Create basis set
let s1 = CGBF {
    origin: Vec3::new(0.0, 0.0, 0.0),
    shell: (0, 0, 0),
    primitives: vec![Primitive { exponent: 1.0, coefficient: 1.0 }],
};
let s2 = CGBF {
    origin: Vec3::new(0.0, 0.0, 1.0),
    shell: (0, 0, 0),
    primitives: vec![Primitive { exponent: 1.0, coefficient: 1.0 }],
};
let basis = vec![s1, s2];

// Compute all unique ERIs in parallel
let eris = compute_eri_tensor_parallel(&basis, ERIMethod::Standard);

// eris is a Vec<(usize, usize, usize, usize, f64)>
// Each element is (i, j, k, l, integral_value)
```

## Available Functions

### 1. `compute_eris_parallel`

Compute a batch of ERIs for specified index quartets in parallel.

```rust
let indices = vec![(0, 0, 0, 0), (0, 1, 0, 1), (1, 1, 1, 1)];
let results = compute_eris_parallel(&basis, &indices, ERIMethod::Standard);
// results[0] corresponds to indices[0], etc.
```

**Best for**: Custom integral lists, screening, or when you know exactly which integrals you need.

### 2. `compute_eri_tensor_parallel`

Compute all unique ERIs exploiting 8-fold permutational symmetry.

```rust
let eris = compute_eri_tensor_parallel(&basis, ERIMethod::Standard);
// Returns Vec<(i, j, k, l, value)> with i>=j, k>=l, ij>=kl
```

**Best for**: Computing complete ERI tensors for quantum chemistry calculations.

**Symmetry**: Exploits (ij|kl) = (ji|kl) = (ij|lk) = (ji|lk) = (kl|ij) = (lk|ij) = (kl|ji) = (lk|ji)

### 3. `compute_eri_tensor_full_parallel`

Compute all N^4 integrals without symmetry reduction.

```rust
let tensor = compute_eri_tensor_full_parallel(&basis, ERIMethod::Standard);
// Returns Vec<Vec<Vec<Vec<f64>>>> indexed as [i][j][k][l]
```

**⚠️ Warning**: This computes N^4 integrals and can be extremely large!
- N=10: 10,000 integrals (~80 KB)
- N=50: 6.25 million integrals (~50 MB)
- N=100: 100 million integrals (~800 MB)

**Best for**: Direct indexing needs or when symmetry exploitation is not desired.

## Choosing an ERI Method

Three methods are available via the `ERIMethod` enum:

```rust
pub enum ERIMethod {
    Standard,         // Taketa-Huzinaga-O-ohata
    Rys,              // Rys quadrature
    HeadGordonPople,  // Head-Gordon-Pople HRR/VRR
}
```

**Recommendations**:
- **Standard**: Default choice, reliable and accurate for all cases
- **Rys**: Best for high angular momentum (d, f orbitals) when accuracy is acceptable
- **HeadGordonPople**: Fast for s-orbitals, but slow for p/d orbitals (see BENCHMARK_RESULTS.md)

## Performance Characteristics

### When Parallel Helps

Parallel computation provides speedup when:

1. **Large number of integrals**: More integrals = better parallelization
2. **Complex integrals**: Higher angular momentum = more work per integral
3. **Available CPU cores**: More cores = more potential speedup

### Benchmark Results

From `examples/parallel_benchmark.rs`:

| Basis Set | Integrals | Serial Time | Parallel Time | Speedup |
|-----------|-----------|-------------|---------------|---------|
| 2 s-orbitals | 6 | 0.01 ms | 0.20 ms | 0.07x ❌ |
| 4 s-orbitals | 55 | 0.02 ms | 0.10 ms | 0.20x ❌ |
| 8 s-orbitals | 666 | 0.19 ms | 0.63 ms | 0.30x ❌ |
| 4s + 4p mixed | 666 | 0.33 ms | 0.23 ms | 1.46x ✅ |

**Key Findings**:
- Small batches (<100 integrals): Parallel overhead dominates, serial is faster
- Large batches with simple integrals: Thread spawning overhead still significant
- Complex integrals (p/d orbitals): Parallel begins to show advantage

### Overhead Analysis

Parallel computation has two sources of overhead:

1. **Thread spawning**: ~0.1-0.2 ms fixed cost
2. **Work distribution**: Rayon's work-stealing scheduler overhead

For very fast integrals (~100-500 ns each), you need thousands of integrals to amortize the thread spawning cost.

## When to Use Parallel

### ✅ Use Parallel When:
- Computing >1000 integrals
- Working with p, d, or f orbitals (higher angular momentum)
- Building full ERI tensors for molecules with 10+ basis functions
- Integral evaluation is the bottleneck in your calculation

### ❌ Use Serial When:
- Computing <100 integrals
- Working with only s-orbitals
- Small basis sets (<5 basis functions)
- Need predictable, deterministic timing

## Example: Real Molecular Calculation

```rust
use rmolints::parallel::{compute_eri_tensor_parallel, ERIMethod};
use rmolints::common::*;

// Water molecule with 6-31G basis (13 basis functions)
// ... construct basis set ...

let basis = vec![
    // 1s, 2s, 2px, 2py, 2pz for oxygen (5 functions)
    // ...
    // 1s for each hydrogen (2 functions each, 4 total)
    // ...
];

// Compute all unique ERIs in parallel
// For 13 basis functions: 13^4 / 8 = ~3,605 unique integrals
let start = Instant::now();
let eris = compute_eri_tensor_parallel(&basis, ERIMethod::Standard);
let elapsed = start.elapsed();

println!("Computed {} integrals in {:.2} ms", eris.len(), elapsed.as_secs_f64() * 1000.0);
// Expected: ~5-10 ms depending on hardware

// Use the integrals in Hartree-Fock or other calculations
for (i, j, k, l, value) in eris {
    // Build Fock matrix, etc.
}
```

## Thread Configuration

Rayon automatically detects and uses all available CPU cores. To customize:

```rust
// Set number of threads (do this before any parallel operations)
rayon::ThreadPoolBuilder::new()
    .num_threads(4)
    .build_global()
    .unwrap();
```

Or use the `RAYON_NUM_THREADS` environment variable:
```bash
RAYON_NUM_THREADS=4 cargo run --release
```

## Memory Considerations

Parallel computation uses more memory than serial:

- Each thread needs stack space (~2-8 MB per thread)
- Rayon's work-stealing queues add overhead
- Full tensors can be very large (N^4 scaling)

**Recommendations**:
- For large basis sets (N > 50), use symmetry-reduced `compute_eri_tensor_parallel`
- Consider computing integrals in chunks if memory is limited
- Monitor memory usage for N > 100

## Future Optimizations

Potential improvements for parallel performance:

- [ ] **Chunking**: Process integrals in larger chunks to reduce overhead
- [ ] **Rayon chunk size tuning**: Optimize work distribution
- [ ] **Nested parallelism**: Parallelize within individual integrals for very high L
- [ ] **Integral screening**: Skip near-zero integrals
- [ ] **Cache-aware scheduling**: Group spatially nearby integrals
- [ ] **GPU acceleration**: Offload to GPU for very large basis sets

## Comparison with Other Libraries

| Library | Language | Parallelization | Notes |
|---------|----------|-----------------|-------|
| **rmolints** | Rust | Rayon (CPU) | This library |
| Libint2 | C++ | OpenMP | Industry standard |
| PyQuante2 | Python/C | None | Reference implementation |
| PySCF | Python/C | OpenMP | Production quantum chemistry |

rmolints provides similar parallelization to Libint2's OpenMP, but using Rust's safe parallelism model.

## Testing

All parallel functions have corresponding tests ensuring:
- Parallel results match serial results exactly
- All three methods (Standard, Rys, HGP) work in parallel
- Symmetry exploitation is correct

Run parallel tests:
```bash
cargo test parallel
```

Run parallel benchmarks:
```bash
cargo run --release --example parallel_benchmark
```

## Summary

**Parallel computation in rmolints**:
- ✅ Easy to use with Rayon
- ✅ Safe and data-race free (Rust guarantees)
- ✅ Automatic work distribution
- ⚠️ Only beneficial for large batches or complex integrals
- ⚠️ Thread spawning overhead for small problems

For production quantum chemistry calculations with moderate to large basis sets (N ≥ 10), parallel computation provides meaningful speedup and should be used by default.
