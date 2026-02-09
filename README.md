# rmolints - Molecular Integrals Library

A high-performance molecular integrals library in Rust, ported from PyQuante2.

## Features

- **Complete integral library**: One-electron and two-electron integrals
- **Multiple ERI algorithms**: Four methods with different performance characteristics
- **Production-ready basis sets**: STO-3G, 6-31G, 6-31G(d), 6-31G(d,p) for elements H-Ar
- **High accuracy**: Matches PyQuante2 to 1e-5 precision
- **Excellent performance**: HGP method is 2x faster than alternatives
- **Parallel computation**: Multithreaded evaluation using Rayon
- **Well-tested**: All tests passing

## Quick Start

```rust
use rmolints::common::*;
use rmolints::{one_electron, hgp};

// Create s-orbital
let s = CGBF {
    origin: Vec3::new(0.0, 0.0, 0.0),
    shell: (0, 0, 0),
    primitives: vec![Primitive { exponent: 1.0, coefficient: 1.0 }],
};

// One-electron integrals
let overlap = one_electron::overlap(&s, &s);
let kinetic = one_electron::kinetic(&s, &s);

// Two-electron integral (recommended method)
let eri = hgp::electron_repulsion_hgp(&s, &s, &s, &s);
```

## Performance

### Production Basis Set Performance (H2O, 6-31G(d,p), 25 basis functions)

| Method | Time (ms) | Status |
|--------|-----------|--------|
| **HGP** | **30** | 🏆 **Fastest - recommended** |
| Rys (optimized) | 55 | Second best |
| Standard THO | 61 | Baseline |

**Recommendation**: Use HGP for all production code.

### Basis Set Comparison

| Molecule | STO-3G Functions | 6-31G(d,p) Functions | ERI Count Ratio |
|----------|------------------|----------------------|-----------------|
| H2       | 2                | 10                   | 500x            |
| H2O      | 7                | 25                   | 250x            |
| Benzene  | 36               | 120                  | 120x            |

**Note**: 6-31G(d,p) is the minimum recommended basis set for chemistry. ERIs scale as O(N⁴) with basis set size.

## Installation

```bash
# Standard build (stable Rust)
cargo build --release

# Run tests
cargo test --release

# Run benchmarks
cargo run --release --example molecule_benchmark
```

### Optional: SIMD Support

SIMD provides minimal benefit on ARM and is opt-in:

```bash
# Requires nightly Rust
cargo +nightly build --release --features simd
```

## Capabilities

### One-Electron Integrals
- Overlap: ⟨bra|ket⟩
- Kinetic: ⟨bra|T|ket⟩
- Nuclear: ⟨bra|V|ket⟩

### Two-Electron Integrals

Four methods available:
1. **HGP** - 🏆 Fastest, recommended (Head-Gordon-Pople)
2. **Rys** - Good for high angular momentum (quadrature)
3. **Standard** - Reference implementation (THO)
4. **HGP-SIMD** - Experimental SIMD variant (nightly required)

### Supported Systems
- **Molecules**: H2, H2O, NH3, benzene (built-in)
- **Basis sets**: STO-3G, 6-31G, 6-31G(d), **6-31G(d,p)** (recommended)
- **Elements**: H through Ar (Z=1-18) for all basis sets
- **Orbitals**: s, p, d (tested); arbitrary angular momentum (theoretical)

## Usage Examples

### Parallel ERI Computation

```rust
use rmolints::parallel::{compute_eri_tensor_parallel, ERIMethod};
use rmolints::molecule::Molecule;
use rmolints::basis::{build_basis, BasisSet};

let molecule = Molecule::h2o();
// Use production-quality 6-31G(d,p) basis set
let basis = build_basis(&molecule, BasisSet::_631GStarStar);

// Compute all unique ERIs (exploits 8-fold symmetry)
let eris = compute_eri_tensor_parallel(&basis, ERIMethod::HeadGordonPople);
```

### Benchmarking

```bash
# Comprehensive molecular benchmarks
cargo run --release --example molecule_benchmark

# Performance ratio analysis
cargo run --release --example ratio_analysis

# Profiling with flamegraphs
cargo run --release --example profile_hgp
```

## Performance Insights

### Key Learnings

1. **Memory layout matters**: Efficient array indexing with pre-computed strides
2. **Algorithm selection**: HGP beats Rys and Standard on real molecules
3. **SIMD limitations**: Minimal benefit for typical quantum chemistry integrals
4. **Profiling results**: VRR stages 6-7 account for 75% of computation time

### Optimization History

- **HGP implementation**: Highly optimized VRR tensor with pre-computed strides
- **Rys optimizations**: Improved storage and caching (1.9x speedup)
- **Overall**: Current HGP is 2x faster than alternatives

## Testing

```bash
# All tests
cargo test --release

# Module-specific
cargo test --release hgp
cargo test --release parallel
```

All 50 tests pass with 1e-5 precision against PyQuante2 reference values.

## Documentation

See [`docs/DOCS.md`](docs/DOCS.md) for:
- Complete installation instructions
- Detailed performance analysis
- Algorithm descriptions
- Extended usage examples
- Profiling and optimization guide
- Next steps and roadmap

## References

- **PyQuante2**: https://github.com/rpmuller/pyquante2
- **THO**: Taketa, Huzinaga, O-ohata, *J. Phys. Soc. Japan* **21**, 2313 (1966)
- **Rys**: Augspurger et al., *J. Comp. Chem.* **11**(8), 972-977 (1990)
- **HGP**: Head-Gordon & Pople, *J. Chem. Phys.* **89**(9), 5777 (1988)

## License

Modified BSD license (same as PyQuante2)
