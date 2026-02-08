# HGP Boundary Condition Bug Fix

## Problem

Both HGP implementations (original and optimized) crashed when computing integrals for real molecules with mixed s/p orbital basis sets:

- **HGP Original**: "no entry found for key" - HashMap lookup failure
- **HGP Optimized**: "index out of bounds: the len is 1 but the index is 18446744073709551615"

The second error shows the telltale sign: 18446744073709551615 is `usize::MAX`, which is what you get when `-1` wraps around in unsigned arithmetic.

## Root Cause

In the VRR (Vertical Recursion Relation) z-direction recursion for the c-center, both implementations had:

```rust
if k > 0 {
    val += k as f64 / (2.0 * (zeta + eta)) * cache[&(i, j - 1, k, ic, jc, kc, im + 1)];
}
```

**Bug**: The code checks `if k > 0` but then accesses `j - 1`. When `j = 0` and `k > 0`, this evaluates `0 - 1` which wraps to `usize::MAX` in unsigned arithmetic.

## The Fix

Simple one-line change in both files:

**Before**:
```rust
if k > 0 {  // WRONG: checks k but accesses j-1
    val += k as f64 / (2.0 * (zeta + eta)) * vrr[...][j_u - 1][...];
}
```

**After**:
```rust
if j > 0 {  // CORRECT: checks j before accessing j-1
    val += j as f64 / (2.0 * (zeta + eta)) * vrr[...][j_u - 1][...];
}
```

## Files Changed

1. **src/hgp_opt.rs:356** - Fixed boundary check in z-direction VRR recursion
2. **src/hgp.rs:321** - Fixed boundary check in z-direction VRR recursion

Note: The x-direction (checking `i > 0` before `i - 1`) and y-direction (checking `j > 0` before `j - 1`) were already correct.

## Why This Bug Wasn't Caught Earlier

The original micro-benchmarks used simple test cases like:
- Four s-orbitals at the same position
- Four s-orbitals at different positions
- Simple p-orbital combinations

These cases often had `j = 0` when `k = 0`, so the buggy branch wasn't executed. Real molecules with mixed orbital types (like H2O with O having 2s and 2p orbitals) exercise more code paths and exposed the bug.

## Verification

All tests pass:
```bash
cargo test --release
# test result: ok. 48 passed; 0 failed
```

Real molecule benchmarks now complete successfully:
- H2 (2 basis functions): ✅ Works
- H2O (7 basis functions): ✅ Works
- Benzene (36 basis functions): ✅ Works

## Performance After Fix

Both HGP methods now work on real molecules, but remain significantly slower than Standard THO and Rys Quadrature:

| Method | H2O (7 fns) | Benzene (36 fns) | vs Rys |
|--------|-------------|------------------|---------|
| Rys Quadrature | 2.42 ms | 1408 ms | 1.00x ✅ |
| Standard THO | 2.64 ms | 1693 ms | 1.17x |
| HGP Optimized | 3.19 ms | 1994 ms | 1.42x |
| HGP Original | 7.72 ms | 5583 ms | 3.97x ❌ |

**Conclusion**: The bug fix makes HGP methods work correctly, but they remain too slow for production use. Rys Quadrature is the clear winner for real molecular calculations.
