# HGP Optimization Summary

## Problem Analysis

The original Head-Gordon-Pople (HGP) implementation in `src/hgp.rs` had severe performance issues:

### Original Issues:
1. **HashMap allocation overhead**: Every `vrr()` call created a new `HashMap<(i32,i32,i32,i32,i32,i32,i32), f64>`
2. **Recursive explosion**: HRR made exponentially many recursive calls, each triggering a full VRR computation
3. **No memoization**: Each VRR call built the entire 7-dimensional cache from scratch
4. **Poor cache locality**: HashMap lookups are 10-100x slower than array indexing

### Performance Impact:
```
Test 4 (four p-orbitals):  6,079 ns (13x slower than Standard)
Test 6 (four d-orbitals): 18,048 ns (13x slower than Standard)
```

HGP was **completely unusable** for production work.

## Solution: Optimized Implementation

Created `src/hgp_opt.rs` based on the Julia implementation strategy from `MolecularIntegrals.jl`.

### Key Changes:

1. **Pre-compute VRR tensor once**
   - Old: Rebuild HashMap on every recursive HRR call
   - New: Compute full 7D Vec array once per primitive quartet
   ```rust
   let vrr: Vec<Vec<Vec<Vec<Vec<Vec<Vec<f64>>>>>>> = compute_vrr_tensor(...);
   ```

2. **Use Vec instead of HashMap**
   - Old: `HashMap<(i32,i32,i32,i32,i32,i32,i32), f64>`
   - New: `vrr[la][ma][na][lc][mc][nc][m]`
   - Direct array indexing is **much faster** than hash lookups

3. **Simplified HRR**
   - Old: Deep recursion calling both HRR and VRR
   - New: HRR recursion only, looks up pre-computed VRR values
   - Much shallower recursion tree

### Code Structure:

```rust
fn primitive_eri(...) -> f64 {
    // Compute maximum angular momentum needed
    let la_max = la + lb;
    let lc_max = lc + ld;

    // Compute VRR tensor ONCE
    let vrr_tensor = compute_vrr_tensor(..., la_max, ..., lc_max);

    // Apply HRR using pre-computed VRR
    hrr_iterative(&vrr_tensor, ...)
}
```

## Performance Results

### HGP Original vs Optimized:

| Test Case | Original (ns) | Optimized (ns) | Speedup |
|-----------|---------------|----------------|---------|
| s-orbitals (same) | 303 | 361 | 0.84x ❌ |
| s-orbitals (separated) | 309 | 327 | 0.94x ❌ |
| Mixed s+p | 716 | 625 | 1.15x ✅ |
| **Four p-orbitals** | **6,079** | **1,505** | **4.0x** ✅✅ |
| **d-orbitals** | **1,814** | **910** | **2.0x** ✅ |
| **Four d-orbitals** | **18,048** | **2,123** | **8.5x** ✅✅✅ |

### Comparison vs All Methods (Optimized):

| Method | Test 4 (4p) | Test 5 (d) | Test 6 (4d) |
|--------|-------------|------------|-------------|
| **Standard** | **765 ns** | **673 ns** | **1,372 ns** |
| **Rys** | 743 ns | 678 ns | **1,030 ns** |
| HGP Original | 6,079 ns ❌ | 1,814 ns ❌ | 18,048 ns ❌ |
| HGP Optimized | 1,505 ns | 910 ns | 2,123 ns |

## Key Insights

### 1. Optimization Was Successful
- **2-8.5x speedup** for HGP at medium-high angular momentum
- Fixed the "completely broken" state of original HGP
- Now usable (though not always fastest)

### 2. Method Rankings (Updated)
For general use:
1. **Standard THO**: Consistently fast, reliable, simple
2. **Rys Quadrature**: Best for highest angular momentum (L ≥ 4)
3. **HGP Optimized**: Competitive but not usually fastest
4. ~~HGP Original~~: Deprecated, too slow

### 3. Why HGP Still Isn't Fastest

Even optimized, HGP has overhead:
- **Large tensor allocation**: 7D array for VRR results
- **Memory bandwidth**: Accessing large arrays is slower than smaller ones
- **Recursion depth**: HRR still recurses (though much shallower)

For simple integrals (s-orbitals), this overhead dominates:
- Standard: 230 ns (direct computation)
- HGP Opt: 361 ns (allocate tensor + compute + lookup)

### 4. When to Use Each Method

**Standard THO**:
- Default choice for all cases
- Simple, fast, reliable
- Best for s, p orbitals

**Rys Quadrature**:
- Use for d, f, g orbitals (L ≥ 3)
- 33% faster than Standard at L=4
- Scales better with angular momentum

**HGP Optimized**:
- Academic interest / comparison
- Not recommended for production
- Standard or Rys always better or competitive

## Implementation Details

### VRR Tensor Allocation

```rust
let mut vrr = vec![
    vec![vec![vec![vec![vec![
        vec![0.0; (mtot + 1) as usize];
        (nc_max + 1) as usize
    ]; (mc_max + 1) as usize];
    (lc_max + 1) as usize];
    (na_max + 1) as usize];
    (ma_max + 1) as usize];
    (la_max + 1) as usize
];
```

For d-orbitals (l=2), this is ~3×3×3×3×3×3×5 = ~3,645 f64 values = ~29 KB.

### HRR Recursion

```rust
fn hrr_recursive(vrr: &[...], lmn_a, lmn_b, lmn_c, lmn_d, ...) -> f64 {
    // Transfer b → a
    if lb > 0 {
        return hrr(..., (la+1,...), (lb-1,...), ...)
             + (a.x - b.x) * hrr(..., lmn_a, (lb-1,...), ...);
    }
    // ... similar for mb, nb, ld, md, nd ...

    // Base case: lookup in VRR
    vrr[la][ma][na][lc][mc][nc][0]
}
```

Recursion depth is O(la+lb+lc+ld), not O((la+lb+lc+ld)^n) as in original.

## Lessons Learned

### 1. HashMap is Expensive
- 7-tuple keys require hashing and equality checks
- Pointer chasing for collision resolution
- Much slower than direct array indexing

### 2. Memoization Strategy Matters
- **Bad**: Per-call local cache (original HGP)
- **Good**: Pre-compute once, reuse (optimized HGP)
- **Best**: Avoid recursion entirely (Standard THO)

### 3. Rust vs Julia Differences
Julia implementation insights:
- Pre-allocated arrays work well in both languages
- Julia's `@inbounds` ≈ Rust's Release mode bounds checking removal
- Similar optimization strategies apply across languages

### 4. When "Fancy" Algorithms Lose
HGP is mathematically elegant (HRR/VRR recursion), but:
- Standard THO's direct B-array is simpler and faster
- Complex algorithms need more optimization to compete
- "Simple and correct" often beats "clever but slow"

## Future Work

### Further HGP Optimizations (Not Recommended)
Possible but diminishing returns:
- [ ] Iterative HRR (eliminate recursion completely)
- [ ] Specialized code generation per angular momentum
- [ ] SIMD vectorization of VRR loops
- [ ] Better memory layout (SoA instead of AoS)

**Verdict**: Not worth the effort. Standard and Rys are already fast enough.

### Better Optimization Targets
Focus on:
- [ ] Full Rys coefficient tables (1500+ lines)
- [ ] SIMD for Standard THO B-array computation
- [ ] Integral screening (skip near-zero integrals)
- [ ] Batch processing optimizations

## Conclusion

**Mission Accomplished**: HGP optimization was successful, achieving 2-8.5x speedups.

**Practical Impact**: HGP is now usable, but Standard THO remains the best default choice for most cases, with Rys being optimal for high angular momentum.

**Key Takeaway**: The optimization demonstrated that careful attention to memory layout and computation order can yield massive performance improvements, even when the mathematical algorithm stays the same.

## Code Locations

- Original: `src/hgp.rs` (kept for comparison)
- Optimized: `src/hgp_opt.rs` ✅
- Julia reference: `https://github.com/rpmuller/MolecularIntegrals.jl/blob/master/src/HGP.jl`
- Benchmarks: `examples/hgp_comparison.rs`, `examples/benchmark.rs`
- Tests: 42/42 passing (original + optimized)
