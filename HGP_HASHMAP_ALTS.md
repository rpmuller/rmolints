# HashMap Alternatives for HGP VRR Storage

## Executive Summary

The Head-Gordon-Pople (HGP) method computes two-electron repulsion integrals using Vertical Recursion Relations (VRR) that require storing intermediate values indexed by 7 quantum numbers: `(i, j, k, ic, jc, kc, im)`.

**Current implementations:**
- **Original HGP** (`src/hgp.rs`): `HashMap<(i32, i32, i32, i32, i32, i32, i32), f64>` - **4x slower than Standard**
- **HGP Optimized** (`src/hgp_opt.rs`): 7D nested `Vec` - **1.4x slower than Standard**

**Key finding**: Switching from HashMap to nested Vec provides 2.8x speedup, but other alternatives may offer further improvements.

## Problem Statement

The VRR tensor has 7 dimensions with bounded ranges:
- `i, j, k`: Angular momentum components for center A (0 to la/ma/na, typically 0-2)
- `ic, jc, kc`: Angular momentum components for center C (0 to lc/mc/nc, typically 0-2)
- `im`: Auxiliary index (0 to mtot, where mtot = la+ma+na+lc+mc+nc, typically 0-8)

For benzene (p-orbitals, L=1): typical dimensions are roughly 2×2×2×2×2×2×5 = 1,280 elements

**Storage requirements:**
- HashMap: ~50-100 bytes per entry (key + value + overhead)
- Vec: 8 bytes per f64 element
- Total for nested Vec: 1,280 × 8 = 10 KB

## Alternative 1: Flat Array with Computed Index (Row-Major)

### Concept

Instead of 7D nested Vec, use a single 1D Vec and compute linear index:

```rust
struct VRRTensor {
    data: Vec<f64>,
    dims: [usize; 7], // [la_max, ma_max, na_max, lc_max, mc_max, nc_max, mtot_max]
}

impl VRRTensor {
    fn index(&self, i: usize, j: usize, k: usize,
             ic: usize, jc: usize, kc: usize, im: usize) -> usize {
        // Row-major ordering (rightmost index varies fastest)
        im
        + kc * self.dims[6]
        + jc * self.dims[6] * self.dims[5]
        + ic * self.dims[6] * self.dims[5] * self.dims[4]
        + k  * self.dims[6] * self.dims[5] * self.dims[4] * self.dims[3]
        + j  * self.dims[6] * self.dims[5] * self.dims[4] * self.dims[3] * self.dims[2]
        + i  * self.dims[6] * self.dims[5] * self.dims[4] * self.dims[3] * self.dims[2] * self.dims[1]
    }

    #[inline]
    fn get(&self, i: usize, j: usize, k: usize,
           ic: usize, jc: usize, kc: usize, im: usize) -> f64 {
        self.data[self.index(i, j, k, ic, jc, kc, im)]
    }

    #[inline]
    fn set(&mut self, i: usize, j: usize, k: usize,
           ic: usize, jc: usize, kc: usize, im: usize, val: f64) {
        let idx = self.index(i, j, k, ic, jc, kc, im);
        self.data[idx] = val;
    }
}
```

### Advantages
- ✅ **Better cache locality**: Contiguous memory layout
- ✅ **Faster index computation**: 7 multiplications + 6 additions vs 7 pointer dereferences
- ✅ **Smaller memory footprint**: No Vec overhead for inner dimensions
- ✅ **SIMD potential**: Contiguous data enables vectorization

### Disadvantages
- ❌ **Complex index calculation**: More arithmetic per access
- ❌ **Less readable**: Harder to understand access pattern
- ❌ **Must pre-compute strides**: Need to store dimension products

### Optimization: Pre-compute Strides

```rust
struct VRRTensor {
    data: Vec<f64>,
    strides: [usize; 7], // Pre-computed products
}

impl VRRTensor {
    fn new(dims: [usize; 7]) -> Self {
        let total_size = dims.iter().product();
        let mut strides = [1; 7];

        // Compute strides: stride[i] = product of all dims to the right
        for i in (0..6).rev() {
            strides[i] = strides[i + 1] * dims[i + 1];
        }

        VRRTensor {
            data: vec![0.0; total_size],
            strides,
        }
    }

    #[inline(always)]
    fn index(&self, i: usize, j: usize, k: usize,
             ic: usize, jc: usize, kc: usize, im: usize) -> usize {
        i * self.strides[0]
        + j * self.strides[1]
        + k * self.strides[2]
        + ic * self.strides[3]
        + jc * self.strides[4]
        + kc * self.strides[5]
        + im  // stride[6] is always 1
    }
}
```

**Performance**: 6 multiplications + 6 additions = 12 operations
Compare to nested Vec: 7 pointer dereferences + 7 bounds checks

**Expected speedup**: 1.2-1.5x over nested Vec (if inlining works well)

## Alternative 2: Hybrid HashMap + Array

### Concept

Use HashMap for outer dimensions (rarely accessed), flat array for inner dimensions (frequently accessed):

```rust
struct HybridVRR {
    // Outer keys: (i, j, k, ic, jc, kc)
    // Inner value: Vec<f64> for all im values
    cache: HashMap<(i32, i32, i32, i32, i32, i32), Vec<f64>>,
    mtot_max: usize,
}

impl HybridVRR {
    fn get(&self, i: i32, j: i32, k: i32,
           ic: i32, jc: i32, kc: i32, im: usize) -> f64 {
        self.cache[&(i, j, k, ic, jc, kc)][im]
    }

    fn set(&mut self, i: i32, j: i32, k: i32,
           ic: i32, jc: i32, kc: i32, im: usize, val: f64) {
        self.cache.entry((i, j, k, ic, jc, kc))
            .or_insert_with(|| vec![0.0; self.mtot_max + 1])
            [im] = val;
    }
}
```

### Advantages
- ✅ **Sparse storage**: Only allocates (i,j,k,ic,jc,kc) tuples that are actually used
- ✅ **Fast inner loop**: Direct array access for im dimension
- ✅ **Simple to implement**: Minimal changes to existing code

### Disadvantages
- ❌ **Still has HashMap overhead**: 6-tuple keys are expensive
- ❌ **Unpredictable allocation**: Hard to reason about memory usage
- ❌ **Worse cache locality**: Scattered allocations

### Performance Estimate

For benzene (2×2×2×2×2×2 = 64 outer keys, ~5 im values each):
- HashMap lookups: ~64 unique lookups per integral
- Array accesses: ~320 accesses (64 × 5)
- **Expected**: 1.5-2x faster than full HashMap, 1.1x slower than flat array

## Alternative 3: Sparse Tensor with Computed Bounds

### Concept

Since not all combinations are valid (constraint: `im ≤ mtot - i - j - k - ic - jc - kc`), we can pack data more efficiently:

```rust
struct SparseTensor {
    data: Vec<f64>,
    // Store cumulative offsets for each (i,j,k,ic,jc,kc)
    offsets: HashMap<(i32, i32, i32, i32, i32, i32), usize>,
    mtot: i32,
}

impl SparseTensor {
    fn im_max(&self, i: i32, j: i32, k: i32,
              ic: i32, jc: i32, kc: i32) -> i32 {
        self.mtot - i - j - k - ic - jc - kc
    }

    fn get(&self, i: i32, j: i32, k: i32,
           ic: i32, jc: i32, kc: i32, im: i32) -> f64 {
        let offset = self.offsets[&(i, j, k, ic, jc, kc)];
        self.data[offset + im as usize]
    }
}
```

### Memory Savings

Example for benzene (mtot ≈ 4):
- Full tensor: 2×2×2×2×2×2×5 = 1,280 elements
- Sparse tensor: Sum over all valid (i,j,k,ic,jc,kc) of (mtot - sum + 1)
  - Approximately: 64 × 2.5 ≈ 160 elements

**Memory reduction**: 8x smaller!

### Disadvantages
- ❌ **Still uses HashMap**: Offset lookups are slow
- ❌ **Complex indexing**: Hard to compute offsets incrementally
- ❌ **Poor cache locality**: Scattered data

## Alternative 4: Column-Major Ordering

### Concept

Reorder indices so most frequently accessed dimension varies slowest:

```rust
// Current row-major: (i, j, k, ic, jc, kc, im)
// Hottest access pattern in VRR: incrementing im
// Row-major: im varies fastest → good!

// Alternative column-major: (im, kc, jc, ic, k, j, i)
fn index_colmajor(&self, i: usize, j: usize, k: usize,
                  ic: usize, jc: usize, kc: usize, im: usize) -> usize {
    i
    + j  * self.dims[0]
    + k  * self.dims[0] * self.dims[1]
    + ic * self.dims[0] * self.dims[1] * self.dims[2]
    + jc * self.dims[0] * self.dims[1] * self.dims[2] * self.dims[3]
    + kc * self.dims[0] * self.dims[1] * self.dims[2] * self.dims[3] * self.dims[4]
    + im * self.dims[0] * self.dims[1] * self.dims[2] * self.dims[3] * self.dims[4] * self.dims[5]
}
```

### Analysis

Looking at VRR recursion patterns:
1. **Base case loop**: `for im in 0..=mtot` - accesses consecutive im values ✅
2. **Build-up loops**: Typically increment `i`, `j`, `k`, `ic`, `jc`, `kc` one at a time
3. **Dependency**: `vrr[i][j][k][ic][jc][kc][im]` depends on:
   - Same indices with `im+1` (consecutive)
   - Decremented outer index with same `im` (strided)

**Conclusion**: Row-major (im fastest) is optimal for VRR access patterns!

## Alternative 5: Unsafe Direct Indexing

### Concept

Remove bounds checks using unsafe code:

```rust
impl VRRTensor {
    #[inline(always)]
    unsafe fn get_unchecked(&self, i: usize, j: usize, k: usize,
                            ic: usize, jc: usize, kc: usize, im: usize) -> f64 {
        let idx = self.index(i, j, k, ic, jc, kc, im);
        *self.data.get_unchecked(idx)
    }

    #[inline(always)]
    unsafe fn set_unchecked(&mut self, i: usize, j: usize, k: usize,
                            ic: usize, jc: usize, kc: usize, im: usize, val: f64) {
        let idx = self.index(i, j, k, ic, jc, kc, im);
        *self.data.get_unchecked_mut(idx) = val;
    }
}
```

### Safety

Since VRR construction iterates through indices in a controlled manner (nested loops with known bounds), unsafe access is safe IF:
- Tensor dimensions are set correctly at construction
- Index computation is correct
- No out-of-bounds access in user code

### Expected Speedup

Removing bounds checks: 1.05-1.15x speedup (bounds checks are cheap on modern CPUs)

**Risk vs reward**: Not worth it - minimal gain, introduces unsafety

## Performance Comparison Matrix

| Method | Memory (bytes) | Access Time | Complexity | Cache Locality | Sparsity |
|--------|---------------|-------------|------------|----------------|----------|
| **HashMap (original)** | 50-100 per entry | 50-100 ns | Low | Poor ❌ | Excellent ✅ |
| **7D Nested Vec (current)** | 8 per entry + overhead | 10-15 ns | Low | Medium | None ❌ |
| **Flat Array + Strides** | 8 per entry | 5-8 ns | Medium | Excellent ✅ | None ❌ |
| **Hybrid HashMap+Array** | 20-30 per entry | 20-30 ns | Medium | Medium | Good |
| **Sparse Tensor** | 8 per valid entry | 40-60 ns | High | Medium | Excellent ✅ |
| **Column-Major** | 8 per entry | 5-8 ns | Medium | Good | None ❌ |
| **Unsafe Unchecked** | 8 per entry | 4-6 ns | Medium | Excellent ✅ | None ❌ |

## Actual Benchmark Results

**Test**: `examples/vrr_tensor_comparison.rs` - measures pure tensor access patterns

| Test Case | Nested Vec | Flat Array | HashMap | Flat Speedup |
|-----------|-----------|------------|---------|--------------|
| s-orbitals (2×2×2×2×2×2×3) | 0.007s | 0.002s | 0.083s | **2.94x faster** ✅ |
| p-orbitals (3×3×3×3×3×3×5) | 0.041s | 0.015s | 0.821s | **2.79x faster** ✅ |
| d-orbitals (4×4×4×4×4×4×9) | 0.094s | 0.028s | 1.914s | **3.40x faster** ✅ |

**Key findings:**
1. ✅ **Flat array is 2.8-3.4x faster** - exceeds predictions!
2. ❌ **HashMap is 12-20x slower** - confirms it's terrible for this use case
3. 📈 **Speedup increases with size** - cache locality matters more for larger tensors

### Predictions for Real HGP Implementation

The benchmark measures pure access patterns. Real HGP has computation between accesses, so actual speedup will be lower.

**Conservative estimate**: 30-50% of benchmark speedup = 1.4-2.0x faster in practice

For benzene (36 basis functions, 111,055 ERIs):

| Method | Predicted Time | vs Current Opt | vs Standard | vs Rys |
|--------|---------------|----------------|-------------|---------|
| **Current (7D Vec)** | 1994 ms | 1.00x | 1.42x slower | 2.83x slower |
| **Flat Array (conservative)** | 1000-1425 ms | **1.4-2.0x faster** ✅ | 0.69-0.99x | 1.42-2.02x slower |
| **Flat Array (optimistic)** | 665 ms | **3.0x faster** ✅ | **0.46x** ✅ | 1.06x slower |

**Best case scenario**: If we achieve 3x speedup, HGP-Opt could beat Standard THO!
**Realistic scenario**: 1.5x speedup would make HGP-Opt comparable to Standard (~1330ms)

## Recommendation

### ✅ STRONGLY RECOMMEND: Implement Flat Array with Pre-computed Strides

**Why:**
1. **HUGE speedup potential**: Benchmarks show 2.8-3.4x faster tensor access
2. **Conservative estimate**: 1.4-2.0x speedup in real HGP (accounting for computation time)
3. **Optimistic possibility**: Could beat Standard THO if full speedup translates!
4. **Better cache locality**: Contiguous memory layout proven effective
5. **Safe**: No unsafe code needed
6. **Reasonable complexity**: Stride computation is straightforward

**Implementation plan:**

1. Create `VRRTensor` struct in `src/hgp_opt.rs` (see proof-of-concept in `examples/vrr_tensor_comparison.rs`)
2. Replace 7D nested Vec with flat Vec + strides
3. Update all `vrr[i][j][k][ic][jc][kc][im]` to `vrr.get(i,j,k,ic,jc,kc,im)`
4. Benchmark on micro-benchmarks (s, p, d orbitals)
5. Benchmark on real molecules (H2O, benzene)
6. If successful (>30% speedup), make it default HGP-Opt implementation
7. Consider backporting to original HGP (though HashMap replacement is more complex)

**Expected outcome:**

Conservative (1.5x speedup):
- HGP-Opt on benzene: ~1330 ms (from 1994 ms)
- Competitive with Standard THO: 1693 ms
- Still slower than Rys: 1408 ms
- But respectable performance!

Optimistic (2.5x speedup):
- HGP-Opt on benzene: ~800 ms (from 1994 ms)
- **BEATS Standard THO** by 2.1x! ✅
- Competitive with Rys: 1408 ms
- Makes HGP a viable production method!

**Reality check**: Even if we only get 30% speedup (1.3x), that's ~1530ms - still a meaningful improvement and worth implementing for completeness.

### Do NOT implement:

1. **Unsafe unchecked access**: Risk not worth 5-10% gain
2. **Sparse tensor**: Complexity outweighs memory savings (only 8x for 10KB → 1.25KB)
3. **Hybrid HashMap**: Worst of both worlds
4. **Column-major**: Row-major is already optimal for VRR access patterns

## Conclusion

**Major Discovery**: Flat array with pre-computed strides is **2.8-3.4x faster** than nested Vec for tensor access patterns typical of VRR recursions!

### Performance Evolution of HGP

1. **Original HGP (HashMap)**: 5583 ms on benzene (baseline)
2. **HGP-Opt (Nested Vec)**: 1994 ms - **2.8x faster** ✅
3. **Predicted with Flat Array**: 800-1330 ms - **additional 1.5-2.5x faster** ✅
4. **Total improvement**: 4.2-7.0x faster than original!

### Why This Matters

The nested Vec → flat array optimization could:
- Make HGP-Opt competitive with or even beat Standard THO (1693 ms)
- Bring HGP within striking distance of Rys Quadrature (1408 ms)
- Transform HGP from "academic curiosity" to "viable alternative"
- Demonstrate that algorithmic complexity can be overcome with careful data structure design

### Reality Check

Even if real-world speedup is only half of benchmark results (1.4x instead of 2.8x):
- HGP-Opt would drop to ~1425 ms
- Still beats current HGP-Opt by 40%
- Competitive with Standard THO
- Worth implementing!

**Bottom line**: This optimization has the potential to dramatically change the HGP performance story. It's no longer just an academic exercise - it could become a production-worthy method.

## Next Steps

### High Priority (Do This!)
1. ✅ **Implement `VRRTensor` struct** in `src/hgp_opt.rs` using flat array approach
2. ✅ **Benchmark on micro-benchmarks** (s, p, d orbitals) - expect 1.5-2.5x speedup
3. ✅ **Benchmark on real molecules** (H2O, benzene) - this is the real test!
4. ✅ **Update BENCHMARK_RESULTS.md** with new timings
5. ✅ **If >30% speedup**: Make it default HGP-Opt implementation

### Medium Priority (If Above Succeeds)
6. **Backport to original HGP** - replace HashMap with flat array
7. **Profile to find next bottleneck** - HRR might become limiting factor
8. **Consider SIMD optimizations** - flat array enables vectorization

### Low Priority (Academic Interest)
9. **Try unsafe unchecked access** - measure if bounds checks matter
10. **Experiment with column-major** - though row-major is likely optimal
11. **Write research note** - document this optimization for publication

## Key Takeaway

**"The difference between theory and practice is larger in practice than in theory."**

We predicted 1.2x speedup. Benchmarks show 3x speedup. The real-world result will be somewhere in between, but the proof-of-concept demonstrates this optimization is **absolutely worth pursuing**.

The flat array approach transforms HGP from a slow, pedagogical implementation to a potentially competitive production method. This is a textbook example of how data structure choice can make or break algorithmic performance.

**Rys Quadrature is still the king**, but with flat array optimization, HGP-Opt could earn a legitimate place as "the method that's almost as fast and easier to understand."
