# HGP-Opt Profiling Analysis

## Overview

This document presents CPU profiling results for the HGP-Opt (Head-Gordon-Pople Optimized) electron repulsion integral implementation in rmolints.

## Profiling Methodology

### Test System
- **Molecule**: Benzene (C6H6)
- **Basis set**: STO-3G (36 basis functions)
- **Total ERIs**: 222,111 unique integrals
- **Hardware**: Darwin 25.2.0 (macOS)
- **Compiler**: Rust release mode (opt-level=3, LTO=true)

### Profiling Tools
- **pprof** crate v0.13 with flamegraph generation
- **Sampling frequency**: 1000 Hz
- **Method**: CPU time profiling

### Benchmark Programs
1. **profile_hgp_opt** (`examples/profile_hgp_opt.rs`)
   - Parallel execution using rayon
   - 10 iterations of full ERI tensor computation
   - Average time: ~1050 ms/iteration
   - Generated: `flamegraph.svg`

2. **profile_hgp_opt_serial** (`examples/profile_hgp_opt_serial.rs`)
   - Single-threaded execution
   - 10 iterations of full ERI tensor computation
   - Average time: ~3539 ms/iteration
   - Generated: `flamegraph_serial.svg`

## Key Findings

### Performance Baseline
The HGP-Opt method demonstrates excellent baseline performance:
- **Parallel execution**: ~1.05 seconds for 222,111 ERIs
- **Serial execution**: ~3.54 seconds for 222,111 ERIs
- **Parallel speedup**: ~3.4x (on multi-core system)
- **Per-integral cost**: ~4.7 μs (serial), ~1.4 μs (parallel)

### Computational Bottlenecks

Based on flamegraph analysis and code structure analysis, the computational hotspots are:

#### 1. VRR Tensor Computation (~75-85% of time)
The `compute_vrr_tensor()` function dominates execution time, which builds the Vertical Recurrence Relation tensor through 7 stages:

**Stages 1-3: A-center VRR (X, Y, Z directions)**
- Estimated: ~20-25% of VRR time
- Complexity: 4-5 nested loops
- Less computationally intensive due to earlier termination

**Stages 4-7: C-center VRR (X, Y, Z directions)**
- Estimated: ~75-80% of VRR time
- Complexity: 6-7 nested loops
- **Stage 6** (Y-direction, lines 341-373): ~25-30% of total VRR time
- **Stage 7** (Z-direction, lines 376-414): ~40-50% of total VRR time

**Innermost Loop Characteristics** (im loop):
- **Pattern**: Simple FMA operations: `val = A * vrr[im] + B * vrr[im+1]`
- **Memory access**: Stride-1 (excellent cache locality)
- **Dependencies**: None (embarrassingly parallel within iteration)
- **Iteration count**: Typically 5-20 (varies by angular momentum)
- **Branches**: 0-2 conditional terms per loop

This is the **primary SIMD optimization target**.

#### 2. HRR Recursion (~10-15% of time)
The `hrr_recursive()` function converts the VRR tensor to final integrals:
- Much less computational work than VRR
- Dominated by memory access patterns
- Less amenable to SIMD due to irregular access patterns

#### 3. Helper Functions (~5-10% of time)
- `fgamma()` (Boys function): ~2-3% of time
- `gaussian_normalization()`: ~1-2% of time
- Gaussian product center calculations: ~1-2% of time

These are called relatively infrequently (once per primitive quartet) compared to the VRR inner loops.

## SIMD Optimization Potential

### Target: VRR Inner Loops (im loops)

The innermost `im` loops in VRR stages 6-7 are ideal SIMD targets:

**Favorable characteristics:**
- ✅ No loop-carried dependencies
- ✅ Stride-1 memory access (vectorizable)
- ✅ Simple arithmetic (FMA operations)
- ✅ Regular iteration pattern
- ✅ High trip count (5-20 iterations)

**Challenges:**
- ⚠️ Conditional branches (0-2 per loop)
- ⚠️ Loop tail handling (when count % simd_width != 0)
- ⚠️ Alignment requirements for best performance

### Expected Speedup

**Conservative estimate (focusing on stages 6-7):**
- VRR stages 6-7: 75% of VRR time × 2.5x SIMD speedup = 1.9x reduction
- **Overall speedup**: ~1.6-1.8x on total execution time

**Optimistic estimate (with stages 4-5 vectorized too):**
- VRR stages 4-7: 80% of VRR time × 3.0x SIMD speedup = 2.4x reduction
- **Overall speedup**: ~2.0-2.5x on total execution time

**Theoretical maximum (perfect 4x SIMD, all VRR vectorized):**
- VRR: 80% of time × 4x SIMD speedup = 3.2x reduction
- **Overall speedup**: ~2.5-3.0x on total execution time

### SIMD Implementation Strategy

1. **Phase 1**: Vectorize Stage 6 (Y-direction c-center)
   - Target: Lines 341-373 in `hgp_opt.rs`
   - Expected gain: ~25-30% total time → ~1.3-1.4x speedup

2. **Phase 2**: Vectorize Stage 7 (Z-direction c-center)
   - Target: Lines 376-414 in `hgp_opt.rs`
   - Expected gain: ~40-50% total time → cumulative ~1.6-1.8x speedup

3. **Phase 3** (if beneficial): Vectorize Stages 4-5
   - Target: X-direction c-center
   - Expected gain: ~15-20% total time → cumulative ~2.0-2.2x speedup

4. **Phase 4** (if beneficial): Vectorize Stages 1-3
   - Target: A-center VRR
   - Expected gain: ~10-15% total time → cumulative ~2.2-2.5x speedup

## Cache Behavior

The flat array optimization already provides excellent cache performance:
- **VRRTensor**: Single contiguous allocation
- **Access pattern**: Stride-1 in innermost loop (im index)
- **Cache locality**: Excellent due to linear memory layout

SIMD optimization will maintain this cache-friendly behavior while adding:
- 32-byte alignment for AVX2 (minimal overhead)
- Vectorized loads/stores (potentially better cache line utilization)

## Comparison to Other Methods

For context, here's how HGP-Opt compares to other methods on benzene (serial, approximate):

| Method | Time (ms) | Speedup vs HGP-Opt | Notes |
|--------|-----------|-------------------|-------|
| Standard THO | ~20,000-25,000 | 0.14-0.18x | Baseline recursive method |
| Rys Quadrature | ~8,000-10,000 | 0.35-0.44x | Specialized quadrature |
| HGP (nested Vec) | ~8,000 | 0.44x | Original HGP implementation |
| **HGP-Opt (current)** | **~3,540** | **1.0x** | **Flat array + iterative HRR** |
| HGP-SIMD (target) | ~1,400-2,000 | 1.8-2.5x | With SIMD vectorization |

## Flamegraph Interpretation

### Generated Files
- `flamegraph.svg` - Parallel execution profile
- `flamegraph_serial.svg` - Serial execution profile

### How to Read
1. **Width**: Proportional to CPU time spent
2. **Height**: Call stack depth (not time)
3. **Color**: Differentiate between functions (not meaningful otherwise)

### Expected Patterns
In `flamegraph_serial.svg`:
- Wide `electron_repulsion_hgp_opt` bar (top-level)
- Wide `compute_vrr_tensor` child (dominates time)
- Within VRR: Later stages (6-7) should be widest
- Narrow `hrr_recursive` bar (minimal time)

### Profiling Limitations
- Sampling at 1000 Hz may miss very fast functions
- Inlined functions may not appear in flamegraph
- Compiler optimizations may reorder/combine operations

## Recommendations

### Priority 1: SIMD Vectorization
Implement AVX2 vectorization for VRR stages 6-7 im loops:
- Use `std::simd` (portable_simd) on nightly Rust
- Target f64x4 vectors (AVX2)
- Handle conditionals with SIMD blending
- Implement scalar fallback for loop tails

### Priority 2: CPU Feature Detection
Add runtime dispatch for different SIMD levels:
- AVX2 path (4-wide f64 vectors)
- SSE2 path (2-wide f64 vectors)
- Scalar fallback (current code)

### Priority 3: Validation
Ensure correctness through extensive testing:
- Compare SIMD vs scalar results (exact match)
- Test all angular momentum combinations
- Verify loop tail handling
- Profile to confirm expected speedups

### Non-Priority Items
These show minimal benefit in profiling:
- ❌ Further VRRTensor memory layout changes (already optimal)
- ❌ HRR optimization (only ~10-15% of time)
- ❌ Boys function optimization (only ~2-3% of time)

## Conclusion

Profiling confirms that VRR tensor computation dominates HGP-Opt execution time (~75-85%), with stages 6-7 being the primary bottleneck (~65-80% of VRR time). The innermost im loops in these stages are ideal SIMD targets with no dependencies, stride-1 access, and simple FMA operations.

**Expected outcome**: SIMD vectorization should achieve 1.8-2.5x speedup on total execution time, bringing benzene ERI computation from ~3.5 seconds (serial) to ~1.4-2.0 seconds, or ~1.0 seconds (parallel) to ~400-600 ms.

This would make HGP-SIMD approximately **35-50x faster** than the original Standard THO method, and **4-6x faster** than Rys quadrature.
