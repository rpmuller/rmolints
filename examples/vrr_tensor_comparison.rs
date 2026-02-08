/// Proof-of-concept: Compare different VRR tensor storage strategies
///
/// This example demonstrates the performance difference between:
/// 1. 7D nested Vec (current HGP-Opt approach)
/// 2. Flat array with pre-computed strides
/// 3. HashMap (original HGP approach)

use std::collections::HashMap;
use std::time::Instant;

// ============================================================================
// Approach 1: 7D Nested Vec (Current)
// ============================================================================

type NestedVec7D = Vec<Vec<Vec<Vec<Vec<Vec<Vec<f64>>>>>>>;

fn create_nested_vec(dims: [usize; 7]) -> NestedVec7D {
    vec![
        vec![
            vec![
                vec![
                    vec![
                        vec![
                            vec![0.0; dims[6]];
                        dims[5]];
                    dims[4]];
                dims[3]];
            dims[2]];
        dims[1]];
    dims[0]]
}

#[inline(always)]
fn nested_get(arr: &NestedVec7D, i: usize, j: usize, k: usize,
              ic: usize, jc: usize, kc: usize, im: usize) -> f64 {
    arr[i][j][k][ic][jc][kc][im]
}

#[inline(always)]
fn nested_set(arr: &mut NestedVec7D, i: usize, j: usize, k: usize,
              ic: usize, jc: usize, kc: usize, im: usize, val: f64) {
    arr[i][j][k][ic][jc][kc][im] = val;
}

// ============================================================================
// Approach 2: Flat Array with Pre-computed Strides
// ============================================================================

struct FlatTensor {
    data: Vec<f64>,
    strides: [usize; 7],
    dims: [usize; 7],
}

impl FlatTensor {
    fn new(dims: [usize; 7]) -> Self {
        let total_size: usize = dims.iter().product();
        let mut strides = [1; 7];

        // Compute strides: stride[i] = product of all dims to the right
        for i in (0..6).rev() {
            strides[i] = strides[i + 1] * dims[i + 1];
        }

        println!("  Flat tensor dimensions: {:?}", dims);
        println!("  Computed strides: {:?}", strides);
        println!("  Total size: {} elements ({:.2} KB)", total_size, total_size * 8 / 1024);

        FlatTensor {
            data: vec![0.0; total_size],
            strides,
            dims,
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

    #[inline(always)]
    fn get(&self, i: usize, j: usize, k: usize,
           ic: usize, jc: usize, kc: usize, im: usize) -> f64 {
        self.data[self.index(i, j, k, ic, jc, kc, im)]
    }

    #[inline(always)]
    fn set(&mut self, i: usize, j: usize, k: usize,
           ic: usize, jc: usize, kc: usize, im: usize, val: f64) {
        let idx = self.index(i, j, k, ic, jc, kc, im);
        self.data[idx] = val;
    }
}

// ============================================================================
// Approach 3: HashMap (Original)
// ============================================================================

type HashMapTensor = HashMap<(usize, usize, usize, usize, usize, usize, usize), f64>;

fn hashmap_get(map: &HashMapTensor, i: usize, j: usize, k: usize,
               ic: usize, jc: usize, kc: usize, im: usize) -> f64 {
    *map.get(&(i, j, k, ic, jc, kc, im)).unwrap_or(&0.0)
}

fn hashmap_set(map: &mut HashMapTensor, i: usize, j: usize, k: usize,
               ic: usize, jc: usize, kc: usize, im: usize, val: f64) {
    map.insert((i, j, k, ic, jc, kc, im), val);
}

// ============================================================================
// Benchmark Functions
// ============================================================================

fn benchmark_nested_vec(dims: [usize; 7], iterations: usize) -> f64 {
    let mut arr = create_nested_vec(dims);
    let start = Instant::now();

    for _ in 0..iterations {
        // Simulate VRR-like access pattern
        for i in 0..dims[0] {
            for j in 0..dims[1] {
                for k in 0..dims[2] {
                    for ic in 0..dims[3] {
                        for jc in 0..dims[4] {
                            for kc in 0..dims[5] {
                                for im in 0..dims[6].saturating_sub(1) {
                                    // Read two values (im and im+1)
                                    let v1 = nested_get(&arr, i, j, k, ic, jc, kc, im);
                                    let v2 = nested_get(&arr, i, j, k, ic, jc, kc, im + 1);
                                    // Write result
                                    nested_set(&mut arr, i, j, k, ic, jc, kc, im, v1 + v2);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    start.elapsed().as_secs_f64()
}

fn benchmark_flat_tensor(dims: [usize; 7], iterations: usize) -> f64 {
    let mut tensor = FlatTensor::new(dims);
    let start = Instant::now();

    for _ in 0..iterations {
        for i in 0..dims[0] {
            for j in 0..dims[1] {
                for k in 0..dims[2] {
                    for ic in 0..dims[3] {
                        for jc in 0..dims[4] {
                            for kc in 0..dims[5] {
                                for im in 0..dims[6].saturating_sub(1) {
                                    let v1 = tensor.get(i, j, k, ic, jc, kc, im);
                                    let v2 = tensor.get(i, j, k, ic, jc, kc, im + 1);
                                    tensor.set(i, j, k, ic, jc, kc, im, v1 + v2);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    start.elapsed().as_secs_f64()
}

fn benchmark_hashmap(dims: [usize; 7], iterations: usize) -> f64 {
    let mut map: HashMapTensor = HashMap::new();
    let start = Instant::now();

    for _ in 0..iterations {
        for i in 0..dims[0] {
            for j in 0..dims[1] {
                for k in 0..dims[2] {
                    for ic in 0..dims[3] {
                        for jc in 0..dims[4] {
                            for kc in 0..dims[5] {
                                for im in 0..dims[6].saturating_sub(1) {
                                    let v1 = hashmap_get(&map, i, j, k, ic, jc, kc, im);
                                    let v2 = hashmap_get(&map, i, j, k, ic, jc, kc, im + 1);
                                    hashmap_set(&mut map, i, j, k, ic, jc, kc, im, v1 + v2);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    start.elapsed().as_secs_f64()
}

// ============================================================================
// Main
// ============================================================================

fn main() {
    println!("\n=== VRR Tensor Storage Comparison ===\n");

    let test_cases = vec![
        ("s-orbitals (small)", [2, 2, 2, 2, 2, 2, 3], 10000),
        ("p-orbitals (medium)", [3, 3, 3, 3, 3, 3, 5], 5000),
        ("d-orbitals (large)", [4, 4, 4, 4, 4, 4, 9], 1000),
    ];

    for (name, dims, iterations) in test_cases {
        println!("\n--- Test: {} ---", name);
        println!("Dimensions: {:?}", dims);
        println!("Iterations: {}", iterations);
        println!();

        // Warmup
        println!("Warming up...");
        let _ = benchmark_nested_vec(dims, 1);
        let _ = benchmark_flat_tensor(dims, 1);
        let _ = benchmark_hashmap(dims, 1);
        println!();

        // Benchmark
        println!("Benchmarking 7D Nested Vec...");
        let nested_time = benchmark_nested_vec(dims, iterations);
        println!("  Time: {:.3} seconds", nested_time);

        println!("\nBenchmarking Flat Array + Strides...");
        let flat_time = benchmark_flat_tensor(dims, iterations);
        println!("  Time: {:.3} seconds", flat_time);

        println!("\nBenchmarking HashMap...");
        let hashmap_time = benchmark_hashmap(dims, iterations);
        println!("  Time: {:.3} seconds", hashmap_time);

        // Results
        println!("\n{}", "=".repeat(50));
        println!("RESULTS:");
        println!("{}", "=".repeat(50));
        println!("Nested Vec:  {:.3}s  (baseline)", nested_time);
        println!("Flat Array:  {:.3}s  ({:.2}x {})",
                 flat_time,
                 nested_time / flat_time,
                 if flat_time < nested_time { "FASTER ✅" } else { "slower" });
        println!("HashMap:     {:.3}s  ({:.2}x {})",
                 hashmap_time,
                 nested_time / hashmap_time,
                 if hashmap_time > nested_time { "slower ❌" } else { "FASTER" });
        println!("{}", "=".repeat(50));

        let speedup = nested_time / flat_time;
        if speedup > 1.0 {
            println!("\n🎉 Flat array is {:.1}% faster than nested Vec!", (speedup - 1.0) * 100.0);
        } else {
            println!("\n⚠️  Flat array is {:.1}% slower than nested Vec", (1.0 - speedup) * 100.0);
        }
    }

    println!("\n\n=== Summary ===\n");
    println!("Expected results:");
    println!("  • HashMap: 2-4x slower (hash overhead, poor cache locality)");
    println!("  • Flat Array: 1.2-1.5x faster (better cache locality, no pointer chasing)");
    println!("  • Effect increases with tensor size (cache matters more)");
    println!("\nConclusion:");
    println!("  If flat array is consistently faster, it should replace nested Vec in HGP-Opt!");
}
