use rmolints::basis::build_sto3g_basis;
use rmolints::molecule::Molecule;
use rmolints::parallel::{compute_eri_tensor_parallel, ERIMethod};
use std::time::Instant;

fn main() {
    println!("\n=== HGP-Opt Profiling ===");
    println!("Generating CPU profiling data and flamegraph for HGP-Opt implementation\n");

    // Use benzene as the test molecule - large enough to show bottlenecks
    // but small enough to profile in reasonable time
    let molecule = Molecule::benzene();
    let basis = build_sto3g_basis(&molecule);
    let n = basis.len();
    let n_unique = n * (n + 1) * (n * (n + 1) / 2 + 1) / 8;

    println!("Molecule: Benzene (C6H6)");
    println!("Basis functions: {}", n);
    println!("Unique ERIs: {}\n", n_unique);

    // Warmup run to ensure code is compiled and caches are warm
    println!("Warmup run...");
    let _ = compute_eri_tensor_parallel(&basis, ERIMethod::HeadGordonPopleOpt);

    println!("Starting profiling run...\n");

    // Start profiler
    let guard = pprof::ProfilerGuardBuilder::default()
        .frequency(1000) // Sample at 1000 Hz
        .blocklist(&["libc", "libgcc", "pthread", "vdso"]) // Exclude system libraries
        .build()
        .unwrap();

    // Run HGP-Opt multiple times to get enough samples
    let iterations = 10;
    let start = Instant::now();

    for i in 0..iterations {
        let iter_start = Instant::now();
        let eris = compute_eri_tensor_parallel(&basis, ERIMethod::HeadGordonPopleOpt);
        let iter_elapsed = iter_start.elapsed().as_secs_f64() * 1000.0;

        println!("Iteration {}/{}: {:.2} ms ({} ERIs computed)",
                 i + 1, iterations, iter_elapsed, eris.len());
    }

    let total_elapsed = start.elapsed().as_secs_f64() * 1000.0;
    let avg_elapsed = total_elapsed / iterations as f64;

    println!("\nTotal time: {:.2} ms", total_elapsed);
    println!("Average per iteration: {:.2} ms\n", avg_elapsed);

    // Stop profiler and generate flamegraph
    println!("Generating flamegraph...");

    if let Ok(report) = guard.report().build() {
        // Write flamegraph to file
        let file = std::fs::File::create("flamegraph.svg").unwrap();
        report.flamegraph(file).unwrap();
        println!("✅ Flamegraph saved to: flamegraph.svg");

        // Print top functions by time
        println!("\n=== Top Functions by Sample Count ===\n");

        // Get the profiling data
        let mut samples: Vec<(String, usize)> = Vec::new();

        // Collect samples from the report
        for (name, sample_count) in report.data.iter() {
            let name_str = format!("{:?}", name);
            // Filter for rmolints functions only
            if name_str.contains("rmolints") || name_str.contains("hgp") {
                samples.push((name_str, (*sample_count).max(0) as usize));
            }
        }

        // Sort by sample count (descending)
        samples.sort_by(|a, b| b.1.cmp(&a.1));

        // Print top 20
        for (i, (name, count)) in samples.iter().take(20).enumerate() {
            println!("{}. {} - {} samples", i + 1, name, count);
        }
    } else {
        println!("❌ Failed to generate profiling report");
    }

    println!("\n=== Profiling Complete ===");
    println!("\nTo view the flamegraph:");
    println!("  1. Open flamegraph.svg in a web browser");
    println!("  2. Click on sections to zoom in");
    println!("  3. Look for wide sections (high CPU time)");
    println!("\nExpected bottlenecks:");
    println!("  - compute_vrr_tensor (should dominate)");
    println!("  - VRR stages 6-7 (innermost loops)");
    println!("  - Inner im loops with FMA operations");
}
