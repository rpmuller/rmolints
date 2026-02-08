use rmolints::basis::build_sto3g_basis;
use rmolints::hgp_opt::electron_repulsion_hgp_opt;
use rmolints::molecule::Molecule;
use std::time::Instant;

fn main() {
    println!("\n=== HGP-Opt Serial Profiling ===");
    println!("Profiling single-threaded HGP-Opt to identify bottlenecks\n");

    // Use benzene as test molecule
    let molecule = Molecule::benzene();
    let basis = build_sto3g_basis(&molecule);
    let n = basis.len();

    println!("Molecule: Benzene (C6H6)");
    println!("Basis functions: {}", n);

    // Warmup
    println!("\nWarmup run...");
    let mut count = 0;
    for i in 0..n {
        for j in 0..=i {
            for k in 0..n {
                for l in 0..=k {
                    if i * (i + 1) / 2 + j >= k * (k + 1) / 2 + l {
                        let _ = electron_repulsion_hgp_opt(
                            &basis[i], &basis[j], &basis[k], &basis[l]
                        );
                        count += 1;
                    }
                }
            }
        }
    }
    println!("Computed {} ERIs\n", count);

    println!("Starting profiling run (10 iterations)...\n");

    // Start profiler
    let guard = pprof::ProfilerGuardBuilder::default()
        .frequency(1000) // Sample at 1000 Hz
        .blocklist(&["libc", "libgcc", "pthread", "vdso"])
        .build()
        .unwrap();

    // Run multiple iterations for better profiling data
    let iterations = 10;
    let start = Instant::now();

    for iter in 0..iterations {
        let iter_start = Instant::now();
        let mut iter_count = 0;

        for i in 0..n {
            for j in 0..=i {
                for k in 0..n {
                    for l in 0..=k {
                        if i * (i + 1) / 2 + j >= k * (k + 1) / 2 + l {
                            let _ = electron_repulsion_hgp_opt(
                                &basis[i], &basis[j], &basis[k], &basis[l]
                            );
                            iter_count += 1;
                        }
                    }
                }
            }
        }

        let iter_elapsed = iter_start.elapsed().as_secs_f64() * 1000.0;
        println!("Iteration {}/{}: {:.2} ms ({} ERIs)",
                 iter + 1, iterations, iter_elapsed, iter_count);
    }

    let total_elapsed = start.elapsed().as_secs_f64() * 1000.0;
    let avg_elapsed = total_elapsed / iterations as f64;

    println!("\nTotal time: {:.2} ms", total_elapsed);
    println!("Average per iteration: {:.2} ms\n", avg_elapsed);

    // Stop profiler and generate flamegraph
    println!("Generating flamegraph...");

    if let Ok(report) = guard.report().build() {
        // Write flamegraph
        let file = std::fs::File::create("flamegraph_serial.svg").unwrap();
        report.flamegraph(file).unwrap();
        println!("✅ Flamegraph saved to: flamegraph_serial.svg");

        // Print simplified profiling summary
        println!("\n=== Profiling Summary ===\n");
        println!("Profiling data captured with {} total data points", report.data.len());
        println!("Review flamegraph_serial.svg for detailed bottleneck analysis");
    } else {
        println!("❌ Failed to generate profiling report");
    }

    println!("\n=== Profiling Complete ===");
    println!("\nTo view the flamegraph:");
    println!("  open flamegraph_serial.svg");
}
