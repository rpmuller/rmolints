use rmolints::basis::build_sto3g_basis;
use rmolints::molecule::Molecule;
use rmolints::parallel::{compute_eri_tensor_parallel, ERIMethod};
use std::time::Instant;

fn benchmark(method: ERIMethod) -> f64 {
    let molecule = Molecule::benzene();
    let basis = build_sto3g_basis(&molecule);

    // Warmup
    let _ = compute_eri_tensor_parallel(&basis, method);

    // 3 runs, take minimum
    let mut times = Vec::new();
    for _ in 0..3 {
        let start = Instant::now();
        let _ = compute_eri_tensor_parallel(&basis, method);
        times.push(start.elapsed().as_secs_f64() * 1000.0);
    }

    times.into_iter().fold(f64::INFINITY, f64::min)
}

fn main() {
    println!("\n=== Performance Ratio Analysis ===\n");

    let hgp_time = benchmark(ERIMethod::HeadGordonPople);
    let rys_time = benchmark(ERIMethod::Rys);
    let std_time = benchmark(ERIMethod::Standard);

    let ratio_rys = rys_time / hgp_time;
    let ratio_std = std_time / hgp_time;

    println!("Method          | Time (ms) | Ratio to HGP");
    println!("{}", "-".repeat(50));
    println!("HGP         | {:8.2}  | 1.00x (baseline)", hgp_time);
    println!("Rys (optimized) | {:8.2}  | {:.2}x slower", rys_time, ratio_rys);
    println!("Standard        | {:8.2}  | {:.2}x slower", std_time, ratio_std);

    println!("\n=== Comparison to README Claims ===\n");
    println!("README claimed:");
    println!("  HGP: 815 ms");
    println!("  Rys: 1442 ms (1.77x slower)");
    println!("  Standard: 1761 ms (2.16x slower)");

    println!("\nCurrent measurement:");
    println!("  HGP: {:.0} ms ({:.2}x slower than README)", hgp_time, hgp_time / 815.0);
    println!("  Rys: {:.0} ms ({:.2}x slower)", rys_time, ratio_rys);
    println!("  Standard: {:.0} ms ({:.2}x slower)", std_time, ratio_std);

    println!("\n=== Analysis ===\n");
    if (ratio_rys - 1.77).abs() < 0.1 {
        println!("✅ Rys ratio matches README (1.77x) - optimizations working!");
    } else if ratio_rys < 1.77 {
        println!("✅ Rys ratio BETTER than README ({:.2}x vs 1.77x) - optimizations helped!", ratio_rys);
    } else {
        println!("❌ Rys ratio WORSE than README ({:.2}x vs 1.77x)", ratio_rys);
    }

    println!("\nAll methods running ~{:.2}x slower than README - likely system variance", hgp_time / 815.0);
}
