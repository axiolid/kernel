//! Print the host's detected shape. Used to confirm detection returns real
//! values rather than only self-consistent ones.

fn main() {
    let features = axiolid_backend_cpu::CpuFeatures::detect();
    let topology = axiolid_backend_cpu::CpuTopology::detect();
    println!("features: {features:?}");
    println!("best isa: {:?}", features.best());
    println!("logical cpus: {:?}", topology.logical_cpus);
    println!("heterogeneous: {:?}", topology.heterogeneous_cores);
    for cache in &topology.caches {
        println!(
            "L{} {} KiB, line {} B, shared by {}",
            cache.level,
            cache.bytes / 1024,
            cache.line_bytes,
            cache.shared_by
        );
    }
}
