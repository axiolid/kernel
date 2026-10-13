//! Measure the batch paths against the budget they admit (#259).
//!
//! A batch has no bound in its input triangles: `union_many` intermediates
//! outgrow their operands, and with `parallel-batch` stolen work puts more
//! booleans in flight than there are workers. So the provider charges each
//! step when its operands exist and refuses the first step that does not fit
//! (`src/budget.rs`). What must hold is that a budget the batch admits is
//! never exceeded, and that a refused batch stopped before exceeding it.
//!
//! For each workload this finds the smallest budget the batch admits on the
//! schedule's worker count (by bisection, outside any measured window, on a
//! warm pool of the same width), then measures one run at that budget and
//! one a byte below it under the forced worst-case schedule of `scratch.rs`.
//! That compares the peak with the provider's largest step charge exactly,
//! without exposing the charge.
//!
//! The workloads grow: overlapping solids, so every union has more triangles
//! than either operand, and a subject cut by overlapping tools, so it gains
//! triangles with every group. Each is measured at several sizes.

use axiolid_contracts::{ExecutionOptions, GeomError, GeomResult};
use axiolid_core::Tolerance;
use axiolid_mesh::{AttributeChannel, Blend, TriMesh};
use axiolid_mesh_boolean_contract::{BooleanOutcome, MeshBoolean};

use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;

use crate::scratch::{box_at, scheduled, subdivide, Schedule};

/// Which batch path a workload runs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Batch {
    /// `MeshBoolean::union_many`.
    UnionMany,
    /// `MeshBoolean::subtract_many`.
    SubtractMany,
}

/// One batch input.
pub struct Workload {
    /// What it is, for messages.
    pub name: String,
    /// The path it runs.
    pub batch: Batch,
    /// `union_many`'s solids, or `subtract_many`'s subject followed by its
    /// tools.
    pub meshes: Vec<TriMesh>,
}

impl Workload {
    /// Input triangles, all meshes together.
    pub fn elements(&self) -> usize {
        self.meshes.iter().map(TriMesh::triangle_count).sum()
    }

    fn run(&self, budget: usize) -> GeomResult<BooleanOutcome> {
        let provider = BoolmeshBoolean::new();
        let options = ExecutionOptions::new(Tolerance::METRE).with_memory_budget(budget);
        match self.batch {
            Batch::UnionMany => provider.union_many(&self.meshes, &options),
            Batch::SubtractMany => {
                provider.subtract_many(&self.meshes[0], &self.meshes[1..], &options)
            }
        }
    }
}

/// One measured workload.
pub struct BatchSample {
    /// The workload's name.
    pub name: String,
    /// Input triangles.
    pub elements: usize,
    /// Smallest budget the batch admits on this schedule's width.
    pub admitted: usize,
    /// Peak scratch of the run at `admitted`.
    pub peak: usize,
    /// Peak scratch of the run refused at `admitted - 1`.
    pub refused_peak: usize,
    /// Peak scratch of the run refused at [`SMALL_BUDGET`].
    pub small_peak: usize,
}

/// A budget below the copies and the bookkeeping of most workloads, which a
/// batch must refuse before allocating them.
pub const SMALL_BUDGET: usize = 1024;

/// `n` unit boxes in a row, each overlapping its neighbour by 40% and
/// stepped up and down by a tenth, subdivided `levels` times; with
/// `channel`, every solid carries a per-vertex id channel.
fn overlapping_row(n: usize, levels: usize, channel: bool) -> Vec<TriMesh> {
    (0..n)
        .map(|i| {
            let x = 0.6 * i as f64;
            let lift = 0.1 * (i % 2) as f64;
            let mut mesh = subdivide(
                &box_at([x, lift, lift], [x + 1.0, 1.0 + lift, 1.0 + lift]),
                levels,
            );
            if channel {
                let ids = vec![i as f64; mesh.positions.len()];
                mesh.attributes
                    .push(AttributeChannel::new("id", ids, 1, Blend::Nearest));
            }
            mesh
        })
        .collect()
}

/// A slab `n` long, subdivided `levels` times, followed by `n` tools. With
/// `overlapping`, each tool overlaps the next three, so the tools fall into
/// several groups and the subject is cut, and grows, once per group;
/// otherwise the tools are disjoint and fuse into one group.
fn slab_with_tools(n: usize, levels: usize, overlapping: bool) -> Vec<TriMesh> {
    let length = n as f64;
    let mut meshes = vec![subdivide(
        &box_at([0.0, 0.0, 0.0], [length, 1.0, 1.0]),
        levels,
    )];
    for i in 0..n {
        let (x0, x1) = if overlapping {
            (0.3 * i as f64 + 0.1, 0.3 * i as f64 + 1.1)
        } else {
            (i as f64 + 0.2, i as f64 + 0.7)
        };
        meshes.push(box_at([x0, 0.2, 0.5], [x1, 0.8, 1.5]));
    }
    meshes
}

/// The measured batch workloads, at several sizes.
pub fn workloads() -> Vec<Workload> {
    let mut out = Vec::new();
    for (n, levels) in [(2, 0), (5, 0), (8, 0), (16, 0), (8, 1), (4, 2)] {
        out.push(Workload {
            name: format!("union_many of {n} overlapping boxes, {levels} subdivisions"),
            batch: Batch::UnionMany,
            meshes: overlapping_row(n, levels, false),
        });
    }
    out.push(Workload {
        name: "union_many of 32 disjoint boxes".to_owned(),
        batch: Batch::UnionMany,
        meshes: (0..32)
            .map(|i| {
                let x = 2.0 * i as f64;
                box_at([x, 0.0, 0.0], [x + 1.0, 1.0, 1.0])
            })
            .collect(),
    });
    out.push(Workload {
        name: "union_many of 2 overlapping boxes with an id channel".to_owned(),
        batch: Batch::UnionMany,
        meshes: overlapping_row(2, 0, true),
    });
    out.push(Workload {
        name: "union_many of 5 overlapping boxes with an id channel".to_owned(),
        batch: Batch::UnionMany,
        meshes: overlapping_row(5, 0, true),
    });
    for (n, levels, overlapping) in [
        (4, 0, true),
        (8, 0, true),
        (8, 1, true),
        (8, 0, false),
        (16, 1, false),
    ] {
        let tools = if overlapping {
            "overlapping"
        } else {
            "disjoint"
        };
        out.push(Workload {
            name: format!("subtract_many of {n} {tools} tools from a slab, {levels} subdivisions"),
            batch: Batch::SubtractMany,
            meshes: slab_with_tools(n, levels, overlapping),
        });
    }
    out
}

/// Whether `workload` admits `budget`; any failure other than a budget
/// refusal is a broken fixture.
fn admits(workload: &Workload, budget: usize) -> bool {
    match workload.run(budget) {
        Ok(_) => true,
        Err(GeomError::BudgetExceeded { resource: "memory" }) => false,
        Err(other) => panic!("{} failed: {other:?}", workload.name),
    }
}

/// Smallest budget `workload` admits when run on `workers` threads.
///
/// Also measures the run refused at [`SMALL_BUDGET`], on the same warm
/// workers: their start-up belongs to the per-worker term, which a batch
/// refused before its first boolean never reaches, and a cold pool would
/// charge it to this window.
fn smallest_admitted(workload: &Workload, schedule: Schedule) -> (usize, usize) {
    let search = || {
        let mut high = 1usize << 20;
        while !admits(workload, high) {
            high = high.checked_mul(2).expect("a finite admitted budget");
        }
        let mut low = 0;
        while low + 1 < high {
            let middle = low + (high - low) / 2;
            if admits(workload, middle) {
                high = middle;
            } else {
                low = middle;
            }
        }
        // After the search, so every worker has already made its first-use
        // allocations (its deque's epoch registration, say).
        let (outcome, small_peak) =
            crate::scratch::peak_of(|| workload.run(SMALL_BUDGET).map(drop));
        assert!(
            matches!(outcome, Err(GeomError::BudgetExceeded { .. })),
            "{} admitted {SMALL_BUDGET} bytes",
            workload.name
        );
        (high, small_peak)
    };
    match schedule {
        Schedule::Inline => search(),
        #[cfg(any(feature = "parallel", feature = "parallel-batch"))]
        Schedule::ColdPool { workers } => {
            // Joined before returning, like the measurement pools: a worker
            // exiting later would free its bookkeeping inside the next
            // measured window and hide that much of a real peak.
            let mut handles = Vec::with_capacity(workers);
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(workers)
                .spawn_handler(|thread| {
                    handles.push(std::thread::Builder::new().spawn(|| thread.run())?);
                    Ok(())
                })
                .build()
                .expect("search pool");
            let found = pool.install(search);
            drop(pool);
            for handle in handles {
                handle.join().expect("search worker");
            }
            found
        }
        #[cfg(not(any(feature = "parallel", feature = "parallel-batch")))]
        Schedule::ColdPool { .. } => unreachable!("no rayon in this build"),
    }
}

/// Measure every workload under `schedule` at the smallest budget it admits
/// and one byte below it.
pub fn measure(schedule: Schedule) -> Vec<BatchSample> {
    workloads()
        .into_iter()
        .map(|workload| {
            let (admitted, small_peak) = smallest_admitted(&workload, schedule);
            let (outcome, peak) = scheduled(schedule, || workload.run(admitted).map(drop));
            if let Err(error) = outcome {
                panic!(
                    "{} refused the budget it admitted ({admitted} bytes) under {schedule:?}: \
                     {error:?}",
                    workload.name
                );
            }
            let (outcome, refused_peak) =
                scheduled(schedule, || workload.run(admitted - 1).map(drop));
            assert!(
                matches!(outcome, Err(GeomError::BudgetExceeded { .. })),
                "{} admitted {} bytes under {schedule:?} but not in the search",
                workload.name,
                admitted - 1
            );
            BatchSample {
                name: workload.name.clone(),
                elements: workload.elements(),
                admitted,
                peak,
                refused_peak,
                small_peak,
            }
        })
        .collect()
}
