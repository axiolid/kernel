//! Booleans whose result is much larger than their operands (#308).
//!
//! The declared bound counts input triangles, so it admits a boolean on its
//! operands alone. Where the operands cross many times the result has far
//! more triangles than they do: two perpendicular combs of `k` teeth have
//! `32k + 24` triangles and cross `k^2` times, and their union has about
//! `21k^2`. The provider therefore charges the result inside the solve, once
//! the intersection pass has counted the crossings, and refuses with
//! `BudgetExceeded` before allocating it when the operands' bound plus the
//! result's charge does not fit (`src/budget.rs`).
//!
//! For each workload this measures the peak with no budget, then finds the
//! smallest budget the boolean admits on the schedule's worker count (by
//! bisection, outside any measured window, on a warm pool of the same
//! width) and measures one run at that budget and one a byte below it,
//! under the forced worst-case schedule of `scratch.rs`. It also runs the
//! boolean at its operands' declared bound, which is what a caller that
//! trusts the declaration admits: either it fits, or it is refused within
//! it. The workloads are the combs at `k = 2..32` under all four operations,
//! with and without a per-vertex channel, and random soups of crossing
//! slabs.

use axiolid_contracts::{ExecutionOptions, GeomError, GeomResult};
use axiolid_core::{BooleanOperator, Tolerance};
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_contract::{BooleanOutcome, MeshBoolean};

use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;

use crate::scratch::{scheduled, Schedule};
use crate::shapes::{box_at, combs, pins, soup, subdivide, with_channel};

/// Teeth of the measured combs.
pub const TEETH: [usize; 5] = [2, 4, 8, 16, 32];

/// Random soups measured.
pub const SOUPS: u64 = 24;

/// One crossing workload.
pub struct Crossing {
    /// What it is, for messages.
    pub name: String,
    /// The boolean it runs.
    pub operation: BooleanOperator,
    /// Subject.
    pub subject: TriMesh,
    /// Tool.
    pub tool: TriMesh,
}

impl Crossing {
    /// Input triangles, both operands.
    pub fn elements(&self) -> usize {
        self.subject.triangle_count() + self.tool.triangle_count()
    }

    /// Run the boolean, with `budget` if any.
    pub fn run(&self, budget: Option<usize>) -> GeomResult<BooleanOutcome> {
        let mut options = ExecutionOptions::new(Tolerance::METRE);
        if let Some(bytes) = budget {
            options = options.with_memory_budget(bytes);
        }
        BoolmeshBoolean::new().boolean(&self.subject, &self.tool, self.operation, &options)
    }
}

/// Two overlapping boxes subdivided 0 to 3 times, as `scratch.rs` measures
/// them, under all four operations: ordinary booleans, whose result is about
/// the size of their operands.
pub fn box_workloads() -> Vec<Crossing> {
    let mut out = Vec::new();
    for levels in 0..=3 {
        let subject = subdivide(&box_at([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]), levels);
        let tool = subdivide(&box_at([0.5, 0.5, 0.5], [1.5, 1.5, 1.5]), levels);
        for operation in BooleanOperator::ALL {
            out.push(Crossing {
                name: format!("boxes subdivided {levels} times, {operation:?}"),
                operation,
                subject: subject.clone(),
                tool: tool.clone(),
            });
        }
    }
    out
}

/// The pierced slabs: `m` of 4, 8 and 16 under the three primitive
/// operations. Not `SymmetricDifference`: its last step subtracts the
/// intersection from the union, whose pins share faces with it, and at 16
/// by 16 that leaves an edge of four faces the solve refuses as an operand.
pub fn pin_workloads() -> Vec<Crossing> {
    let mut out = Vec::new();
    for m in [4, 8, 16] {
        let (subject, tool) = pins(m);
        for operation in [
            BooleanOperator::Union,
            BooleanOperator::Intersection,
            BooleanOperator::Difference,
        ] {
            out.push(Crossing {
                name: format!("slab pierced by {m}x{m} pins, {operation:?}"),
                operation,
                subject: subject.clone(),
                tool: tool.clone(),
            });
        }
    }
    out
}

/// The comb workloads: every `k` of [`TEETH`] under all four operations,
/// and the largest two with a channel of width two.
pub fn comb_workloads() -> Vec<Crossing> {
    let mut out = Vec::new();
    for k in TEETH {
        let (subject, tool) = combs(k);
        for operation in BooleanOperator::ALL {
            out.push(Crossing {
                name: format!("combs of {k} teeth, {operation:?}"),
                operation,
                subject: subject.clone(),
                tool: tool.clone(),
            });
        }
    }
    for k in [16, 32] {
        let (subject, tool) = combs(k);
        for operation in [BooleanOperator::Union, BooleanOperator::SymmetricDifference] {
            out.push(Crossing {
                name: format!("combs of {k} teeth with a channel, {operation:?}"),
                operation,
                subject: with_channel(subject.clone(), 2),
                tool: with_channel(tool.clone(), 2),
            });
        }
    }
    out
}

/// The random soups: `count` seeds, each under one operation in turn.
pub fn soup_workloads(count: u64) -> Vec<Crossing> {
    (0..count)
        .map(|seed| {
            let (subject, tool) = soup(seed);
            let operation = BooleanOperator::ALL[seed as usize % BooleanOperator::ALL.len()];
            Crossing {
                name: format!("soup {seed}, {operation:?}"),
                operation,
                subject,
                tool,
            }
        })
        .collect()
}

/// Every measured workload: combs, soups, pierced slabs and, for contrast,
/// ordinary overlapping boxes.
pub fn workloads() -> Vec<Crossing> {
    let mut out = comb_workloads();
    out.extend(soup_workloads(SOUPS));
    out.extend(pin_workloads());
    out.extend(box_workloads());
    out
}

/// One workload's peak with no budget.
pub struct Unbudgeted {
    /// The workload's name.
    pub name: String,
    /// Input triangles.
    pub elements: usize,
    /// Result triangles.
    pub output: usize,
    /// Peak scratch.
    pub peak: usize,
}

/// Peak of every workload with no budget, under `schedule`.
pub fn unbudgeted(workloads: &[Crossing], schedule: Schedule) -> Vec<Unbudgeted> {
    workloads
        .iter()
        .map(|workload| {
            let (outcome, peak) = scheduled(schedule, || {
                workload.run(None).map(|o| o.mesh.triangle_count())
            });
            let output = outcome.unwrap_or_else(|e| panic!("{} failed: {e:?}", workload.name));
            Unbudgeted {
                name: workload.name.clone(),
                elements: workload.elements(),
                output,
                peak,
            }
        })
        .collect()
}

/// One workload measured against budgets.
pub struct Budgeted {
    /// The workload's name.
    pub name: String,
    /// Input triangles.
    pub elements: usize,
    /// Smallest budget the boolean admits on this schedule's width.
    pub admitted: usize,
    /// Peak of the run at `admitted`.
    pub peak: usize,
    /// Peak of the run refused at `admitted - 1`.
    pub refused_peak: usize,
    /// The operands' declared bound on this width.
    pub declared: usize,
    /// Peak of the run at `declared`, and whether it was refused.
    pub declared_peak: usize,
    /// Whether the run at `declared` was refused.
    pub declared_refused: bool,
}

/// Whether `workload` admits `budget`; any failure other than a budget
/// refusal is a broken fixture.
fn admits(workload: &Crossing, budget: usize) -> bool {
    match workload.run(Some(budget)) {
        Ok(_) => true,
        Err(GeomError::BudgetExceeded { resource: "memory" }) => false,
        Err(other) => panic!("{} failed: {other:?}", workload.name),
    }
}

/// Smallest budget `workload` admits on `schedule`'s width, searched on a
/// warm pool of that width outside any measured window.
fn smallest_admitted(workload: &Crossing, schedule: Schedule) -> usize {
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
        high
    };
    match schedule {
        Schedule::Inline => search(),
        #[cfg(any(feature = "parallel", feature = "parallel-batch"))]
        Schedule::ColdPool { workers } => {
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

/// Measure every workload under `schedule` at the smallest budget it
/// admits, a byte below it, and at its operands' declared bound.
pub fn budgeted(workloads: &[Crossing], schedule: Schedule) -> Vec<Budgeted> {
    let requirement = BoolmeshBoolean::new().scratch_requirement();
    workloads
        .iter()
        .map(|workload| {
            let admitted = smallest_admitted(workload, schedule);
            let (outcome, peak) = scheduled(schedule, || workload.run(Some(admitted)).map(drop));
            if let Err(error) = outcome {
                panic!(
                    "{} refused the budget it admitted ({admitted} bytes) under {schedule:?}: \
                     {error:?}",
                    workload.name
                );
            }
            let (outcome, refused_peak) =
                scheduled(schedule, || workload.run(Some(admitted - 1)).map(drop));
            assert!(
                matches!(outcome, Err(GeomError::BudgetExceeded { .. })),
                "{} admitted {} bytes under {schedule:?} but not in the search",
                workload.name,
                admitted - 1
            );
            let declared = requirement
                .upper_bound_bytes_on(workload.elements(), schedule.workers())
                .expect("boolmesh declares a finite bound");
            let (outcome, declared_peak) =
                scheduled(schedule, || workload.run(Some(declared)).map(drop));
            let declared_refused = match outcome {
                Ok(()) => false,
                Err(GeomError::BudgetExceeded { resource: "memory" }) => true,
                Err(other) => panic!("{} failed: {other:?}", workload.name),
            };
            Budgeted {
                name: workload.name.clone(),
                elements: workload.elements(),
                admitted,
                peak,
                refused_peak,
                declared,
                declared_peak,
                declared_refused,
            }
        })
        .collect()
}
