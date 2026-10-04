//! Measure `boolmesh` peak scratch so its declared bound is evidence-based.
//!
//! ADR 0017 section 4 says `ScratchRequirement::Unbounded` is a declared
//! deficiency, not a resting state. Replacing it requires a *measured* bound,
//! not a guessed one, so this wraps the global allocator with a counter and
//! measures peak bytes across representative workloads.
//!
//! Shared by the `scratch_probe` binary, which prints the table, and the
//! `scratch_bound` test, which fails when the provider's declared bound is
//! below a measured peak (#110). Run the table with:
//! ```text
//! cargo run --release -p axiolid-mesh-boolean-boolmesh --bin scratch_probe --all-features
//! ```
//!
//! # Forced worst-case schedule (#226)
//!
//! With rayon enabled a boolean's peak depends on the schedule, and an idle
//! machine shows the cheapest one. Two costs only appear under load:
//!
//! - Worker start-up. Each rayon worker allocates its bookkeeping (deque,
//!   registry slot, thread-locals) on its own thread, after the pool is
//!   built. On a loaded machine those allocations land late, inside
//!   whichever boolean is running, and stay live. This was the flake: a
//!   24-triangle `SymmetricDifference` charged with most of the 20-thread
//!   global pool's start-up.
//! - Concurrent halves. `boolean03` forks with `rayon::join`; when another
//!   worker steals a half, both halves' scratch is live at once instead of
//!   one after the other.
//!
//! [`Schedule::ColdPool`] forces both instead of waiting for load to produce
//! them: every boolean runs in a pool built inside the measured window, all
//! of whose workers have finished starting before the boolean begins, and
//! every allocation stalls the allocating thread briefly so idle workers
//! steal the other half of each join. The pool is retired and its threads
//! joined before the next sample, so no thread exit frees memory inside a
//! later window.

use std::alloc::{GlobalAlloc, Layout, System};
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

use axiolid_contracts::ExecutionOptions;
use axiolid_core::{BooleanOperator, Point3, Tolerance};
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_contract::MeshBoolean;

use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;

static LIVE: AtomicUsize = AtomicUsize::new(0);
static PEAK: AtomicUsize = AtomicUsize::new(0);
/// While set, every allocation busy-waits [`STALL_NANOS`] after it succeeds.
static STALL: AtomicBool = AtomicBool::new(false);

/// Long enough for a parked rayon worker to wake and steal (tens of
/// microseconds), short enough that 1,536 triangles still measure in about a
/// second. A busy wait rather than a sleep: the stall costs the same CPU time
/// on a loaded machine instead of whatever the scheduler grants.
const STALL_NANOS: u128 = 20_000;

/// Counting allocator: tracks live bytes and the high-water mark.
struct Counting;

unsafe impl GlobalAlloc for Counting {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        let pointer = unsafe { System.alloc(layout) };
        if !pointer.is_null() {
            let live = LIVE.fetch_add(layout.size(), Ordering::Relaxed) + layout.size();
            PEAK.fetch_max(live, Ordering::Relaxed);
            if STALL.load(Ordering::Relaxed) {
                let start = std::time::Instant::now();
                while start.elapsed().as_nanos() < STALL_NANOS {
                    std::hint::spin_loop();
                }
            }
        }
        pointer
    }

    unsafe fn dealloc(&self, pointer: *mut u8, layout: Layout) {
        LIVE.fetch_sub(layout.size(), Ordering::Relaxed);
        unsafe { System.dealloc(pointer, layout) }
    }
}

#[global_allocator]
static ALLOCATOR: Counting = Counting;

fn box_at(min: [f64; 3], max: [f64; 3]) -> TriMesh {
    let [x0, y0, z0] = min;
    let [x1, y1, z1] = max;
    let positions: Vec<Point3> = vec![
        [x0, y0, z0].into(),
        [x1, y0, z0].into(),
        [x1, y1, z0].into(),
        [x0, y1, z0].into(),
        [x0, y0, z1].into(),
        [x1, y0, z1].into(),
        [x1, y1, z1].into(),
        [x0, y1, z1].into(),
    ];
    let indices = vec![
        0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 1, 2, 6, 1, 6, 5, 2, 3, 7, 2, 7, 6,
        3, 0, 4, 3, 4, 7,
    ];
    TriMesh::new(positions, indices)
}

/// Subdivide each triangle four ways, `levels` times, to grow triangle count.
fn subdivide(mesh: &TriMesh, levels: usize) -> TriMesh {
    let mut current = mesh.clone();
    for _ in 0..levels {
        let mut positions = current.positions.clone();
        let mut indices = Vec::new();
        for triangle in current.indices.chunks_exact(3) {
            let [ia, ib, ic] = [triangle[0], triangle[1], triangle[2]];
            let (a, b, c) = (
                current.positions[ia as usize],
                current.positions[ib as usize],
                current.positions[ic as usize],
            );
            let base = positions.len() as u32;
            positions.push((a + b) * 0.5);
            positions.push((b + c) * 0.5);
            positions.push((c + a) * 0.5);
            let (ab, bc, ca) = (base, base + 1, base + 2);
            indices.extend_from_slice(&[ia, ab, ca]);
            indices.extend_from_slice(&[ab, ib, bc]);
            indices.extend_from_slice(&[ca, bc, ic]);
            indices.extend_from_slice(&[ab, bc, ca]);
        }
        current = TriMesh::new(positions, indices);
    }
    current
}

/// How a measured boolean is scheduled.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Schedule {
    /// On the calling thread with no pool: the only schedule a build without
    /// rayon has.
    Inline,
    /// The forced worst case with rayon (module docs): a pool of `workers`
    /// threads started inside the measured window, allocations stalled.
    ColdPool {
        /// Worker threads in the pool.
        workers: usize,
    },
}

impl Schedule {
    /// The schedules this build can run. With rayon: pools of 1 to 64
    /// workers, past this machine's core count, so the per-worker term is
    /// measured rather than extrapolated.
    pub fn all() -> Vec<Self> {
        if cfg!(any(feature = "parallel", feature = "parallel-batch")) {
            [1, 2, 4, 16, 64]
                .into_iter()
                .map(|workers| Self::ColdPool { workers })
                .collect()
        } else {
            vec![Self::Inline]
        }
    }

    /// Worker threads the boolean can run on, for the declared bound.
    pub fn workers(self) -> usize {
        match self {
            Self::Inline => 1,
            Self::ColdPool { workers } => workers,
        }
    }
}

/// One measured boolean.
pub struct Sample {
    pub elements: usize,
    pub peak: usize,
    pub operation: BooleanOperator,
}

/// Peak bytes allocated by `run` above what was live when it started.
fn peak_of<R>(run: impl FnOnce() -> R) -> (R, usize) {
    let baseline = LIVE.load(Ordering::Relaxed);
    PEAK.store(baseline, Ordering::Relaxed);
    let result = run();
    (result, PEAK.load(Ordering::Relaxed) - baseline)
}

/// Run `boolean` under `schedule`, returning its result and peak scratch.
///
/// For [`Schedule::ColdPool`] the pool's construction and its workers'
/// start-up are inside the window: that is what a first boolean on rayon's
/// global pool pays, and what a loaded machine can defer into any later one.
fn scheduled<R: Send>(schedule: Schedule, boolean: impl FnOnce() -> R + Send) -> (R, usize) {
    match schedule {
        Schedule::Inline => peak_of(boolean),
        #[cfg(any(feature = "parallel", feature = "parallel-batch"))]
        Schedule::ColdPool { workers } => {
            // Sized before the window so collecting handles allocates nothing
            // inside it that a real call would not.
            let mut handles = Vec::with_capacity(workers);
            let ((result, pool), peak) = peak_of(|| {
                let pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(workers)
                    .spawn_handler(|thread| {
                        handles.push(std::thread::Builder::new().spawn(|| thread.run())?);
                        Ok(())
                    })
                    .build()
                    .expect("measurement pool");
                // Every worker has run, so its start-up allocations are live
                // before the boolean begins: the superposition a loaded
                // machine can produce, not the staggered one an idle one does.
                pool.broadcast(|_| ());
                STALL.store(true, Ordering::Relaxed);
                let result = pool.install(boolean);
                STALL.store(false, Ordering::Relaxed);
                (result, pool)
            });
            drop(pool);
            for handle in handles {
                handle.join().expect("measurement worker");
            }
            (result, peak)
        }
        #[cfg(not(any(feature = "parallel", feature = "parallel-batch")))]
        Schedule::ColdPool { .. } => unreachable!("no rayon in this build"),
    }
}

/// Peak scratch of every operation on two overlapping boxes, subdivided to
/// 24, 96, 384 and 1,536 input triangles, under `schedule`.
///
/// One discarded boolean runs first, under the same schedule: the first call
/// in a process pays for allocator arena growth, lazy statics and first-touch
/// pages, and used to charge them all to whichever operation happened to be
/// measured first (#110). Peaks are measured above the bytes already live
/// when the call starts, so nothing allocated before it is counted.
pub fn measure(schedule: Schedule) -> Vec<Sample> {
    let provider = BoolmeshBoolean::new();
    let options = ExecutionOptions::new(Tolerance::METRE);
    let warmup = box_at([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    let warmup_tool = box_at([0.5, 0.5, 0.5], [1.5, 1.5, 1.5]);
    let _ = scheduled(schedule, || {
        provider
            .boolean(&warmup, &warmup_tool, BooleanOperator::Union, &options)
            .is_ok()
    });

    let mut samples = Vec::new();
    for levels in 0..=3 {
        let subject = subdivide(&box_at([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]), levels);
        let tool = subdivide(&box_at([0.5, 0.5, 0.5], [1.5, 1.5, 1.5]), levels);
        let elements = subject.triangle_count() + tool.triangle_count();

        for operation in BooleanOperator::ALL {
            let (outcome, peak) = scheduled(schedule, || {
                provider.boolean(&subject, &tool, operation, &options)
            });
            assert!(outcome.is_ok(), "{operation:?} failed at level {levels}");
            drop(outcome);
            samples.push(Sample {
                elements,
                peak,
                operation,
            });
        }
    }
    samples
}
