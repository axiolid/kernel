//! Measure `boolmesh` peak scratch so its declared bound is evidence-based.
//!
//! ADR 0017 section 4 says `ScratchRequirement::Unbounded` is a declared
//! deficiency, not a resting state. Replacing it requires a *measured* bound,
//! not a guessed one, so this wraps the global allocator with a counter and
//! measures peak bytes per input triangle across representative workloads.
//!
//! Shared by the `scratch_probe` binary, which prints the table, and the
//! `scratch_bound` test, which fails when the provider's declared bound is
//! below a measured peak (#110). Run the table with:
//! ```text
//! cargo run --release -p axiolid-mesh-boolean-boolmesh --bin scratch_probe --all-features
//! ```

use std::alloc::{GlobalAlloc, Layout, System};
use std::sync::atomic::{AtomicUsize, Ordering};

use axiolid_contracts::ExecutionOptions;
use axiolid_core::{BooleanOperator, Point3, Tolerance};
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_contract::MeshBoolean;

use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;

static LIVE: AtomicUsize = AtomicUsize::new(0);
static PEAK: AtomicUsize = AtomicUsize::new(0);

/// Counting allocator: tracks live bytes and the high-water mark.
struct Counting;

unsafe impl GlobalAlloc for Counting {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        let pointer = unsafe { System.alloc(layout) };
        if !pointer.is_null() {
            let live = LIVE.fetch_add(layout.size(), Ordering::Relaxed) + layout.size();
            PEAK.fetch_max(live, Ordering::Relaxed);
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

/// One measured boolean.
pub struct Sample {
    pub elements: usize,
    pub peak: usize,
    pub operation: BooleanOperator,
}

/// Peak scratch of every operation on two overlapping boxes, subdivided to
/// 24, 96, 384 and 1,536 input triangles.
///
/// One discarded boolean runs first: the first call in a process pays for
/// allocator arena growth, lazy statics and first-touch pages, and used to
/// charge them all to whichever operation happened to be measured first
/// (#110). Peaks are measured above the bytes already live when the call
/// starts, so nothing allocated before it is counted.
pub fn measure() -> Vec<Sample> {
    let provider = BoolmeshBoolean::new();
    let options = ExecutionOptions::new(Tolerance::METRE);
    let warmup = box_at([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    let _ = provider.boolean(
        &warmup,
        &box_at([0.5, 0.5, 0.5], [1.5, 1.5, 1.5]),
        BooleanOperator::Union,
        &options,
    );

    let mut samples = Vec::new();
    for levels in 0..=3 {
        let subject = subdivide(&box_at([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]), levels);
        let tool = subdivide(&box_at([0.5, 0.5, 0.5], [1.5, 1.5, 1.5]), levels);
        let elements = subject.triangle_count() + tool.triangle_count();

        for operation in BooleanOperator::ALL {
            let baseline = LIVE.load(Ordering::Relaxed);
            PEAK.store(baseline, Ordering::Relaxed);
            let outcome = provider.boolean(&subject, &tool, operation, &options);
            let peak = PEAK.load(Ordering::Relaxed) - baseline;
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
