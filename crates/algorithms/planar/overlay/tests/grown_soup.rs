//! A soup of grown triangles whose union is one polygon (#338): outputs
//! pinned bit for bit, and the memory budget of `overlay_within` and
//! `union_soup_within` measured with a counting allocator.
//!
//! The soup is axioval/engine's: the strips of a faceted bar, each corner
//! grown by a 16-gon, hulled and fanned (`grown`). With the circumscribed
//! 16-gon the hulls' flat sides line up, and before #338 the chain
//! reduction gave up on the whole soup and its subdivision grew far past
//! its output: 64 facets took 0.98 GB, 128 facets 4 GB. The hashes were
//! taken from `axiolid-overlay` 0.3.18, before #338.
//!
//! For the budget, each scene's unbudgeted call's heap high-water mark
//! above what was held before it is `P`. Under every budget tried, from an
//! eighth of `P` to many times it, a call either gives the unbudgeted
//! result bit for bit or is refused with `BudgetExceeded`, and in both
//! cases its heap high-water mark stays within the budget: the charges
//! come before the allocations and bound them. A budget of `ADMITS` times
//! `P` is always admitted, so the charges are not wildly above what is
//! held.
//!
//! The tests take one lock, so that none allocates while another is
//! measured.

use std::alloc::{GlobalAlloc, Layout, System};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;

use axiolid_contracts::ExecutionOptions;
use axiolid_core::{Frame2, Point2, Tolerance};
use axiolid_overlay::{
    overlay, overlay_with, overlay_within, union_soup, union_soup_with, union_soup_within,
    FillRule, OverlayError, OverlayInput, OverlayOperation, OverlayOptions, OverlayRefusal,
    Polygon, Ring,
};

/// The bytes now allocated and the most allocated since the last reset.
struct Counting;

static HELD: AtomicUsize = AtomicUsize::new(0);
static PEAK: AtomicUsize = AtomicUsize::new(0);

unsafe impl GlobalAlloc for Counting {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        let p = unsafe { System.alloc(layout) };
        if !p.is_null() {
            let now = HELD.fetch_add(layout.size(), Ordering::SeqCst) + layout.size();
            PEAK.fetch_max(now, Ordering::SeqCst);
        }
        p
    }

    unsafe fn dealloc(&self, p: *mut u8, layout: Layout) {
        unsafe { System.dealloc(p, layout) };
        HELD.fetch_sub(layout.size(), Ordering::SeqCst);
    }

    unsafe fn realloc(&self, p: *mut u8, layout: Layout, size: usize) -> *mut u8 {
        // Counted as the allocation of the new block before the old one
        // is freed, as a moving reallocation holds both.
        let q = unsafe { System.realloc(p, layout, size) };
        if !q.is_null() {
            let now = HELD.fetch_add(size, Ordering::SeqCst) + size;
            PEAK.fetch_max(now, Ordering::SeqCst);
            HELD.fetch_sub(layout.size(), Ordering::SeqCst);
        }
        q
    }
}

#[global_allocator]
static ALLOCATOR: Counting = Counting;

/// Held by each test, so that measurements do not overlap.
static SERIAL: Mutex<()> = Mutex::new(());

/// `run`'s value and its heap high-water mark above what was held before.
fn measured<T>(run: impl FnOnce() -> T) -> (T, usize) {
    let before = HELD.load(Ordering::SeqCst);
    PEAK.store(before, Ordering::SeqCst);
    let value = run();
    (value, PEAK.load(Ordering::SeqCst) - before)
}

/// A budget of this many times the unbudgeted high-water mark is always
/// admitted.
const ADMITS: usize = 4;

fn tolerance() -> Tolerance {
    Tolerance::new(1.0e-9, 1.0e-9).expect("a tolerance")
}

fn ccw(mut p: Vec<(f64, f64)>) -> Ring {
    let a = (p[1].0 - p[0].0) * (p[2].1 - p[0].1) - (p[1].1 - p[0].1) * (p[2].0 - p[0].0);
    if a < 0.0 {
        p.reverse();
    }
    Ring {
        points: p.into_iter().map(|(x, y)| Point2::new(x, y)).collect(),
    }
}

/// #338's strips: two triangles per facet of a bar of `n` facets.
fn strips(n: usize, y0: f64, radius: f64) -> Vec<[(f64, f64); 3]> {
    let y = |i: usize| {
        let t = std::f64::consts::TAU * (i % n) as f64 / n as f64 + 0.1;
        y0 + radius * t.cos()
    };
    let mut out = Vec::new();
    for i in 0..n {
        let (a, b) = (y(i), y(i + 1));
        if a == b {
            continue;
        }
        out.push([(0.0, a), (0.0, b), (6.0, b)]);
        out.push([(0.0, a), (6.0, b), (6.0, a)]);
    }
    out
}

/// #338's grown clip soup: each strip triangle's corners grown by a
/// 16-gon, hulled and fanned from the hull's first vertex.
fn grown(n: usize, circumscribed: bool) -> Vec<Ring> {
    let step = std::f64::consts::PI / 8.0;
    let (offset, r) = if circumscribed {
        (step / 2.0, 0.051 / (step / 2.0).cos())
    } else {
        (0.0, 0.049)
    };
    let disc: Vec<(f64, f64)> = (0..16)
        .map(|k| {
            let a = offset + step * k as f64;
            (r * a.cos(), r * a.sin())
        })
        .collect();
    let cross = |o: (f64, f64), a: (f64, f64), b: (f64, f64)| {
        (a.0 - o.0) * (b.1 - o.1) - (a.1 - o.1) * (b.0 - o.0)
    };
    let mut out = Vec::new();
    for t in strips(n, 0.01, 0.12) {
        let mut pts: Vec<(f64, f64)> = t
            .iter()
            .flat_map(|p| disc.iter().map(move |d| (p.0 + d.0, p.1 + d.1)))
            .collect();
        pts.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.total_cmp(&b.1)));
        pts.dedup();
        let mut h: Vec<(f64, f64)> = Vec::new();
        for pass in 0..2 {
            let floor = h.len() + 1;
            let it: Vec<_> = if pass == 1 {
                pts.iter().rev().skip(1).copied().collect()
            } else {
                pts.clone()
            };
            for p in it {
                while h.len() >= floor.max(2) && cross(h[h.len() - 2], h[h.len() - 1], p) <= 0.0 {
                    h.pop();
                }
                h.push(p);
            }
        }
        h.pop();
        for w in h[1..].windows(2) {
            out.push(ccw(vec![h[0], w[0], w[1]]));
        }
    }
    out
}

/// `k` squares of side 0.3 at seeded places in the unit square: a soup
/// that crosses itself about `k^2` times.
fn squares(k: usize) -> Vec<Ring> {
    let mut seed = 0x9e37_79b9_7f4a_7c15_u64;
    let mut next = || {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        (seed >> 11) as f64 / (1u64 << 53) as f64
    };
    (0..k)
        .map(|_| {
            let (x, y) = (next(), next());
            ccw(vec![(x, y), (x + 0.3, y), (x + 0.3, y + 0.3), (x, y + 0.3)])
        })
        .collect()
}

fn input(rings: Vec<Ring>) -> OverlayInput {
    OverlayInput {
        frame: Frame2 {
            origin: Point2::ZERO,
            x: Point2::X,
            y: Point2::Y,
        },
        polygons: rings
            .into_iter()
            .map(|outer| Polygon {
                outer,
                holes: Vec::new(),
            })
            .collect(),
    }
}

/// One scene's call under an optional budget.
type Scene<'a> = Box<dyn Fn(Option<usize>) -> Result<Vec<Polygon>, OverlayRefusal> + 'a>;

fn execution(budget: Option<usize>) -> ExecutionOptions {
    let options = ExecutionOptions::new(tolerance());
    match budget {
        Some(bytes) => options.with_memory_budget(bytes),
        None => options,
    }
}

fn check(name: &str, scene: &Scene<'_>, unbudgeted: &[Polygon]) {
    let (free, peak) = measured(|| scene(None));
    assert_eq!(
        free.as_deref(),
        Ok(unbudgeted),
        "{name}: no budget, no change"
    );
    let mut budgets: Vec<usize> = (0..12).map(|k| (peak >> 3) << k).collect();
    budgets.extend([1, 1 << 10, peak - 1, peak, ADMITS * peak, usize::MAX]);
    for budget in budgets {
        let (outcome, held) = measured(|| scene(Some(budget)));
        assert!(
            held <= budget,
            "{name}: held {held} bytes under a budget of {budget}"
        );
        match outcome {
            Ok(polygons) => assert_eq!(polygons, unbudgeted, "{name}: admitted under {budget}"),
            Err(OverlayRefusal::BudgetExceeded {
                needed,
                budget: refused,
            }) => {
                assert_eq!(refused, budget, "{name}");
                assert!(needed > budget, "{name}: {needed} charged under {budget}");
                assert!(
                    budget < ADMITS * peak,
                    "{name}: {budget} refused, {peak} held without a budget"
                );
            }
            Err(other) => panic!("{name}: {other:?} under {budget}"),
        }
    }
}

#[test]
fn a_budget_bounds_the_heap_and_refuses_by_name() {
    let _serial = SERIAL
        .lock()
        .unwrap_or_else(|poisoned| poisoned.into_inner());
    let t = tolerance();
    let scenes: Vec<(String, Vec<Ring>)> = vec![
        ("inscribed 32".into(), grown(32, false)),
        ("circumscribed 32".into(), grown(32, true)),
        ("circumscribed 64".into(), grown(64, true)),
        ("squares 120".into(), squares(120)),
        ("squares 300".into(), squares(300)),
    ];
    for (name, rings) in &scenes {
        let unbudgeted = union_soup(rings, t).expect("the union");
        let union: Scene<'_> = Box::new(|budget| {
            union_soup_within(rings, OverlayOptions::new(), &execution(budget))
                .map(|reported| reported.polygons)
        });
        check(&format!("union of {name}"), &union, &unbudgeted);

        let clip = input(rings.clone());
        let subject = subject(64);
        let operation = OverlayOperation::Intersection;
        let unbudgeted = overlay_with(
            &subject,
            &clip,
            operation,
            FillRule::NonZero,
            t,
            OverlayOptions::new(),
        )
        .expect("the intersection")
        .result
        .polygons;
        let intersection: Scene<'_> = Box::new(|budget| {
            overlay_within(
                &subject,
                &clip,
                operation,
                FillRule::NonZero,
                OverlayOptions::new(),
                &execution(budget),
            )
            .map(|reported| reported.result.polygons)
        });
        check(
            &format!("intersection with {name}"),
            &intersection,
            &unbudgeted,
        );
    }
}

/// An invalid operand is refused by its own name under a budget that
/// holds its operands; a budget too small even for them refuses first.
#[test]
fn an_invalid_operand_is_refused_as_without_a_budget() {
    let _serial = SERIAL
        .lock()
        .unwrap_or_else(|poisoned| poisoned.into_inner());
    let bow = Ring {
        points: vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 1.0),
            Point2::new(1.0, 0.0),
            Point2::new(0.0, 1.0),
        ],
    };
    let soup = [bow];
    assert_eq!(
        union_soup_with(&soup, tolerance(), OverlayOptions::new()),
        Err(OverlayError::SelfIntersection)
    );
    for budget in [None, Some(1 << 20), Some(usize::MAX)] {
        assert_eq!(
            union_soup_within(&soup, OverlayOptions::new(), &execution(budget)),
            Err(OverlayRefusal::Invalid(OverlayError::SelfIntersection)),
            "{budget:?}"
        );
    }
    assert!(matches!(
        union_soup_within(&soup, OverlayOptions::new(), &execution(Some(1))),
        Err(OverlayRefusal::BudgetExceeded { budget: 1, .. })
    ));
}

/// FNV-1a over 64-bit words of every output coordinate's bits, as
/// `benches/grown_soup.rs` prints it.
fn fingerprint(polygons: &[Polygon]) -> u64 {
    let mut h = 0xcbf2_9ce4_8422_2325_u64;
    let mut eat = |v: u64| {
        h ^= v;
        h = h.wrapping_mul(0x0100_0000_01b3);
    };
    for p in polygons {
        for r in std::iter::once(&p.outer).chain(&p.holes) {
            eat(r.points.len() as u64);
            for q in &r.points {
                eat(q.x.to_bits());
                eat(q.y.to_bits());
            }
        }
    }
    h
}

fn subject(n: usize) -> OverlayInput {
    input(
        strips(n, 0.0, 0.1)
            .into_iter()
            .map(|t| ccw(t.to_vec()))
            .collect(),
    )
}

/// The union, the intersection with the subject strips and the union of
/// the soup with itself, against 0.3.18's outputs.
#[test]
fn grown_soups_are_unchanged() {
    let _serial = SERIAL
        .lock()
        .unwrap_or_else(|poisoned| poisoned.into_inner());
    let t = tolerance();
    let cases: [(&str, usize, bool, u64); 7] = [
        ("soup", 64, false, 0xb04a_7664_9cfb_d6a7),
        ("soup", 32, true, 0xc166_503c_25d4_c781),
        ("soup", 64, true, 0x563b_3df6_a5b1_6e53),
        ("inter", 64, false, 0xb99d_fd16_f9f5_7fe9),
        ("inter", 32, true, 0x34fe_2012_a484_3b0b),
        ("inter", 64, true, 0x9055_6315_d5e5_1973),
        ("self", 32, true, 0xc166_503c_25d4_c781),
    ];
    for (op, n, circumscribed, pinned) in cases {
        let rings = grown(n, circumscribed);
        let clip = input(rings.clone());
        let out = match op {
            "soup" => union_soup(&rings, t).expect("the union"),
            "inter" => {
                overlay(
                    &subject(n),
                    &clip,
                    OverlayOperation::Intersection,
                    FillRule::NonZero,
                    t,
                )
                .expect("the intersection")
                .polygons
            }
            _ => {
                overlay(&clip, &clip, OverlayOperation::Union, FillRule::NonZero, t)
                    .expect("the self-union")
                    .polygons
            }
        };
        assert_eq!(
            fingerprint(&out),
            pinned,
            "{op} of {n} facets, circumscribed {circumscribed}: {out:?}"
        );
    }
}

/// 128 facets, circumscribed: 0.3.18 held about 4 GB for the union and
/// 1.9 GB for the intersection. Both now fit a budget of 256 MB, with the
/// same outputs.
#[test]
fn a_large_grown_soup_fits_a_small_budget() {
    let _serial = SERIAL
        .lock()
        .unwrap_or_else(|poisoned| poisoned.into_inner());
    let budget = 256 << 20;
    let rings = grown(128, true);
    let union = union_soup_within(&rings, OverlayOptions::new(), &execution(Some(budget)))
        .expect("within the budget");
    assert_eq!(fingerprint(&union.polygons), 0x563b_3df6_a5b1_6e53);
    let intersection = overlay_within(
        &subject(128),
        &input(rings),
        OverlayOperation::Intersection,
        FillRule::NonZero,
        OverlayOptions::new(),
        &execution(Some(budget)),
    )
    .expect("within the budget");
    assert_eq!(
        fingerprint(&intersection.result.polygons),
        0x9055_6315_d5e5_1973
    );
}
