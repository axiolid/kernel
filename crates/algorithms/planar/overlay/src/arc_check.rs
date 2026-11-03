//! Exact validity of an arc polygon (#352).
//!
//! [`validate_arc_ring`] checks one ring's degeneracies within a tolerance
//! and leaves open whether the ring crosses itself or whether the holes sit
//! inside the outer boundary. [`check_arc_polygon`] decides those exactly.
//!
//! # What is decided, and how
//!
//! Every contact and nesting decision is an exact sign over the given `f64`
//! input (ADR 0070): the edges' meeting points come from `crossings`, point
//! identity from `same_point`, and containment from the exact
//! [`ArcArrangement`]. No tolerance enters any of them, so a scene gives
//! the same answer in millimetres and in metres, and a tangency that holds
//! exactly is a contact while one that misses by `2^-40` is not. The
//! tolerance only feeds the per-ring [`validate_arc_ring`].
//!
//! - **Self-contact.** Two edges of one ring may share only the vertex
//!   their ring order puts between them: consecutive edges share their
//!   joint, and the closing edge shares vertex 0 with the first. Any other
//!   shared point is a [`ArcPolygonDefect::SelfContact`]: a crossing, a
//!   pinch, a tangency or a fold-back. A fold-back has consecutive edges
//!   lying over each other, and the overlap's far end is such a point.
//! - **Nesting.** Run only when every ring passes [`validate_arc_ring`] and
//!   has no self-contact, because the arrangement requires simple rings.
//!   Pieces carried by a hole alone must lie inside the outer ring and
//!   outside every other hole. A piece carried by two rings is a
//!   [`ArcPolygonDefect::BoundaryOverlap`].
//! - **Contacts.** An isolated point where two rings meet, with the regions
//!   on each side still consistent, is an [`ArcRingContact`], not a defect:
//!   a hole touching its outer ring at a vertex, or two tangent holes, is
//!   valid under the OGC rule and under the straight overlay's own
//!   `hole_outside`, as long as the contacts do not disconnect the interior. [`ArcPolygonCheck::is_strict`] refuses them. The list
//!   holds every meeting point of two rings, also those of rings that carry
//!   a defect between them (a crossing is a meeting point too).
//!
//! - **Disconnected interior.** Contacts can cut the interior apart even
//!   though each is an isolated point: a diamond hole touching all four
//!   sides of a rectangle, or a closed chain of holes. Such a polygon is not
//!   valid under the OGC rule, whose interior must be connected. The rings and
//!   the distinct contact points form a bipartite graph; a cycle in it is a
//!   [`ArcPolygonDefect::DisconnectedInterior`]. It is judged only when no
//!   other nesting defect is present, since a crossing also makes cycles.
//!
//! Winding is not read: each ring is taken as the region it encloses.
//!
//! Output order is deterministic: defects by ring, then edge pair, then the
//! nesting defects in the order overlaps, holes outside, holes overlapping,
//! disconnected interiors;
//! contacts by ring pair, then arrangement vertex index.

use std::collections::BTreeSet;

use axiolid_core::{Point2, Tolerance};

use crate::arc::{validate_arc_ring, ArcRing};
use crate::arc_overlay::ArcPolygon;
use crate::arrangement::ArcArrangement;
use crate::exact_arc::boxes::BoxTree;
use crate::exact_arc::{crossings, edges_of, same_point};
use crate::OverlayError;

/// A ring of an [`ArcPolygon`]: the outer boundary or hole `k`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum ArcRingRef {
    Outer,
    Hole(usize),
}

impl ArcRingRef {
    fn of(index: usize) -> Self {
        match index {
            0 => Self::Outer,
            k => Self::Hole(k - 1),
        }
    }
}

/// Why an [`ArcPolygon`] does not bound a valid region.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq)]
pub enum ArcPolygonDefect {
    /// A [`validate_arc_ring`] refusal of one ring.
    Ring {
        ring: ArcRingRef,
        error: OverlayError,
    },
    /// Two edges of one ring share a point other than the vertex their ring
    /// order puts between them: a crossing, pinch, tangency or fold-back.
    /// `edges` ascending; `at` is the first such point, rounded once.
    SelfContact {
        ring: ArcRingRef,
        edges: [usize; 2],
        at: Point2,
    },
    /// Two rings share a stretch of boundary. `rings` ascending.
    BoundaryOverlap { rings: [ArcRingRef; 2] },
    /// Contacts between rings cut the interior into separate pieces: the
    /// rings and the contact points form a cycle. `rings` ascending, every
    /// ring of the cut component.
    DisconnectedInterior { rings: Vec<ArcRingRef> },
    /// The exact arrangement that decides nesting refused the rings, so
    /// nesting is undecided and the polygon is not called valid.
    NestingUndecided { error: OverlayError },
    /// Part of a hole lies outside the outer ring.
    HoleOutsideOuter { hole: usize },
    /// Two holes' interiors meet (crossing or nested). `holes` ascending.
    HolesOverlap { holes: [usize; 2] },
}

/// A point where two different rings meet.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ArcRingContact {
    /// The two rings, ascending.
    pub rings: [ArcRingRef; 2],
    /// The meeting point, rounded once.
    pub at: Point2,
}

/// The result of [`check_arc_polygon`].
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Default)]
pub struct ArcPolygonCheck {
    pub defects: Vec<ArcPolygonDefect>,
    pub contacts: Vec<ArcRingContact>,
}

impl ArcPolygonCheck {
    /// No defect. Isolated contacts between rings are allowed.
    pub fn is_valid(&self) -> bool {
        self.defects.is_empty()
    }

    /// No defect and no contact between rings.
    pub fn is_strict(&self) -> bool {
        self.defects.is_empty() && self.contacts.is_empty()
    }
}

/// Exact validity of an arc polygon: rings simple, holes nested.
///
/// See the module documentation for what is decided. `tolerance` only feeds
/// [`validate_arc_ring`]. Nesting is checked only when every ring passes it
/// and has no self-contact.
pub fn check_arc_polygon(polygon: &ArcPolygon, tolerance: Tolerance) -> ArcPolygonCheck {
    let rings: Vec<&ArcRing> = std::iter::once(&polygon.outer)
        .chain(&polygon.holes)
        .collect();
    let mut check = ArcPolygonCheck::default();
    let mut simple = true;
    for (index, ring) in rings.iter().enumerate() {
        let name = ArcRingRef::of(index);
        if let Err(error) = validate_arc_ring(ring, tolerance) {
            // A zero net area can be a figure eight, whose crossing is the
            // more useful finding; the other refusals leave no edges to test.
            let zero_area = error == OverlayError::ZeroArea;
            check
                .defects
                .push(ArcPolygonDefect::Ring { ring: name, error });
            simple = false;
            if !zero_area {
                continue;
            }
        }
        if let Some((edges, at)) = self_contact(ring) {
            check.defects.push(ArcPolygonDefect::SelfContact {
                ring: name,
                edges,
                at,
            });
            simple = false;
        }
    }
    if simple {
        nesting(&rings, tolerance, &mut check);
    }
    check
}

/// The first pair of edges of `ring` that meet other than at their shared
/// joint, with the first such point.
fn self_contact(ring: &ArcRing) -> Option<([usize; 2], Point2)> {
    let edges = edges_of(ring);
    let n = edges.len();
    let tree = BoxTree::new(edges.iter().map(|e| e.bounds).collect());
    let (mut hits, mut stack) = (Vec::new(), Vec::new());
    for i in 0..n {
        tree.query_into(|b| b.overlaps(&edges[i].bounds), &mut hits, &mut stack);
        for &j in hits.iter().filter(|&&j| j > i) {
            let next = j == i + 1;
            let closing = i == 0 && j == n - 1;
            let stray = crossings(&edges[i], &edges[j]).into_iter().find(|x| {
                !(next && same_point(x, &edges[i].p1)) && !(closing && same_point(x, &edges[0].p0))
            });
            if let Some(x) = stray {
                return Some(([i, j], x.rounded()));
            }
        }
    }
    None
}

/// Hole nesting and ring contacts, on rings that are simple one by one.
fn nesting(rings: &[&ArcRing], tolerance: Tolerance, check: &mut ArcPolygonCheck) {
    let owned: Vec<ArcRing> = rings.iter().map(|&r| r.clone()).collect();
    let arrangement = match ArcArrangement::new(&owned, tolerance) {
        Ok(arrangement) => arrangement,
        Err(error) => {
            // Rings passed `validate_arc_ring` above, so this is the
            // arrangement's own refusal; nesting stays undecided rather
            // than the polygon being called valid.
            check
                .defects
                .push(ArcPolygonDefect::NestingUndecided { error });
            return;
        }
    };
    let mut overlaps = BTreeSet::new();
    let mut outside = BTreeSet::new();
    let mut crowded = BTreeSet::new();
    let mut at_vertex: Vec<Vec<usize>> = vec![Vec::new(); arrangement.vertices().len()];
    let mut holding = Vec::new();
    for edge in arrangement.edges() {
        let mut carriers: Vec<usize> = edge.sources.iter().map(|s| s.ring).collect();
        carriers.sort_unstable();
        carriers.dedup();
        for &ring in &carriers {
            for v in [edge.from, edge.to] {
                if !at_vertex[v].contains(&ring) {
                    at_vertex[v].push(ring);
                }
            }
        }
        match carriers[..] {
            [] | [0] => {}
            [only] => {
                // The outer ring does not carry this piece, so it holds both
                // sides or neither.
                if !edge.inside_left(0) {
                    outside.insert(only - 1);
                }
                edge.holding(true, &mut holding);
                for &other in holding.iter().filter(|&&r| r != 0 && r != only) {
                    crowded.insert([only.min(other) - 1, only.max(other) - 1]);
                }
            }
            _ => {
                for (k, &a) in carriers.iter().enumerate() {
                    for &b in &carriers[k + 1..] {
                        overlaps.insert([a, b]);
                    }
                }
            }
        }
    }
    let disconnected = if overlaps.is_empty() && outside.is_empty() && crowded.is_empty() {
        cut_components(&at_vertex, rings.len())
    } else {
        Vec::new()
    };
    // (ring pair, vertex) of every meeting point.
    let mut meetings: BTreeSet<([usize; 2], usize)> = BTreeSet::new();
    for (vertex, present) in at_vertex.iter_mut().enumerate() {
        present.sort_unstable();
        for (k, &a) in present.iter().enumerate() {
            for &b in &present[k + 1..] {
                meetings.insert(([a, b], vertex));
            }
        }
    }
    check.defects.extend(
        overlaps
            .into_iter()
            .map(|[a, b]| ArcPolygonDefect::BoundaryOverlap {
                rings: [ArcRingRef::of(a), ArcRingRef::of(b)],
            }),
    );
    check.defects.extend(
        outside
            .into_iter()
            .map(|hole| ArcPolygonDefect::HoleOutsideOuter { hole }),
    );
    check.defects.extend(
        crowded
            .into_iter()
            .map(|holes| ArcPolygonDefect::HolesOverlap { holes }),
    );
    check.defects.extend(
        disconnected
            .into_iter()
            .map(|rings| ArcPolygonDefect::DisconnectedInterior { rings }),
    );
    check
        .contacts
        .extend(meetings.into_iter().map(|([a, b], vertex)| ArcRingContact {
            rings: [ArcRingRef::of(a), ArcRingRef::of(b)],
            at: arrangement.vertices()[vertex],
        }));
}

/// The ring sets of the components of the ring-and-contact-point graph that
/// hold a cycle. `at_vertex[v]` lists the rings present at arrangement vertex
/// `v`; a vertex with two or more rings is a contact point. Union-find over
/// the ring-point edges: an edge whose ends are already joined closes a cycle.
/// The graph is bipartite, so three rings meeting at one point are a star,
/// while two contacts between one pair of rings are a cycle.
fn cut_components(at_vertex: &[Vec<usize>], ring_count: usize) -> Vec<Vec<ArcRingRef>> {
    fn find(parent: &mut [usize], mut x: usize) -> usize {
        while parent[x] != x {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        x
    }
    // Nodes: rings first, then one per arrangement vertex.
    let mut parent: Vec<usize> = (0..ring_count + at_vertex.len()).collect();
    let mut closing = Vec::new();
    for (vertex, present) in at_vertex.iter().enumerate() {
        if present.len() < 2 {
            continue;
        }
        let point = ring_count + vertex;
        for &ring in present {
            let (a, b) = (find(&mut parent, ring), find(&mut parent, point));
            if a == b {
                closing.push(ring);
            } else {
                parent[a] = b;
            }
        }
    }
    let roots: BTreeSet<usize> = closing
        .into_iter()
        .map(|ring| find(&mut parent, ring))
        .collect();
    roots
        .into_iter()
        .map(|root| {
            (0..ring_count)
                .filter(|&ring| find(&mut parent, ring) == root)
                .map(ArcRingRef::of)
                .collect()
        })
        .collect()
}
