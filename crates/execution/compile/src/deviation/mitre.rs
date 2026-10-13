//! The deviation of a swept disk swept piece by piece (#232, #245, #248).
//!
//! `axiolid_construct::pipe` proves every part of the pipe within the
//! budget except where a corner is mitred beside an arc. There the exact
//! solid is each tube continued to the mitre plane and cut by it, and its
//! boundary near the corner is, on each side, the tube's wall up to the
//! plane, and in the plane the region between the two sections (where one
//! tube ends and the other does not begin). A segment's wall is proven
//! (its ring is its own mitre ring); the rest is certified here against
//! the mesh by the branch and bound of [`crate::certify`]:
//!
//! - An arc's band, `T(b, psi)`: the point of the corner's section at
//!   azimuth `b` and the wall's radius `w`, turned by `psi` in `[-back,
//!   ahead]` about the bend's axis. `|T_bb| = w`, `|T_b psi| <= w`,
//!   `|T_psi psi| <= R + w`. A cell is kept where the height above the
//!   plane may be `<= 0`: the height is 1-Lipschitz in space and the patch
//!   moves at most `w h_b + (R + w) h_psi` over a cell of half-widths `h`.
//! - The planar face, a square of the plane round the corner covering both
//!   sections, parameterised linearly. A cell is in the face where exactly
//!   one tube holds its points: a point is in a segment's tube when its
//!   distance from the axis line is within the disk's radii, and in a
//!   band's when its distance from the bend's circle is and its azimuth
//!   about the bend's axis lies within the band. Both distances are
//!   1-Lipschitz and the azimuth `1 / rho`-Lipschitz (`rho` the distance
//!   from the axis), so a cell whose centre decides every test with room
//!   for its half-diagonal is in or out for certain; any other cell is
//!   kept whole as a boundary cell.
//!
//! The certified value is reported under
//! [`super::MITRE_BESIDE_AN_ARC`], next to the proven budget for the rest.

use axiolid_construct::pipe::{MitreBand, MitreSide, PipePlan};
use axiolid_core::{Point3, Scalar, Vec3};
use axiolid_mesh::TriMesh;

use super::{Deviation, DeviationBound, DeviationPath, MITRE_BESIDE_AN_ARC};
use crate::certify::{certify_thin_quads, Cell, Coverage, Patch};
use crate::directrix::PiecePath;

/// Fraction of the target a quad may be twisted and still be covered as
/// one flat region; its thickness is added to the bound. A band's quads
/// twist more than a tube's (its vertices turn by different angles), and a
/// quarter of the target lets a sharp corner's band settle within the work
/// cap at 0.1 mm, where a twentieth does not.
const FLATNESS: Scalar = 0.25;

/// The deviation of a pipe swept along `path` with chord budget `built`,
/// reported against the requested `chord`.
pub(super) fn of_pipe(
    path: &PiecePath,
    radius: Scalar,
    inner: Option<Scalar>,
    fillet: Option<Scalar>,
    built: Scalar,
    chord: Scalar,
    mesh: &TriMesh,
) -> axiolid_contracts::GeomResult<Deviation> {
    let plan = axiolid_construct::pipe::pipe_plan(
        &path.pieces,
        path.closed,
        radius,
        inner,
        fillet,
        built,
    )?;
    let detail = if plan.closed {
        "closed loop of segments and arcs"
    } else {
        "segments and arcs"
    };
    let mut deviation = Deviation::one(
        DeviationPath::SweptDisk,
        detail,
        DeviationBound::Proven(chord),
    );
    if !plan.arc_mitres.is_empty() {
        let bound = certify_mitres(&plan, radius, inner, mesh, chord).map_or(
            DeviationBound::Unbounded("swept disk: a mitre beside an arc could not be certified"),
            DeviationBound::Certified,
        );
        deviation.add(DeviationPath::SweptDisk, MITRE_BESIDE_AN_ARC, bound);
    }
    Ok(deviation)
}

/// A certified bound on the distance from the exact surface near every
/// corner mitred beside an arc to `mesh`.
fn certify_mitres(
    plan: &PipePlan,
    radius: Scalar,
    inner: Option<Scalar>,
    mesh: &TriMesh,
    target: Scalar,
) -> Option<Scalar> {
    let walls: Vec<Scalar> = [Some(radius), inner].into_iter().flatten().collect();
    let mut bands: Vec<BandWall> = Vec::new();
    let mut faces: Vec<Face> = Vec::new();
    for mitre in &plan.arc_mitres {
        for side in &mitre.sides {
            if let MitreSide::Arc(band) = side {
                for &wall in &walls {
                    bands.push(BandWall::new(*band, wall));
                }
            }
        }
        faces.push(Face::new(mitre.sides, mitre.normal, radius, inner)?);
    }
    let mut patches: Vec<&dyn Patch> = Vec::new();
    let mut cells = Vec::new();
    for band in &bands {
        patches.push(band);
        let (back, ahead) = (band.band.back, band.band.ahead);
        grid(
            &mut cells,
            patches.len() - 1,
            (0.0, core::f64::consts::TAU),
            16,
            (-back, ahead),
            4,
        );
    }
    for face in &faces {
        patches.push(face);
        let d = face.half;
        grid(&mut cells, patches.len() - 1, (-d, d), 8, (-d, d), 8);
    }
    certify_thin_quads(&patches, cells, mesh, target, target, FLATNESS * target)
}

/// `nx` by `ny` cells over the box `x` by `y` of patch `patch`.
fn grid(
    cells: &mut Vec<Cell>,
    patch: usize,
    x: (Scalar, Scalar),
    nx: usize,
    y: (Scalar, Scalar),
    ny: usize,
) {
    let at =
        |r: (Scalar, Scalar), k: usize, n: usize| r.0 + (r.1 - r.0) * k as Scalar / n as Scalar;
    for i in 0..nx {
        for j in 0..ny {
            cells.push(Cell {
                patch,
                x: (at(x, i, nx), at(x, i + 1, nx)),
                y: (at(y, j, ny), at(y, j + 1, ny)),
            });
        }
    }
}

/// One wall (outer or bore) of an arc's band, up to the mitre plane.
struct BandWall {
    band: MitreBand,
    /// The wall's radius.
    wall: Scalar,
    /// Unit vector from the bend's centre to the corner.
    e: Vec3,
    /// The bend's radius.
    big: Scalar,
}

impl BandWall {
    fn new(band: MitreBand, wall: Scalar) -> Self {
        let big = band.bend_radius();
        Self {
            band,
            wall,
            e: (band.corner - band.centre) / big,
            big,
        }
    }

    /// The point of the corner's section at azimuth `b` on the wall.
    fn section(&self, b: Scalar) -> Point3 {
        let (sin, cos) = b.sin_cos();
        self.band.corner + (self.e * cos + self.band.axis * sin) * self.wall
    }
}

impl Patch for BandWall {
    fn jet(&self, b: Scalar, psi: Scalar) -> Option<(Point3, Vec3, Vec3)> {
        let (sin, cos) = b.sin_cos();
        let point = self.band.turned(self.section(b), psi);
        let along = (self.band.axis * cos - self.e * sin) * self.wall;
        let db = self.band.turned(self.band.centre + along, psi) - self.band.centre;
        let dpsi = self.band.axis.cross(point - self.band.centre);
        (point.is_finite() && db.is_finite()).then_some((point, db, dpsi))
    }

    fn second(&self, _b: (Scalar, Scalar), _psi: (Scalar, Scalar)) -> Option<[Scalar; 3]> {
        Some([self.wall, self.wall, self.big + self.wall])
    }

    fn coverage(&self, b: (Scalar, Scalar), psi: (Scalar, Scalar)) -> Coverage {
        let (bc, pc) = (0.5 * (b.0 + b.1), 0.5 * (psi.0 + psi.1));
        let height = self.band.height(self.band.turned(self.section(bc), pc));
        let moves = self.wall * 0.5 * (b.1 - b.0) + (self.big + self.wall) * 0.5 * (psi.1 - psi.0);
        if height + moves <= 0.0 {
            Coverage::Inside
        } else if height - moves > 0.0 {
            Coverage::Outside
        } else {
            Coverage::Boundary
        }
    }
}

/// The planar face of a mitre beside an arc: the points of the plane in
/// exactly one of the two tubes, inside the square of half-width `half`
/// round the incoming side's corner.
struct Face {
    origin: Point3,
    m1: Vec3,
    m2: Vec3,
    half: Scalar,
    sides: [MitreSide; 2],
    radii: (Scalar, Scalar),
}

impl Face {
    fn new(
        sides: [MitreSide; 2],
        normal: Vec3,
        radius: Scalar,
        inner: Option<Scalar>,
    ) -> Option<Self> {
        let corner = |side: &MitreSide| match side {
            MitreSide::Segment { corner, .. } => *corner,
            MitreSide::Arc(band) => band.corner,
            _ => Point3::splat(Scalar::NAN),
        };
        let origin = corner(&sides[0]);
        // Every section point lies within this of its own corner: a
        // segment's ellipse within `r / cos`, a band's within `r` of its
        // circle's points up to `back` or `ahead` from the corner.
        let mut half: Scalar = 0.0;
        for side in &sides {
            let extent = match side {
                MitreSide::Segment { direction, .. } => radius / direction.dot(normal).abs(),
                MitreSide::Arc(band) => {
                    radius + 2.0 * band.bend_radius() * (0.5 * band.back.max(band.ahead)).sin()
                }
                _ => return None,
            };
            half = half.max(extent + (corner(side) - origin).length());
        }
        let seed = if normal.x.abs() < 0.9 {
            Vec3::X
        } else {
            Vec3::Y
        };
        let m1 = (seed - normal * normal.dot(seed)).normalize();
        let m2 = normal.cross(m1);
        (half.is_finite() && half > 0.0).then_some(Self {
            origin,
            m1,
            m2,
            half: half * (1.0 + 1e-9),
            sides,
            radii: (inner.unwrap_or(0.0), radius),
        })
    }

    /// Whether every point within `d` of `p` is in the tube of `side`
    /// (`Some(true)`), none is (`Some(false)`), or it is not certain.
    fn holds(&self, side: &MitreSide, p: Point3, d: Scalar) -> Option<bool> {
        let (inner, outer) = self.radii;
        // Without a bore every distance down to the axis is in.
        let solid = inner <= 0.0;
        let within = |dist: Scalar| -> Option<bool> {
            if (solid || dist - d >= inner) && dist + d <= outer {
                Some(true)
            } else if (!solid && dist + d < inner) || dist - d > outer {
                Some(false)
            } else {
                None
            }
        };
        match side {
            MitreSide::Segment { corner, direction } => {
                let q = p - *corner;
                within((q - *direction * q.dot(*direction)).length())
            }
            MitreSide::Arc(band) => {
                let big = band.bend_radius();
                let e = (band.corner - band.centre) / big;
                let t = band.axis.cross(e);
                let q = p - band.centre;
                let h = q.dot(band.axis);
                let rho = (q - band.axis * h).length();
                let dist = (rho - big).hypot(h);
                let ring = within(dist);
                if ring == Some(false) {
                    return Some(false);
                }
                // The azimuth about the bend's axis from the corner's.
                let room = rho - d;
                if room <= 0.0 {
                    return None;
                }
                let azimuth = q.dot(t).atan2(q.dot(e));
                let spread = d / room;
                let window = if azimuth + spread < -band.back || azimuth - spread > band.ahead {
                    Some(false)
                } else if azimuth - spread >= -band.back && azimuth + spread <= band.ahead {
                    Some(true)
                } else {
                    None
                };
                match (ring, window) {
                    (_, Some(false)) => Some(false),
                    (Some(true), Some(true)) => Some(true),
                    _ => None,
                }
            }
            _ => None,
        }
    }
}

impl Patch for Face {
    fn jet(&self, x: Scalar, y: Scalar) -> Option<(Point3, Vec3, Vec3)> {
        Some((self.origin + self.m1 * x + self.m2 * y, self.m1, self.m2))
    }

    fn second(&self, _x: (Scalar, Scalar), _y: (Scalar, Scalar)) -> Option<[Scalar; 3]> {
        Some([0.0; 3])
    }

    fn coverage(&self, x: (Scalar, Scalar), y: (Scalar, Scalar)) -> Coverage {
        let centre = self.origin + self.m1 * (0.5 * (x.0 + x.1)) + self.m2 * (0.5 * (y.0 + y.1));
        let d = (0.5 * (x.1 - x.0)).hypot(0.5 * (y.1 - y.0));
        match (
            self.holds(&self.sides[0], centre, d),
            self.holds(&self.sides[1], centre, d),
        ) {
            (Some(a), Some(b)) if a != b => Coverage::Inside,
            (Some(_), Some(_)) => Coverage::Outside,
            _ => Coverage::Boundary,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Second partials of a band's wall by central differences, never
    /// above the patch's bounds.
    #[test]
    fn a_band_wall_bounds_its_second_partials() {
        // A leg along +x, then a bend of radius 0.3 heading -y: a right
        // angle beside the arc.
        let plan = axiolid_construct::pipe::pipe_plan(
            &[
                axiolid_construct::pipe::PathPiece::Segment {
                    start: Point3::new(-1.0, 0.0, 0.0),
                    end: Point3::ZERO,
                },
                axiolid_construct::pipe::PathPiece::Arc {
                    centre: Point3::new(0.3, 0.0, 0.0),
                    axis: Vec3::Z,
                    start: Point3::ZERO,
                    angle: 1.5,
                },
            ],
            false,
            0.05,
            None,
            None,
            1e-3,
        )
        .unwrap();
        assert_eq!(plan.arc_mitres.len(), 1);
        for mitre in &plan.arc_mitres {
            for side in &mitre.sides {
                let MitreSide::Arc(band) = side else { continue };
                let wall = BandWall::new(*band, 0.05);
                let [a, b, c] = wall.second((0.0, 1.0), (-0.1, 0.1)).unwrap();
                let h = 1e-6;
                for i in 0..=10 {
                    for j in 0..=10 {
                        let s = i as Scalar / 10.0;
                        let t = -0.1 + 0.2 * j as Scalar / 10.0;
                        let (_, xp, yp) = wall.jet(s + h, t).unwrap();
                        let (_, xm, ym) = wall.jet(s - h, t).unwrap();
                        let (_, _, yq) = wall.jet(s, t + h).unwrap();
                        let (_, _, yr) = wall.jet(s, t - h).unwrap();
                        assert!(((xp - xm) / (2.0 * h)).length() <= a * 1.0001);
                        assert!(((yp - ym) / (2.0 * h)).length() <= b * 1.0001);
                        assert!(((yq - yr) / (2.0 * h)).length() <= c * 1.0001);
                    }
                }
            }
        }
    }
}
