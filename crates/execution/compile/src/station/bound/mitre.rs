//! The bound of a cell beside a mitre (#335, ADR 0082 amendment).

use axiolid_contracts::GeomResult;
use axiolid_core::{Point3, Scalar, Vec3};
use axiolid_model::SeamSide;
use axiolid_reference::station::{FrameLawBounds, FrameReading, Mitre};

use super::super::seam_slack;
use super::{CellBound, SectionLaw, UNCERTIFIED_FRAMES};

/// What `SectionLaw::cell` knows of a cell beside a mitre.
pub(super) struct MitreInputs<'a> {
    pub(super) lo: Scalar,
    pub(super) hi: Scalar,
    pub(super) at_lo: &'a [Point3],
    pub(super) at_hi: &'a [Point3],
    /// Per vertex, the chord, breakpoint and read terms.
    pub(super) eta: &'a [Scalar],
    pub(super) slide: &'a [Scalar],
    /// The frame bounds of the pieces next to the cell's ends.
    pub(super) first: &'a FrameLawBounds,
    pub(super) last: &'a FrameLawBounds,
    pub(super) reading: FrameReading,
    /// The longitudinal offset at the cell's ends.
    pub(super) long: [Scalar; 2],
}

/// One side of a mitre, as the cell beside it reads it.
struct SideCut {
    /// The point of the mitre plane the exact solid is cut at: the seam
    /// point moved by the longitudinal offset along the normal.
    q0: Point3,
    /// The plane's normal, pointing out of the cell.
    normal: Vec3,
    /// The cell's end tangent, pointing at the plane.
    toward: Vec3,
    cos: Scalar,
    /// The cell's own exact section at the knot.
    exact: Vec<Point3>,
    /// That section extruded along `toward` onto the plane.
    projected: Vec<Point3>,
    /// Per vertex, the mesh's mitre row's distance from `projected`.
    gap: Vec<Scalar>,
    /// The mesh's mitre row's largest distance from the plane.
    off_plane: Scalar,
    /// How far the plane's point may lie from the exact one.
    point_error: Scalar,
    /// How far the end tangent may lie from the exact one.
    tangent_error: Scalar,
    /// How far the plane's normal may lie from the exact one.
    normal_error: Scalar,
}

impl SideCut {
    fn height(&self, x: Point3) -> Scalar {
        (x - self.q0).dot(self.normal)
    }

    /// The bound on how far the plane read from the frames lies from the
    /// exact one at the points involved: `far` is the cell's other row.
    fn plane_error(&self, far: &[Point3], eta_max: Scalar) -> Scalar {
        let reach = far
            .iter()
            .chain(&self.exact)
            .chain(&self.projected)
            .map(|x| (*x - self.q0).length())
            .fold(0.0, Scalar::max)
            + eta_max;
        self.point_error + reach * self.normal_error
    }

    fn lift(&self) -> Scalar {
        self.exact
            .iter()
            .map(|c| self.height(*c).abs())
            .fold(0.0, Scalar::max)
    }
}

impl SectionLaw<'_, '_> {
    /// The frame bounds over a span of `length` on the other side of the
    /// knot `at`, clipped to the next seam, station or run end.
    fn across(
        &self,
        at: Scalar,
        side: SeamSide,
        length: Scalar,
        reading: FrameReading,
    ) -> Option<FrameLawBounds> {
        let mut marks: Vec<Scalar> = self
            .seams
            .iter()
            .chain(self.distances.iter())
            .copied()
            .collect();
        marks.retain(|m| (m - at).abs() > seam_slack(at));
        let (first, last) = (self.distances[0], self.distances[self.distances.len() - 1]);
        match side {
            // The cell is incoming: the other side starts at the knot.
            SeamSide::Incoming => {
                let end = marks
                    .iter()
                    .copied()
                    .filter(|m| *m > at)
                    .fold((at + length).min(last), Scalar::min);
                self.frames.bounds(at, end, reading)
            }
            _ => {
                let start = marks
                    .iter()
                    .copied()
                    .filter(|m| *m < at)
                    .fold((at - length).max(first), Scalar::max);
                self.frames.bounds(start, at, reading)
            }
        }
    }

    /// The cell's side of `mitre`, `side` the cell's own: `Err` names why
    /// the frames beside it are not bounded.
    fn side_cut(
        &self,
        inputs: &MitreInputs<'_>,
        mitre: &Mitre,
        side: SeamSide,
    ) -> GeomResult<Result<SideCut, &'static str>> {
        let incoming = side == SeamSide::Incoming;
        let (at, row, near, long) = if incoming {
            (inputs.hi, inputs.at_hi, inputs.last, inputs.long[1])
        } else {
            (inputs.lo, inputs.at_lo, inputs.first, inputs.long[0])
        };
        let (toward, normal, base) = if incoming {
            (mitre.incoming.tangent, mitre.normal, mitre.incoming.point)
        } else {
            (-mitre.outgoing.tangent, -mitre.normal, mitre.outgoing.point)
        };
        let Some(other) = self.across(at, side, inputs.hi - inputs.lo, inputs.reading) else {
            return Ok(Err(UNCERTIFIED_FRAMES));
        };
        let cos = toward.dot(normal);
        let tangent_error = near.axis_error + near.measure * near.tangent[0];
        let other_error = other.axis_error + other.measure * other.tangent[0];
        let mut cut = SideCut {
            // The shift is along the mitre's own normal, whichever side.
            q0: base + mitre.normal * long,
            normal,
            toward,
            cos,
            exact: self.exact_at(at, side)?,
            projected: Vec::new(),
            gap: Vec::new(),
            off_plane: 0.0,
            point_error: near.point_error + near.measure * near.point[0],
            tangent_error,
            normal_error: (tangent_error + other_error) / cos,
        };
        cut.projected = cut
            .exact
            .iter()
            .map(|&c| c - cut.toward * (cut.height(c) / cut.cos))
            .collect();
        cut.gap = row
            .iter()
            .zip(&cut.projected)
            .map(|(m, p)| (*m - *p).length())
            .collect();
        cut.off_plane = row
            .iter()
            .map(|m| cut.height(*m).abs())
            .fold(0.0, Scalar::max);
        Ok(Ok(cut))
    }

    /// A cell on a straight piece with a constant section, beside one or
    /// both of its mitres: every vertex runs on a line through the plane,
    /// so the exact wall is the quad up to the plane and the mesh's mitre
    /// row differs from it by the rows' gap alone.
    pub(super) fn straight_mitre_cell(
        &self,
        inputs: &MitreInputs<'_>,
        lo: Option<&Mitre>,
        hi: Option<&Mitre>,
    ) -> GeomResult<Result<CellBound, &'static str>> {
        let count = inputs.at_lo.len();
        let mut gap = vec![0.0 as Scalar; count];
        let (mut off_plane, mut plane_error, mut extension): (Scalar, Scalar, Scalar) =
            (0.0, 0.0, 0.0);
        let eta_max = inputs.eta.iter().copied().fold(0.0, Scalar::max);
        for (mitre, side, far) in [
            (lo, SeamSide::Outgoing, inputs.at_hi),
            (hi, SeamSide::Incoming, inputs.at_lo),
        ] {
            let Some(mitre) = mitre else { continue };
            let cut = match self.side_cut(inputs, mitre, side)? {
                Ok(cut) => cut,
                Err(reason) => return Ok(Err(reason)),
            };
            for (g, d) in gap.iter_mut().zip(&cut.gap) {
                *g = g.max(*d);
            }
            off_plane = off_plane.max(cut.off_plane);
            plane_error = plane_error.max(cut.plane_error(far, eta_max));
            extension = extension.max(cut.lift() / cut.cos * cut.tangent_error);
        }
        let mut worst: Scalar = 0.0;
        for &(k, l) in self.edges {
            let quad = axiolid_construct::loft::quad_deviation(
                inputs.at_lo[k],
                inputs.at_lo[l],
                inputs.at_hi[k],
                inputs.at_hi[l],
            );
            let gap = gap[k].max(gap[l]);
            let slide = inputs.slide[k].max(inputs.slide[l]) + extension;
            worst = worst.max(slide + (quad + gap).max(gap + off_plane + plane_error));
        }
        let floor = gap.iter().copied().fold(0.0, Scalar::max);
        let worst = worst * (1.0 + 1e-12);
        Ok(if worst.is_finite() {
            Ok(CellBound {
                value: worst,
                floor,
            })
        } else {
            Err(UNCERTIFIED_FRAMES)
        })
    }
}
