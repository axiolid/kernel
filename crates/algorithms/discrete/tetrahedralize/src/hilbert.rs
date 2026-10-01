//! Insertion order along a 3D Hilbert curve.
//!
//! Incremental insertion locates each point by walking from the cell the
//! previous insertion created. Consecutive points that are close in space
//! keep that walk short; a Hilbert order makes them so. The order affects
//! only speed: the perturbed triangulation is unique.

use axiolid_core::Point3;

/// Bits per axis of the quantised grid: 3 * 21 = 63 bits of curve index.
const BITS: u32 = 21;

/// Indices of `points` (restricted to `subset`) in Hilbert-curve order.
pub(crate) fn order(points: &[Point3], subset: &[u32]) -> Vec<u32> {
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for &i in subset {
        let p = points[i as usize];
        for (axis, value) in [p.x, p.y, p.z].into_iter().enumerate() {
            lo[axis] = lo[axis].min(value);
            hi[axis] = hi[axis].max(value);
        }
    }
    let extent = (0..3).map(|axis| hi[axis] - lo[axis]).fold(0.0, f64::max);
    let cells = f64::from((1u32 << BITS) - 1);
    let scale = if extent > 0.0 && extent.is_finite() {
        cells / extent
    } else {
        0.0
    };
    let quantise = |value: f64, axis: usize| {
        // Saturating float-to-int conversion keeps this total even if the
        // scaled value rounds past the grid.
        ((value - lo[axis]) * scale).clamp(0.0, cells) as u32
    };
    let mut keyed: Vec<(u64, u32)> = subset
        .iter()
        .map(|&i| {
            let p = points[i as usize];
            let key = index([quantise(p.x, 0), quantise(p.y, 1), quantise(p.z, 2)]);
            (key, i)
        })
        .collect();
    keyed.sort_unstable();
    keyed.into_iter().map(|(_, i)| i).collect()
}

/// Hilbert index of a grid cell (Skilling, "Programming the Hilbert curve",
/// 2004): transform the axes to the curve's transposed form, then
/// interleave the bits.
fn index(mut x: [u32; 3]) -> u64 {
    let top = 1u32 << (BITS - 1);
    let mut q = top;
    while q > 1 {
        let p = q - 1;
        for i in 0..3 {
            if x[i] & q != 0 {
                x[0] ^= p;
            } else {
                let t = (x[0] ^ x[i]) & p;
                x[0] ^= t;
                x[i] ^= t;
            }
        }
        q >>= 1;
    }
    x[1] ^= x[0];
    x[2] ^= x[1];
    let mut t = 0;
    q = top;
    while q > 1 {
        if x[2] & q != 0 {
            t ^= q - 1;
        }
        q >>= 1;
    }
    for value in &mut x {
        *value ^= t;
    }
    let mut key = 0u64;
    for bit in (0..BITS).rev() {
        for value in x {
            key = (key << 1) | u64::from((value >> bit) & 1);
        }
    }
    key
}

#[cfg(test)]
mod tests {
    use super::index;

    /// Consecutive Hilbert indices are face-adjacent grid cells: the property
    /// that makes the order local.
    #[test]
    fn consecutive_indices_are_adjacent_cells() {
        let n = 8u32;
        let mut cells: Vec<(u64, [u32; 3])> = (0..n)
            .flat_map(|x| (0..n).flat_map(move |y| (0..n).map(move |z| [x, y, z])))
            .map(|c| (index(c), c))
            .collect();
        cells.sort_unstable();
        for pair in cells.windows(2) {
            let (a, b) = (pair[0].1, pair[1].1);
            let distance: u32 = (0..3).map(|i| a[i].abs_diff(b[i])).sum();
            assert_eq!(distance, 1, "{a:?} -> {b:?}");
        }
    }
}
