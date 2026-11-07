//! A static tree of boxes, for the subdivision's broad phase (#173).
//!
//! Only skips work: a query returns every item whose box passes the test
//! (and possibly items near it), and every decision after it is exact. A
//! box with a NaN coordinate cannot be sorted or tested, so it is returned
//! by every query instead.

use super::edge::Bounds;

/// Items per leaf.
const LEAF: usize = 8;

struct Node {
    bounds: Bounds,
    /// Leaf: `items[start..end]`. Inner: children at `start` and `end`.
    start: usize,
    end: usize,
    leaf: bool,
}

pub(crate) struct BoxTree {
    nodes: Vec<Node>,
    /// Item indices, grouped by leaf.
    items: Vec<usize>,
    boxes: Vec<Bounds>,
    /// Items whose box cannot be tested.
    always: Vec<usize>,
}

impl BoxTree {
    pub(crate) fn new(boxes: Vec<Bounds>) -> Self {
        let (mut items, always): (Vec<usize>, Vec<usize>) =
            (0..boxes.len()).partition(|&i| !boxes[i].has_nan());
        let mut nodes = Vec::new();
        if !items.is_empty() {
            let count = items.len();
            build(&boxes, &mut items, 0, count, &mut nodes);
        }
        Self {
            nodes,
            items,
            boxes,
            always,
        }
    }

    /// Every item whose box passes `test`, ascending. `test` must pass a
    /// node's box whenever it passes a box inside it.
    ///
    /// The items go into `out`, replacing what it held, with `stack` as
    /// scratch: callers asking once per piece reuse both (#218).
    pub(crate) fn query_into(
        &self,
        test: impl Fn(&Bounds) -> bool,
        out: &mut Vec<usize>,
        stack: &mut Vec<usize>,
    ) {
        out.clear();
        out.extend_from_slice(&self.always);
        if !self.nodes.is_empty() {
            stack.clear();
            stack.push(0);
            while let Some(index) = stack.pop() {
                let node = &self.nodes[index];
                if !test(&node.bounds) {
                    continue;
                }
                if node.leaf {
                    out.extend(
                        self.items[node.start..node.end]
                            .iter()
                            .copied()
                            .filter(|&item| test(&self.boxes[item])),
                    );
                } else {
                    stack.push(node.start);
                    stack.push(node.end);
                }
            }
        }
        // Ascending, so callers see items in a fixed order.
        out.sort_unstable();
    }

    /// The first item whose box passes `test` and which `accept` takes, in
    /// a fixed order (the untestable items, then the leaves depth-first,
    /// left to right), or `None` when no item does
    /// (#344). Stops at the first, so a search for any witness costs the
    /// path to it, not every candidate. `test` must pass a node's box
    /// whenever it passes a box inside it.
    pub(crate) fn first(
        &self,
        test: impl Fn(&Bounds) -> bool,
        mut accept: impl FnMut(usize) -> bool,
    ) -> Option<usize> {
        if let Some(&item) = self.always.iter().find(|&&item| accept(item)) {
            return Some(item);
        }
        if self.nodes.is_empty() {
            return None;
        }
        let mut stack = vec![0];
        while let Some(index) = stack.pop() {
            let node = &self.nodes[index];
            if !test(&node.bounds) {
                continue;
            }
            if node.leaf {
                let found = self.items[node.start..node.end]
                    .iter()
                    .copied()
                    .find(|&item| test(&self.boxes[item]) && accept(item));
                if found.is_some() {
                    return found;
                }
            } else {
                stack.push(node.end);
                stack.push(node.start);
            }
        }
        None
    }
}

/// Build the subtree over `items[start..end]`; returns its node index.
fn build(
    boxes: &[Bounds],
    items: &mut [usize],
    start: usize,
    end: usize,
    nodes: &mut Vec<Node>,
) -> usize {
    let bounds = Bounds::hull(items[start..end].iter().map(|&i| boxes[i]))
        .expect("a subtree holds at least one item");
    let index = nodes.len();
    nodes.push(Node {
        bounds,
        start,
        end,
        leaf: true,
    });
    if end - start <= LEAF {
        return index;
    }
    // Split at the median centre along the longer side.
    let wide = bounds.wider_than_tall();
    let centre = |i: usize| {
        let (x, y) = boxes[i].centre();
        if wide {
            x
        } else {
            y
        }
    };
    let mid = start + (end - start) / 2;
    items[start..end].select_nth_unstable_by(mid - start, |&a, &b| centre(a).total_cmp(&centre(b)));
    let left = build(boxes, items, start, mid, nodes);
    let right = build(boxes, items, mid, end, nodes);
    nodes[index] = Node {
        bounds,
        start: left,
        end: right,
        leaf: false,
    };
    index
}

#[cfg(test)]
mod tests {
    use axiolid_core::Point2;

    use super::*;

    fn at(x: f64, y: f64) -> Bounds {
        Bounds::around(&[Point2::new(x, y)], 0.0)
    }

    #[test]
    fn the_first_box_query_stops_at_the_first_accepted_item() {
        let boxes: Vec<Bounds> = (0..200).map(|i| at(f64::from(i), 0.0)).collect();
        let tree = BoxTree::new(boxes);
        let probe = Bounds::around(&[Point2::new(-1.0, -1.0), Point2::new(300.0, 1.0)], 0.0);
        // Every box passes: the search must stop at the first accepted item.
        let mut asked = 0;
        let found = tree.first(
            |b| b.overlaps(&probe),
            |item| {
                asked += 1;
                item % 7 == 3
            },
        );
        let item = found.expect("an accepted item");
        assert_eq!(item % 7, 3);
        assert!(asked < 200, "{asked} items were asked");
        // Nothing accepted: every box-passing item is asked, and none found.
        let mut seen = Vec::new();
        assert_eq!(
            tree.first(
                |b| b.overlaps(&probe),
                |item| {
                    seen.push(item);
                    false
                }
            ),
            None
        );
        seen.sort_unstable();
        assert_eq!(seen, (0..200).collect::<Vec<_>>());
        // The box test prunes: only the items near 50 are ever asked.
        let near = Bounds::around(&[Point2::new(49.5, -1.0), Point2::new(51.5, 1.0)], 0.0);
        let mut asked = Vec::new();
        assert_eq!(
            tree.first(
                |b| b.overlaps(&near),
                |item| {
                    asked.push(item);
                    false
                }
            ),
            None
        );
        asked.sort_unstable();
        assert_eq!(asked, vec![50, 51]);
    }
}
