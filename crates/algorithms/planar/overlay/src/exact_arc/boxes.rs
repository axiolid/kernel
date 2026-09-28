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
    pub(crate) fn query(&self, test: impl Fn(&Bounds) -> bool) -> Vec<usize> {
        let mut out = self.always.clone();
        if !self.nodes.is_empty() {
            let mut stack = vec![0];
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
        out
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
