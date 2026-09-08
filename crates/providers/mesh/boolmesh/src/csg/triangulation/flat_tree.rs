//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

use crate::csg::triangulation::Pt;
use crate::csg::{Real, Vec2};

pub fn compute_flat_tree(pts: &mut [Pt]) {
    if pts.len() <= 8 {
        return;
    }
    compute_flat_tree_impl(pts, true);
}

fn compute_flat_tree_impl(pts: &mut [Pt], sort_x: bool) {
    let eq = std::cmp::Ordering::Equal;
    if sort_x {
        pts.sort_by(|a, b| a.pos.x.partial_cmp(&b.pos.x).unwrap_or(eq));
    } else {
        pts.sort_by(|a, b| a.pos.y.partial_cmp(&b.pos.y).unwrap_or(eq));
    }

    if pts.len() < 2 {
        return;
    }

    let (l, mr) = pts.split_at_mut(pts.len() / 2);
    if !mr.is_empty() {
        let (_, r) = mr.split_first_mut().unwrap();
        compute_flat_tree_impl(l, !sort_x);
        compute_flat_tree_impl(r, !sort_x);
    }
}

pub fn compute_query_flat_tree<F>(pts: &[Pt], rect: &Rect, mut func: F)
where
    F: FnMut(&Pt),
{
    for p in pts.iter() {
        if rect.contains(&p.pos) {
            func(p);
        }
    }

    //if pts.len() <= 8 {
    //    for p in pts.iter() { if rect.contains(&p.pos) { func(p);} }
    //} else {
    //    query_two_d_tree(pts, rect.clone(), func);
    //}
}

#[derive(Clone)]
pub struct Rect {
    pub min: Vec2,
    pub max: Vec2,
}

impl Rect {
    pub fn default() -> Self {
        Self {
            min: Vec2::MAX,
            max: Vec2::MIN,
        }
    }

    pub fn new(a: &Vec2, b: &Vec2) -> Self {
        Self {
            min: a.min(*b),
            max: a.max(*b),
        }
    }

    pub fn contains(&self, p: &Vec2) -> bool {
        p.x >= self.min.x && p.x <= self.max.x && p.y >= self.min.y && p.y <= self.max.y
    }

    pub fn size(&self) -> Vec2 {
        self.max - self.min
    }

    pub fn scale(&self) -> Real {
        let a_min = self.min.x.abs().max(self.min.y.abs());
        let a_max = self.max.x.abs().max(self.max.y.abs());
        a_min.max(a_max)
    }

    pub fn union(&mut self, p: Vec2) {
        self.min = self.min.min(p);
        self.max = self.max.max(p);
    }
}
