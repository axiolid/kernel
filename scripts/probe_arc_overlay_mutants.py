"""Mutation probe for the exact arc overlay: each mutant must be caught."""
import pathlib, shutil, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
SRC = ROOT / "crates/algorithms/planar/overlay/src"

MUTANTS = [
    ("turn rank inverted", "exact_arc.rs",
     "        Sign::Positive => 0,\n        Sign::Negative => 2,",
     "        Sign::Positive => 2,\n        Sign::Negative => 0,"),
    ("tie-break inverted", "exact_arc.rs",
     "                        }) == Sign::Positive\n                    }\n                    _ => false,",
     "                        }) == Sign::Negative\n                    }\n                    _ => false,"),
    ("half-open rule dropped", "exact_arc.rs",
     "        let up = ya != Sign::Positive && yb == Sign::Positive;",
     "        let up = ya == Sign::Negative && yb == Sign::Positive;"),
    ("output rounding cleanup disabled", "arc_overlay.rs",
     "            (b - a).length() <= tolerance.linear()\n        });",
     "            (b - a).length() < 0.0\n        });"),
    ("shared pieces kept twice", "exact_arc.rs",
     "(Op::Union | Op::Intersection, S::SharedSame) if subject => Some(false),",
     "(Op::Union | Op::Intersection, S::SharedSame) => Some(false),"),
    ("arc side ignored in crossing count", "exact_arc.rs",
     "    if chord != bulge && chord != Sign::Zero {\n        return chord;\n    }",
     "    if true {\n        return chord;\n    }"),
    ("box shortcut in point identity unsound", "exact_arc/point.rs",
     "    if a.bx.disjoint(b.bx) || a.by.disjoint(b.by) {\n        return false;\n    }",
     "    if !a.bx.disjoint(b.bx) {\n        return false;\n    }"),
    ("difference keeps clip outside", "exact_arc.rs",
     "        (Op::Difference, S::Inside) if !subject => Some(true),",
     "        (Op::Difference, S::Outside) if !subject => Some(true),"),
    # Broad phase (bounding boxes): a box that misses part of its edge
    # skips a real crossing or a real shared edge.
    ("arc box ignores the bulge", "exact_arc/edge.rs",
     "    Bounds::around(&[from, to], bulge.abs() * chord / 2.0)",
     "    Bounds::around(&[from, to], 0.0)"),
    ("arc box uses the half sagitta", "exact_arc/edge.rs",
     "    Bounds::around(&[from, to], bulge.abs() * chord / 2.0)",
     "    Bounds::around(&[from, to], bulge.abs() * chord / 4.0)"),
    ("box overlap test wrong on y", "exact_arc/edge.rs",
     "self.y0 <= other.y1 && other.y0 <= self.y1",
     "self.y0 <= other.y1 && other.y0 <= self.y0"),
    ("link index window too narrow", "exact_arc.rs",
     "            .take_while(|entry| entry.0 <= hi)",
     "            .take_while(|entry| entry.0 < lo)"),
]

def run():
    return subprocess.run(
        ["cargo", "test", "-q", "--release", "-p", "axiolid-overlay",
         "--test", "arc_exact_oracle", "--test", "arc_overlay"],
        cwd=ROOT, capture_output=True, text=True, timeout=900,
    ).returncode

survivors = []
for name, rel, old, new in MUTANTS:
    path = SRC / rel
    original = path.read_text()
    assert original.count(old) == 1, f"anchor for '{name}' not unique/found"
    path.write_text(original.replace(old, new))
    try:
        try:
            code = run()
        except subprocess.TimeoutExpired:
            code = -1  # a hang is a detection
    finally:
        path.write_text(original)
    status = "killed" if code != 0 else "SURVIVED"
    print(f"{status:8} {name}", flush=True)
    if code == 0:
        survivors.append(name)
print(f"{len(MUTANTS) - len(survivors)}/{len(MUTANTS)} killed")
sys.exit(1 if survivors else 0)

