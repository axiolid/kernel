# Parameterised section lowering

Status: accepted
Date: 2026-09-12

## Context

`Profile::Section` was the last refusing profile variant, and unlike the
others it had NO prior art: it was refused on the exact path and on the
tessellating path alike. It carries nine variants of parameterised structural
section -- I, asymmetric I, L, T, U, C, Z, trapezium -- each with optional
root fillets, edge radii and tapers.

## Decision

Lower each variant to a closed counter-clockwise ring of CORNERS, each
carrying an optional radius, and route that ring through one shared function
that inserts a tangent arc at every rounded corner. The result is a
`ContourProfile`, which the contour path (ADR 0053) already extrudes exactly,
including genuine cylindrical walls for the fillets.

Concave root fillets and convex toe radii are the SAME operation: the sign of
the turn decides which way the arc bends. Giving them separate code paths
would let them drift apart.

## Why the fillets are built, not dropped

For a rolled steel section the web-to-flange root fillet is real material.
Measured on a 0.4 x 0.3 I-section with an 0.021 root radius, the four fillets
carry **2.40%** of the cross-sectional area:

    area without fillets   0.015382000
    area with fillets      0.015760558

Each fillet adds `r^2 (1 - pi/4)` over the sharp corner it replaces. That
term was verified by Monte-Carlo integration to 1.1e-4 relative before being
used, and the assembled outline then matched the closed form to 9e-9.

Dropping the fillets yields a section that still looks like an I and whose
area, second moment and mass are all wrong.

## Consequences

- `Profile::Section` extrudes exactly. All eight profile variants now do.
- A stated TAPER is refused by name. A sloped flange moves the fillet
  tangency onto an inclined face, which is a different construction; building
  the parallel-flange outline anyway would silently return the wrong section.
  `Some(0.0)` is a declared parallel flange and is accepted -- only a
  non-zero slope is a taper.
- `None` and `Some(0.0)` stay distinct where the source distinguishes them:
  an absent L width means an EQUAL angle, an absent top flange thickness
  means the bottom value. Both are mutation-tested.
- `C` is treated as thin-walled: the boundary follows the wall all the way
  round including the lips, so the enclosed area is material rather than the
  channel envelope.
- Radii that do not fit between two corners sharing an edge are refused,
  reusing the pairwise check the fillet and chamfer paths already use.

## Findings

### A right-angle-only bug, found by testing the router's generality

The arc centre was placed at `r / sin(turn / 2)` from the corner. The correct
distance is `r / cos(turn / 2)`: `turn` is the EXTERIOR deflection, so the
interior angle is `pi - turn` and the centre sits `r / sin(interior / 2)`
away.

At a right-angle corner `sin(45) == cos(45)`, so the two agree EXACTLY. Every
rounded corner reachable through `SectionProfile` today is a right angle, so
all sixteen area tests passed against the wrong formula, and the
corresponding mutant SURVIVED.

It was caught by testing the corner router directly at a 116.565-degree turn,
where the arc came out with radius 0.01748 instead of 0.02. The invariant
asserted is TANGENCY -- the distance from the arc centre to each adjacent
edge line must equal the radius -- which holds at any corner angle and pins
the formula rather than the symptom.

The lesson: when a shared routine is general but every current caller
exercises one special case, test the routine's generality directly. The area
tests could not have found this, however many variants they covered.

### Closed-form areas need an exact integrator

`exact_properties` is planar-only, so a filleted section cannot go through
it. Areas are computed by Green's theorem over the exact contour, with a
circular arc contributing its closed-form `r^2 (sweep - sin sweep)` term.
Chord-sampling the arcs there would have made the fillet tests measure the
approximation rather than the geometry.

## Mutation evidence

8/8 killed:

- fillets dropped entirely
- arc bends the wrong way (convex root fillet)
- setback ignores the turn angle
- centre distance uses sin (the right-angle-only bug above)
- absent L width becomes zero instead of an equal angle
- absent top flange thickness becomes zero
- taper silently built parallel
- oversized radius check removed
