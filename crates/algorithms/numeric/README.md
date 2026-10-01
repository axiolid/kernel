# axiolid-numeric

General numeric routines that geometry algorithms build on: bracketed root
finding (Brent) and every real root of a polynomial in an interval,
adaptive Gauss-Kronrod quadrature and fixed Gauss-Legendre rules, LU,
Cholesky and column-pivoted QR solves, linear least squares with and
without equality constraints, Levenberg-Marquardt for nonlinear least
squares, and Brent's bounded scalar minimiser. Each result carries a
bracket, an error estimate or a condition estimate, and each iterative
routine says whether it met its tolerance. Non-finite input, singular or
rank-deficient matrices and unbracketed roots are refused by name. Values
are plain `f64`; there are no geometry types and no dependencies.

```bash
cargo add axiolid-numeric
```

- API documentation: [docs.rs/axiolid-numeric](https://docs.rs/axiolid-numeric)
- Reference page: [axiolid.github.io/kernel](https://axiolid.github.io/kernel/reference/crates/axiolid-numeric)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)

## Design notes

Polynomial roots are isolated through the roots of the derivative, and a
sign is trusted only where the value exceeds a rigorous bound on its own
rounding error. Where it does not (a multiple root, a cluster closer than
`f64` separates) the answer is an unresolved region with the most roots it
can hold, not a guessed count. For exact decisions on integer polynomials
use `axiolid-exact`.

Older inline routines elsewhere in the workspace (bisection and an 8-point
Gauss-Legendre rule in `axiolid-evaluate`, the interpolation solve in
`axiolid-nurbs`) are unchanged; new code should use this crate.
