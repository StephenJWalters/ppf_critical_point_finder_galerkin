# Critical Point Finder for Plane Poiseuille Flow (Orr–Sommerfeld)

A Fortran program that computes the **critical Reynolds number**, **critical wavenumber**, and corresponding wavespeed for the linear stability of plane Poiseuille flow, by solving the Orr–Sommerfeld equation with a Chebyshev spectral (NSC, "no-slip Chebyshev") method and locating the neutral point via 2D Newton iteration in `(Re, alpha)`, using a fully analytic Jacobian.

At the critical point the flow transitions from linearly stable to linearly unstable. This code starts with the classical result (Orszag 1971):

```
Re_c    ≈ 5772.22
alpha_c ≈ 1.02056
c_r     ≈ 0.264
```
and finds the true critical parameters with at least 28 significant digits.

## Background

The Orr–Sommerfeld equation governs the linear stability of parallel shear flows. Discretising it with a Chebyshev-Galerkin (NSC) basis that satisfies the no-slip boundary conditions exactly turns the eigenvalue problem into a generalized matrix eigenvalue problem

```
A(alpha, Re) q = lambda B(alpha) q,      c = lambda / Re
```

The critical point is where the growth rate `ci = Im(c)` and its alpha-derivative `dci/dalpha` simultaneously vanish — a 2D root-finding problem in `(Re, alpha)` solved here with Newton's method.

### Analytic Jacobian

Because `A` is polynomial in `alpha` and linear in `Re`, and `B` is polynomial in `alpha` only, every derivative the Newton iteration needs (`dA/dalpha`, `dA/dRe`, `d²A/dalpha²`, `d²A/dalpha dRe`, `dB/dalpha`, `d²B/dalpha²`) is available in closed form — no finite differencing is required. Standard eigenvalue-perturbation theory then gives:

- **First derivatives** of the eigenvalue via the left/right eigenvector pair:
  `dlambda/dp = w^H (dA/dp − lambda dB/dp) q / (w^H B q)`
- **Second derivatives** via the "bordered system" technique from bifurcation/continuation   theory: `dq/dalpha` and `dw/dalpha` are obtained by solving the singular systems `(A − lambda B) dq/dalpha = RHS` and its adjoint, each bordered with one extra row/column to remove the null-space ambiguity. The result is provably independent of which particular solution of the singular system is chosen.

This replaces an earlier version that used a 9-point finite-difference stencil (9 full eigenvalue solves per Newton step) with 1 direct eigensolve, 1 adjoint eigensolve, and 2 small bordered linear solves — converging to the same critical point to ~13 significant figures, in fewer iterations, and substantially faster.

All arithmetic is done in **quadruple precision** (`selected_real_kind(33)`) to support extremely tight convergence tolerances and to cleanly separate genuine truncation error (from the number of Chebyshev modes `N`) from floating-point roundoff.

## Files

| File | Description |
|---|---|
| `critical_point_finder_galerkin.f90` | Main program and module (`finder_mod`) |

## Requirements

- `gfortran` (or any Fortran 2008-compatible compiler)
- No external libraries required (LU factorisation/solve routines are self-contained)

## Building

```bash
gfortran critical_point_finder_galerkin.f90 -O3 -o finder
```

## Running

```bash
./finder
```

The program runs the solver at several resolutions `N = 256, 32, 64, 128, 144, 192` (Chebyshev modes) — an initial high-resolution reference run followed by a spectral convergence sweep — and prints, for each `N`:

```
N   Re_c   alpha_c   c_r
```

to stdout, along with per-run CPU timings.

It also writes a LaTeX convergence table to `convergence_table.tex`, showing the digits of `Re_c` that agree with the reference (`N = 256`) solution as `N` increases — demonstrating the spectral (exponential) convergence rate of the method.

## Adjusting resolution / precision

- `Nmax`, `ny` (Chebyshev modes / grid points) are set per-run in the `nmax_vec` array in the main program.
- Newton tolerance (`tol`) and max iterations (`max_iter`) are set near the top of the `finder` subroutine.
- Working precision is controlled by the `qp` parameter (`selected_real_kind(33)`) at the top of `finder_mod`.

## Citation

If you use this code in academic work, please cite the accompanying paper:

> *[Add paper title, authors, journal/arXiv reference here]*

## License

MIT
