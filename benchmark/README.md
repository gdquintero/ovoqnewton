# Numerical experiments of the revised manuscript

Code, data and results of Section 4 of *A regularized Newton-type method for
box-constrained order-value optimization*.

## Requirements

- GNU Fortran, LAPACK/BLAS.
- [Algencan 3.1.x](https://www.ime.usp.br/~egbirgin/tango/) with HSL, installed in
  `$HOME/algencan-3.1.2` (or set `ALGENCAN` when calling `make`).
- Python 3 (standard library) for the scripts that run the experiments; NumPy and
  Matplotlib for the analysis and the figures.

## Code (`src/`)

- `ovo_bench.f90`: Algorithm 2.1 with a multistart over perturbed starting points.
  Curvature families `B=0|H|GN|QN` (first-order, shifted Hessian, Gauss-Newton,
  shared BFGS), fixed bound `M=` with `mmode=scale|clip`, adaptive bound `cad=c`
  (`||B|| <= c*sigma`), shift margin `tau=`, and the stopping test of Step 5. All
  options are documented in the header of the file.
- `hyperdual.f90`: forward-mode automatic differentiation of second order, used to
  obtain exact gradients and Hessians of the residuals.
- `models.f90`: the twelve models (polynomial, Bard and ten NIST problems).
- `test_derivs.f90`: checks the derivatives against finite differences.

Build with `cd src && make`. The linker warning about an executable stack is
expected: the Algencan callbacks are internal procedures, implemented by gfortran
with trampolines. Do not link with `-z noexecstack`.

## Data

- `nist/`: the NIST StRD nonlinear regression data sets (public domain), as
  downloaded from https://www.itl.nist.gov/div898/strd/nls/nls_main.shtml.
- `make_instances.py`: contaminates them with 5% and 10% of outliers (fixed seed),
  checks the certified residual sums of squares, and writes `instances/` together
  with the polynomial and Bard problems.

## Experiments and results

| Section / table | Script | Results |
|---|---|---|
| 4.1-4.4, Tables 1-4, Figure 1 | `run_paper_tables.py`, `make_paper_tables.py`, `plot_bard_panels.py` | `results/paper/`, `results/paper_tables.tex` |
| 4.5, Tables 5-6, Figure 2 | `run_benchmark.py`, `analyze_benchmark.py`, `make_bench_tables.py` | `results/raw/`, `results/benchmark.csv`, `results/bench_tables.tex`, `results/summary_main.txt` |
| 4.5, Table 7 (size of the active set) | `run_active.py`, `analyze_active.py` | `results/active/`, `results/active_table.tex` |
| 4.6, Table 8 (bound M) | `run_benchmark.py --variants ...`, `run_adaptive.py`, `make_bench_tables.py` | `results/raw/`, `results/bench_tables.tex` |
| 4.6, shift margin tau | `run_tau.sh` | `results/tau/` |

`merge_raw.py` rebuilds `results/benchmark.csv` from `results/raw/`. Each line of a
raw file is one run: instance, number of discarded observations, curvature
family, M, norm control, tau, delta, sigma_min, starting point, trial, status
(0: stopped at Step 5, 1: iteration limit, 2: CPU-time limit, 3: subproblem solver
did not return), number of trials accepted by the inner-loop safeguard,
iterations, CPU time (s), best value, function evaluations, and number of true
outliers among the discarded observations. Runs were pinned to the performance
cores of an Intel Core i7-13700KF.
