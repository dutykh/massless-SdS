# massless-SdS

**Quasinormal modes of Schwarzschild–de Sitter black holes: spectral branches and the Nariai limit**

This repository contains the computational codes and numerical data accompanying
our manuscript on the quasinormal mode (QNM) spectrum of four-dimensional
Schwarzschild–de Sitter (SdS) black holes. We compute the QNMs of massless
**scalar** ($s=0$), **electromagnetic** ($s=1$), and **gravitational** ($s=2$,
in both odd/axial and even/polar parities) perturbations using a **Chebyshev
spectral method** adapted to the two-horizon structure of the spacetime. After
factoring out the Frobenius asymptotics — ingoing waves at the black hole
horizon and outgoing waves at the cosmological horizon — the radial problem is
recast as a quadratic spectral eigenvalue problem on the finite interval
$[-1,1]$. We track both the oscillatory branches and the purely imaginary
(overdamped) modes as functions of the dimensionless cosmological constant, and
we treat the extremal **Nariai limit** as a separate spectral problem whose
analytic Pöschl–Teller spectrum benchmarks the method.

The manuscript is included here as
[`pdf/SCHWARZSCHILD_DE_SITTER_MASSLESS.pdf`](pdf/SCHWARZSCHILD_DE_SITTER_MASSLESS.pdf).

## Authors

- **Davide Batic** — Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE
- **Anna Chrysostomou** — Laboratoire de Physique Théorique et Hautes Énergies (LPTHE), Sorbonne Université, CNRS, Paris, France
- **Alan S. Cornell** — Department of Physics, University of Johannesburg, Auckland Park, South Africa
- **Denys Dutykh** — Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE

## Physical and numerical setup

The SdS line element uses the lapse $f(r) = 1 - 2M/r - \Lambda r^2/3$ with a
positive cosmological constant $\Lambda$. After spherical-harmonic and Fourier
decomposition, each perturbation reduces to a radial wave equation on the
interhorizon region $r_h < r < r_c$ (event horizon to cosmological horizon),
governed by a spin- and multipole-dependent effective potential:

- **Scalar / electromagnetic** ($s=0,1$): a single master equation with
  $\mathcal{A}_s = \ell(\ell+1)$.
- **Gravitational** ($s=2$): the **odd/axial** sector uses the Regge–Wheeler
  potential and the **even/polar** sector uses the Zerilli-type potential; both
  reduce to $\mathcal{A}_2 = \ell(\ell+1)-2$ in the Nariai limit and are
  expected to be isospectral (Chandrasekhar–Darboux), except for an isolated
  algebraically special axial mode.

The control parameter is the **dimensionless cosmological constant**
$\widetilde{\Lambda} = M^2 \Lambda$, with the de Sitter range $0 < \widetilde{\Lambda} < 1/9$;
the critical value $\widetilde{\Lambda} = 1/9$ is the extremal (Nariai) point where
the two horizons coalesce. The data here span
$\widetilde{\Lambda} \in \{0.02, 0.04, 0.06, 0.08, 0.09, 0.10, 0.11\}$.

The QNM boundary conditions are factored out analytically; the regular part is
expanded in Chebyshev polynomials of the first kind, and the collocation
procedure produces a **quadratic matrix pencil**

$$\big(M_0 + i\,\lambda\, M_1 + \lambda^2 M_2\big)\,\mathbf{a} = 0,$$

where $\lambda = \Omega = M\omega$ is the (rescaled) non-extremal SdS frequency,
and $\lambda = \varpi$ is the Nariai frequency conjugate to the blown-up time
coordinate. With the time convention $e^{-i\omega t}$, damped modes have
$\Im(\omega) < 0$. **Note:** the raw SdS frequencies $\Omega$ should not be
compared directly with the Nariai frequencies $\varpi$; near extremality they
are related by $\varpi = \Omega / [\widetilde{\Lambda}(x_c - x_h)]$, as discussed
in the manuscript.

## Repository layout

```
massless-SdS/
├── pdf/                    Manuscript (PDF)
├── maple/                  Symbolic derivations + spectral matrix assemblers
├── matlab/                 Multiprecision eigenvalue solver + WKB package
├── julia/                  Eigenvalue solvers / post-processing (SdS & Nariai)
├── mathematica/            Benchmarks (WKB, Pöschl–Teller, Zhidenko) and plots
├── data/
│   ├── extreme/            Filtered Nariai QNM tables
│   └── nonextreme/         Non-extremal SdS spectra (reports + raw eigenvalues)
├── assets/                 (reserved)
├── LICENSE                 GNU LGPL v2.1
└── README.md
```

### `maple/`

Symbolic worksheets (Maple `.mw`). The **matrix assemblers** build the
Chebyshev collocation matrices $M_0, M_1, M_2$ for the quadratic pencil and
export them as ASCII matrices (assembled at 200 decimal digits):

- `Matrix_Assembler_s0_s1-Final Version.mw` — scalar and electromagnetic sectors
- `Matrix_Assembler_s2even-Final Version.mw` — gravitational even/polar (Zerilli)
- `Matrix_Assembler_s2odd-Final Version.mw` — gravitational odd/axial (Regge–Wheeler)
- `Matrix_Assembler_all_s_Nariai-Final Version.mw` — all sectors in the Nariai limit

The **supplemental-material worksheets** (`01-...SUPPLEMENTAL MATERIAL.mw`)
contain the full derivations behind the manuscript: the Frobenius/indicial
analysis, the compactification to $[-1,1]$, and the explicit coefficient
functions, for the Nariai geometry and for the non-extremal scalar/EM, even, and
odd gravitational cases.

### `matlab/`

- `chaseeigs.m` — solves the quadratic pencil $M_0 + i\Omega M_1 + \Omega^2 M_2 = 0$
  with the [Advanpix Multiprecision Computing Toolbox](https://www.advanpix.com/)
  `polyeig`. It reads the sector parameters $(s, L, \Lambda)$ from `params.txt`
  and the list of resolutions from `resolutions.txt`, then writes the full
  multiprecision eigenvalues to `results/eigs_<N>.dat`.
- `WKB.m` — a large WKB-formula package (after Konoplya–Zhidenko, arXiv:1904.10333)
  used by the Mathematica benchmark notebook for the higher-order WKB estimates.

### `julia/`

- `QNMSdS_readable_output.jl` — non-extremal post-processor. Reads the
  Maple-exported matrices at three resolutions, linearises the quadratic pencil
  into a generalized eigenvalue problem (solved with `GenericSchur` in extended
  `BigFloat` precision), retains modes that are stable across all three
  resolutions, and writes human-readable convergence reports (general,
  purely-imaginary, and purely-real modes) together with a full-precision dump.
- `QNM_Nariai.jl` — the Nariai (extremal) counterpart. In addition to the
  resolution-stability filtering, it compares the numerical spectrum against the
  exact analytic Pöschl–Teller frequencies and writes both the raw and the
  analytically filtered stable modes.

> The two Julia scripts define a small number of user-editable constants near
> the top — the input/output directory (`DATA_DIR`), the field sector
> (`SPIN`, `ELL`) and file tag, and the list of resolutions (`N_VALS`). Set
> these to match your Maple output before running.

### `mathematica/`

- `CompareViaZhidenko_SdS.nb` — defines the scalar, electromagnetic, and
  gravitational effective potentials and cross-checks the spectral QNMs against
  WKB, Pöschl–Teller, and Zhidenko-type semi-analytic estimates (imports
  `matlab/WKB.m`).
- `QNMplots.nb` — reads the data tables in `data/` and produces the
  complex-plane spectrum figures used in the manuscript (exported to PDF).

## Data

All frequencies use the convention $e^{-i\omega t}$, so damped modes have
$\Im < 0$. Because every spectrum is symmetric about the imaginary axis, the
reports list only $\Re \ge 0$.

### Naming convention

Non-extremal cases are named `s{spin}_L{ell}_La0_{1000·Λ}`, e.g.
`s0_L2_La0_06` is the scalar ($s=0$), $\ell=2$, $\widetilde{\Lambda}=0.06$ case;
`s1_L1_La0_11` is electromagnetic, $\ell=1$, $\widetilde{\Lambda}=0.11$.
(`La0_1` denotes $\widetilde{\Lambda}=0.10$ and `La0_11` denotes $0.11$.)

### `data/extreme/` — Nariai (extremal) spectra

Four whitespace-aligned tables of the **filtered** physical Nariai QNMs, one per
sector:

| File | Sector |
|------|--------|
| `s_0_L_0.txt`     | scalar, $\ell=0$ |
| `s_0_s_1_L_1.txt` | scalar / electromagnetic, $\ell=1$ |
| `s_2_L_2.txt`     | gravitational, $\ell=2$ |
| `S_2_L_3.txt`     | gravitational, $\ell=3$ |

Each file carries a header documenting the sector, the value of $\mathcal{A}_s$,
and the pencil, followed by columns:

```
k    Re(avg)    Im(avg)    err12    err23    err13    maxerr
```

where `k` indexes the mode, `Re(avg)`/`Im(avg)` are the resolution-averaged
frequency $\varpi$ (16 decimal places), and `err12, err23, err13, maxerr` are
the pairwise/maximum differences across the three resolutions (scientific
notation). These reproduce the analytic Pöschl–Teller values to the displayed
precision.

### `data/nonextreme/` — non-extremal SdS spectra

One directory per case. Scalar and electromagnetic cases each contain:

- **`<case>.txt`** — a *QNM Convergence Report*: the date, tolerance, resolution
  set, and total counts, followed by the converged **general** ($\Re\neq0$,
  $\Im\neq0$), **purely imaginary**, and **purely real** modes.
- **`raw/`** — the unfiltered solver output, `eigs_<N>.dat`, with **two
  whitespace-separated columns** `Re(ω) Im(ω)` at full multiprecision. Each file
  contains all $2N$ candidate eigenvalues of the linearised pencil at resolution
  $N$. The primary resolutions are $N \in \{190, 195, 200\}$; the
  $\{290, 300, 310\}$ files are higher-resolution cross-checks.

Gravitational ($s=2$) results live under `data/nonextreme/s2/`:

- **`odd/`** — Regge–Wheeler / axial reports for $\ell=2$ and $\ell=3$ across the
  full $\widetilde{\Lambda}$ range.
- **`even/`** — Zerilli / polar reports ($\ell=2$). These are provided as the
  independent **isospectrality check** against the odd sector (the isolated
  algebraically special axial mode is, as expected, absent from the even
  spectrum).

## Computational pipeline (overview)

1. **Assemble** the Chebyshev collocation matrices $M_0, M_1, M_2$ symbolically
   in **Maple** (`maple/Matrix_Assembler_*.mw`), exporting them at high
   precision.
2. **Solve** the resulting quadratic eigenvalue problem at several spectral
   resolutions, either in **MATLAB** (`matlab/chaseeigs.m`, Advanpix `polyeig`)
   or in **Julia** (`julia/*.jl`, linearisation + `GenericSchur`). This yields
   all candidate eigenvalues — the `data/.../raw/*.dat` files.
3. **Filter** by retaining only the frequencies that match across resolutions
   (relative/absolute tolerance $10^{-4}$); the reported frequency is the
   average over the matched triplet, with the pairwise spreads kept as
   convergence diagnostics. In the Nariai case, modes are additionally compared
   against the exact Pöschl–Teller spectrum. This produces the readable `.txt`
   reports.
4. **Benchmark and visualise** with **Mathematica**: WKB / Pöschl–Teller /
   Zhidenko cross-checks (`CompareViaZhidenko_SdS.nb`) and the complex-plane
   spectrum plots (`QNMplots.nb`).

## Software requirements

- **Maple** — symbolic derivations and matrix assembly.
- **MATLAB** with the **Advanpix Multiprecision Computing Toolbox** — for
  `chaseeigs.m`.
- **Julia** with `DelimitedFiles`, `Printf`, `LinearAlgebra`, `GenericSchur`,
  and `Plots`.
- **Wolfram Mathematica** — for the benchmarks and plots.

## How to cite

If you use these codes or data, please cite:

> D. Batic, A. Chrysostomou, A. S. Cornell, and D. Dutykh,
> *Quasinormal modes of Schwarzschild de Sitter black holes: spectral branches
> and the Nariai limit* (2026).

## License

This project is released under the **GNU Lesser General Public License v2.1**.
See [`LICENSE`](LICENSE) for the full text.
