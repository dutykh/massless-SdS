<!-- Author: Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE) -->

# massless-SdS: quasinormal modes of Schwarzschild–de Sitter black holes and the Nariai limit

![The least damped quasinormal mode of five perturbations of Schwarzschild-de Sitter black holes against the cosmological constant, and the ratio of its real part to its damping rate](docs/figures/massless-sds-modes.png)

The least damped oscillatory quasinormal mode Ω = Mω of five perturbations against the dimensionless cosmological
constant Λ̃ = M²Λ, read from the tables of `data/nonextreme/` (the folder of the scalar monopole at Λ̃ = 0.11 holds no
table of its own, see Contents). Panels (a) and (b) show that the frequency and the damping rate both decrease towards
zero as Λ̃ approaches 1/9 (vertical line), where the event and cosmological horizons merge. Panel (c) shows that their
ratio tends to 2√(A_s − 1/4) (dotted lines), the ratio of the fundamental frequency of the Pöschl–Teller problem of the
Nariai limit, with A_s = ℓ(ℓ+1) for the scalar and electromagnetic fields and ℓ(ℓ+1) − 2 for the gravitational field;
at Λ̃ = 0.11 the four ratios are within 0.2 % of their limits. For the scalar monopole, A_s = 0 and the limit is 0.

massless-SdS holds the codes and data of a study of the quasinormal modes of four-dimensional Schwarzschild–de Sitter
black holes by Davide Batic, Anna Chrysostomou, Alan S. Cornell and Denys Dutykh, for massless scalar, electromagnetic
and gravitational perturbations (the latter in both parities). After the behaviour at the two horizons is factored out,
the radial problem becomes a quadratic eigenvalue problem on −1 ⩽ y ⩽ 1, solved by a Chebyshev spectral method for
0 < Λ̃ < 1/9; the Nariai limit Λ̃ = 1/9 is treated as a separate spectral problem, whose exact Pöschl–Teller spectrum
serves as a benchmark. Maple assembles the matrices, MATLAB (with the Advanpix toolbox) or Julia solves the problem,
and Mathematica notebooks compare the results with WKB estimates and draw the spectra. The manuscript, "Quasinormal
modes of Schwarzschild–de Sitter black holes: spectral branches and the Nariai limit", is included as
`pdf/SCHWARZSCHILD_DE_SITTER_MASSLESS.pdf`; it has not been published, so the software itself is to be cited (see How
to cite).

## The problem

In units with c = G = ℏ = 1 and with the time dependence e^(−iωt), the line element and the radial equation of the
scalar (s = 0) and electromagnetic (s = 1) perturbations are

```text
ds² = −f(r) dt² + dr²/f(r) + r² (dϑ² + sin²ϑ dφ²),      f(r) = 1 − 2M/r − Λr²/3,      Λ > 0,

f d/dr( f dψ/dr ) + ( ω² − U(r) ) ψ = 0,      U(r) = f(r) ( ε f′(r)/r + ℓ(ℓ+1)/r² ),      ε = 1 − s².
```

Gravitational perturbations of odd parity replace U by the Regge–Wheeler potential, and those of even parity by the
Zerilli-type potential of equation (6) of the manuscript:

```text
V_odd(r) = f(r) ( ℓ(ℓ+1)/r² − 6M/r³ ),      ℓ ⩾ 2.
```

With x = r/(2M), Λ̃ = M²Λ and Ω = Mω, the lapse is f(x) = 1 − 1/x − (4Λ̃/3) x². For 0 < Λ̃ < 1/9 it has two positive
zeros, the event horizon x_h and the cosmological horizon x_c,

```text
x_h = cos( (π + α)/3 ) / √Λ̃,      x_c = cos( (π − α)/3 ) / √Λ̃,      α = arccos( 3√Λ̃ ),
```

which merge at x = 3/2 when Λ̃ = 1/9. A quasinormal mode is ingoing at the event horizon and outgoing at the
cosmological horizon. Once these behaviours are factored out, the remaining function is expanded in Chebyshev
polynomials of the first kind on −1 ⩽ y ⩽ 1 and the equation is collocated, which gives the quadratic pencil

```text
( M0 + i Ω M1 + Ω² M2 ) a = 0.
```

It is solved at three resolutions, and a root is kept when the three agree: the nearest roots at the three resolutions
must lie within max(10⁻⁴, 10⁻⁴ · max(1, |Ω|)) of each other, and their average is reported.

In the Nariai limit the radial equation becomes a Pöschl–Teller problem in its own time coordinate, with the frequency
ϖ and the coefficient A_s = ℓ(ℓ+1) for s = 0 and 1 and A_s = ℓ(ℓ+1) − 2 for s = 2. Its quasinormal frequencies are
known exactly (equation (90) of the manuscript), and near extremality they are related to Ω by

```text
ϖ_n = ± √(A_s − 1/4) − i (n + 1/2),   n = 0, 1, 2, ...,            ϖ = Ω / ( Λ̃ (x_c − x_h) ),
```

so that the raw frequencies Ω are not to be compared with ϖ directly.

## Quick start

First, Maple assembles the three matrices `M0`, `M1` and `M2` of one case at each resolution: the worksheets
`maple/Matrix_Assembler_*.mw` (see Contents) set the sector, the multipole ℓ, the cosmological constant and the
resolutions, and write the matrices as ASCII files into the folder that the Julia solver reads, under the names it
expects (as saved, all four worksheets last assembled N = 90, 95 and 100).

Second, the pencil is solved, in one of two ways. In Julia, `julia/QNMSdS_readable_output.jl` reads the matrices
`M0-<N>-<tag>.mat`, `M1-<N>-<tag>.mat` and `M2-<N>-<tag>.mat`, the tag being 1000 Λ̃, linearises the pencil, solves it in
`BigFloat` arithmetic with GenericSchur, applies the three-resolution test and writes the table of retained modes
(`stable_modes-<tag>.txt`, with `_details` and `_full_precision` companions), the spectra and pictures of them;
`julia/QNM_Nariai.jl` does the same for the Nariai problem and also keeps, among the retained roots, only those within
10⁻⁸ of an exact Pöschl–Teller frequency. Both Julia files set their parameters as constants near the top, to be
edited before a run: the folder of the matrices `DATA_DIR` (a folder of the authors' computer by default), the
resolutions `N_VALS` (90, 95 and 100 by default), the tag or the sector (`FILETAG`, `SPIN`, `ELL`) and the precision
(70 digits, the matrices coming from Maple with 50; the constant `FILE_DECIMAL_DIGITS` is not used). In MATLAB with the Advanpix Multiprecision Computing Toolbox,
`matlab/chaseeigs.m` reads the resolutions from `resolutions.txt` and s, ℓ and Λ̃ from the first three lines of
`params.txt` in the working folder, and the matrices from `assemble/M0_<N>.mat`, `assemble/M1_<N>.mat` and
`assemble/M2_<N>.mat` (other names than those of the Julia solver, so the files are to be copied under these names);
it solves the problem with `mp.Digits(N)`, as many decimal digits as the resolution, and writes every root, real and
imaginary parts in two columns, into `results/eigs_<N>.dat`. Its first lines add the toolbox from a fixed folder, which
must be changed to that of the local installation.

Third, the notebooks of `mathematica/` compare and draw the results. `CompareViaZhidenko_SdS.nb` loads the WKB package
from its own folder, so `matlab/WKB.m` must first be copied into `mathematica/`.

## Contents

| Folder | What it holds |
|---|---|
| `maple/` | the matrix assemblers `Matrix_Assembler_s0_s1-Final Version.mw` (scalar and electromagnetic), `Matrix_Assembler_s2odd-Final Version.mw` (odd gravitational), `Matrix_Assembler_s2even-Final Version.mw` (even gravitational) and `Matrix_Assembler_all_s_Nariai-Final Version.mw` (Nariai limit), and the five worksheets `01-...-SUPPLEMENTAL MATERIAL.mw` with the derivations (Frobenius analysis, compactification, coefficients) |
| `matlab/` | `chaseeigs.m`, the multiple-precision solver, and `WKB.m`, a Wolfram Language package (not MATLAB code) of higher-order WKB formulas used by the notebook, credited below |
| `julia/` | `QNMSdS_readable_output.jl`, the solver and the three-resolution test for 0 < Λ̃ < 1/9, and `QNM_Nariai.jl`, its counterpart for the Nariai limit |
| `mathematica/` | `CompareViaZhidenko_SdS.nb`, the potentials and the comparison of the spectral modes with WKB and semi-analytic estimates, and `QNMplots.nb`, which reads tables of modes and draws spectra in the complex plane (for instance the gravitational ℓ = 2 case at Λ̃ = 0.11 and its Nariai counterpart) |
| `data/nonextreme/` | one folder per scalar or electromagnetic case, with a table of modes and the raw spectra in `raw/`; the gravitational tables in `s2/odd/` and `s2/even/` |
| `data/extreme/` | the four tables of the Nariai limit |
| `pdf/` | the manuscript (6 June 2026) |
| `docs/figures/` | the script that draws the figure of this page from the tables of `data/nonextreme/`, and the style file beside it |

The folders of `data/nonextreme/` are named `s<s>_L<ℓ>_La0_<digits>`, the digits being those of Λ̃ after the decimal
point (`La0_06` is Λ̃ = 0.06, `La0_1` is 0.10 and `La0_11` is 0.11). They cover the scalar field with ℓ = 0, 1 and 2 and
the electromagnetic field with ℓ = 1 and 2, at Λ̃ = 0.02, 0.04, 0.06, 0.08, 0.09, 0.10 and 0.11; the odd gravitational
tables cover ℓ = 2 and 3 at the same seven values, and the even ones ℓ = 2 at Λ̃ = 0.02, 0.06, 0.08 and 0.09, for the
comparison of the two parities. Their common modes agree, as the Chandrasekhar–Darboux relation predicts; the odd
tables list a few more overtones and, at Λ̃ = 0.06, 0.08 and 0.09, the purely imaginary mode Ω = −2i, which the
manuscript identifies as an algebraically special mode of the axial sector, absent from the even one. The
tables list the modes with Re Ω ⩾ 0, the spectra being symmetric about the imaginary axis. Six of them, the scalar
ℓ = 0 and ℓ = 1 at Λ̃ = 0.02, 0.04 and 0.06, are reports of QNM Analyser made from the raw spectra beside them at
N = 190, 195 and 200; the others were written by `QNMSdS_readable_output.jl` and give, for each mode, the running index
`k`, the average `Re(avg)` and `Im(avg)` and the three distances `err12`, `err23`, `err13` between the resolutions, with
their maximum `maxerr`. The raw spectra, `raw/eigs_<N>.dat` for N = 190, 195, 200, 290, 300 and 310, are in the format
that `chaseeigs.m` writes and hold every root (2N of them) with as many digits as the resolution.

Three points of these data call for care. At Λ̃ = 0.02, 0.04 and 0.06 the modes of the Julia tables are found in the
raw spectra beside them, but the distances `err12` and `err23` of the tables are not those of these spectra, which
agree more closely: the tables come from another run of the solver, whose spectra are not stored. At Λ̃ = 0.08, 0.09,
0.10 and 0.11 the raw spectra stored for the five cases do not reproduce the tables: at most two of their several
hundred roots pass the three-resolution test, all on the imaginary axis and some with Im Ω > 0, and none of the
tabulated modes is among their roots, so these tables too come from computations whose spectra are not stored.
Finally, the folder `s0_L0_La0_11` holds a copy of the table of ℓ = 1 (`s0_L1_La0_11.txt`, identical to the one in
`s0_L1_La0_11`), so that the table of the scalar monopole at Λ̃ = 0.11 is missing.

An independent computation confirms the tables themselves. A collocation of the same equations written afresh (Gauss
nodes, 90-digit arithmetic, N = 90, 95 and 100) finds every mode of the 19 tables at Λ̃ = 0.08 to 0.11 within its
tabulated error (the largest offset is 3 × 10⁻⁵), so that only the stored spectra are faulty there, and it reproduces
the three tables at Λ̃ = 0.02, their distances included, which therefore come from a run at N = 90, 95 and 100. For the
scalar monopole at Λ̃ = 0.11 it retains 0.0020352 − 0.0189819i, 0.0038667 − 0.0380418i, 0.0054057 − 0.0571828i,
0.0066629 − 0.0763740i and 0.0076830 − 0.0955908i; the missing table is rebuilt by running the assembler and the solver
for ℓ = 0 and Λ̃ = 0.11.

The tables of `data/extreme/` give the Nariai frequencies ϖ: `s_0_L_0.txt` (scalar, ℓ = 0, A_s = 0), `s_0_s_1_L_1.txt`
(scalar or electromagnetic, ℓ = 1, A_s = 2), `s_2_L_2.txt` (gravitational, ℓ = 2, A_s = 4) and `S_2_L_3.txt`
(gravitational, ℓ = 3, A_s = 10), in the format of the Julia tables. Every tabulated frequency lies within 6 × 10⁻⁹ of
an exact value, for instance ±1.3228756555322953 − 0.5i = ±√(7/4) − i/2 for A_s = 2; the filter of `QNM_Nariai.jl`,
which keeps only the roots within 10⁻⁸ of an exact value, makes this agreement partly a matter of construction.

## Requirements

Maple, for the assemblers and the derivations. MATLAB with the Advanpix Multiprecision Computing Toolbox (a commercial
toolbox), for `chaseeigs.m`; or Julia with GenericSchur and Plots (and PlotlyJS for the interactive pictures),
DelimitedFiles, Printf and LinearAlgebra coming with Julia. Wolfram Mathematica, for the notebooks. The figure of this
page is drawn from `docs/figures/` by `python3 plot_readme_figure.py` (Python 3 with NumPy and Matplotlib).

## Third-party code

`matlab/WKB.m` is the WKB package distributed with the article of R. A. Konoplya, A. Zhidenko and A. F. Zinhailo,
Higher order WKB formula for quasinormal modes and grey-body factors: recipes for quick and accurate calculations,
Classical and Quantum Gravity 36, 155002 (2019),

https://doi.org/10.1088/1361-6382/ab2e25

whose header credits the WKB formulas of S. Iyer and C. M. Will (third order), R. A. Konoplya (sixth order) and
J. Matyjasek and M. Opala (thirteenth order). It belongs to its authors, and the licence of this repository does not
cover it.

## How to cite

The manuscript has not been published. If you use these codes or data, please cite the software:

D. Batic, A. Chrysostomou, A. S. Cornell and D. Dutykh, massless-SdS: quasinormal modes of Schwarzschild–de Sitter
black holes and the Nariai limit, software, version 1.0.0 (2026).

https://github.com/dutykh/massless-SdS

The button "Cite this repository" in the right-hand column of this page gives the same reference, read from
`CITATION.cff`.

The quasinormal-mode codes of this family are gathered in the QNMs Hall of Fame:

https://qnms.denys-dutykh.com/

## Licence

GNU Lesser General Public License, version 2.1 (see `LICENSE`), except `matlab/WKB.m` (see Third-party code).

## Authors

Davide Batic (Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi), Anna Chrysostomou
(Laboratoire de Physique Théorique et Hautes Énergies, Sorbonne Université, CNRS, Paris), Alan S. Cornell (Department
of Physics, University of Johannesburg) and Denys Dutykh (Mathematics Department, Khalifa University of Science and
Technology, Abu Dhabi).

https://www.denys-dutykh.com/
