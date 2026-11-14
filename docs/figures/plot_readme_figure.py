#!/usr/bin/env python3
# Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE)
"""Figure of the README: the least damped oscillatory mode of Schwarzschild-de Sitter black holes
against the dimensionless cosmological constant, read from the reports in data/nonextreme/.

From this folder,
    python3 plot_readme_figure.py
writes massless-sds-modes.png. Panels (a) and (b): real part and damping rate of the least damped
mode with Re Omega > 0. Panel (c): the ratio Re Omega / (-Im Omega), with the values 2 sqrt(A - 1/4)
of the Nariai limit (A = l(l+1) for s = 0, 1 and l(l+1) - 2 for s = 2). A folder whose report is
named after another case is skipped. The style is the stamped copy manuscript_style.py beside this script.
"""
import math
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
DATA = HERE.parent.parent / "data" / "nonextreme"
sys.path.insert(0, str(HERE))
import manuscript_style as hs  # noqa: E402

SERIES = (  # folder pattern, report prefix, legend, A of the Nariai limit
    ("s0_L0_La*", "s0_L0", r"scalar, $\ell = 0$", 0),
    ("s0_L1_La*", "s0_L1", r"scalar, $\ell = 1$", 2),
    ("s0_L2_La*", "s0_L2", r"scalar, $\ell = 2$", 6),
    ("s1_L1_La*", "s1_L1", r"electromagnetic, $\ell = 1$", 2),
    ("s2/odd/s2_L2_La*.txt", "s2_L2", r"axial gravitational, $\ell = 2$", 4),
)
NARIAI = 1.0 / 9.0


def lambda_of(name):
    """Dimensionless cosmological constant from a name such as s0_L1_La0_06 or s2_L2_La0_10."""
    digits = re.search(r"La0_(\d+)", name).group(1)
    return float("0." + digits) if len(digits) > 1 else float("0." + digits + "0")


def roots(report):
    """Roots of a report, in the format of QNM Analyser or in that of the Julia solver."""
    text = report.read_text()
    found = []
    if "QNM Convergence Report" in text:
        inside = False
        for line in text.splitlines():
            if re.match(r"^(General|Purely Real|Purely Imaginary) QNMs", line):
                inside = True
                continue
            parts = line.split()
            if inside and len(parts) in (2, 3):
                try:
                    found.append(complex(float(parts[0]), float(parts[1])))
                except ValueError:
                    pass
    else:
        for line in text.splitlines():
            parts = line.split()
            if len(parts) == 7 and parts[0].isdigit():
                found.append(complex(float(parts[1]), float(parts[2])))
    return np.array(found)


def least_damped(pattern, prefix):
    points = []
    for path in sorted(DATA.glob(pattern)):
        reports = [path] if path.is_file() else sorted(path.glob(prefix + "_La*.txt"))
        if not reports:
            continue                      # the folder holds the report of another case
        z = roots(reports[0])
        oscillatory = z[(z.real > 1e-6) & (z.imag < 0)]
        if oscillatory.size:
            points.append((lambda_of(path.name), oscillatory[np.argmax(oscillatory.imag)]))
    points.sort(key=lambda p: p[0])
    return np.array([p[0] for p in points]), np.array([p[1] for p in points])


hs.configure(family="libertine")
plt.rcParams["savefig.dpi"] = 265
plt.rcParams["figure.facecolor"] = hs.PAPER
plt.rcParams["savefig.facecolor"] = hs.PAPER

fig, (ax_re, ax_im, ax_ratio) = plt.subplots(1, 3, figsize=(hs.WIDTH, 2.45))
for k, (pattern, prefix, label, area) in enumerate(SERIES):
    lam, omega = least_damped(pattern, prefix)
    style = hs.series(k)
    style["markersize"] = 2.8
    ax_re.plot(lam, omega.real, lw=0.9, label=label, **style)
    ax_im.plot(lam, -omega.imag, lw=0.9, **style)
    ax_ratio.plot(lam, omega.real / -omega.imag, lw=0.9, **style)
    limit = 2.0 * math.sqrt(area - 0.25) if area > 0 else 0.0
    ax_ratio.axhline(limit, color=hs.SLATE, lw=0.6, linestyle=(0, (1.6, 1.6)), zorder=0)

for axis in (ax_re, ax_im, ax_ratio):
    axis.set_xlim(0.0, 0.1155)
    axis.set_xticks([0.0, 0.05, 0.1])
    axis.axvline(NARIAI, color=hs.MIST, lw=0.7, zorder=0)
    axis.set_xlabel(r"$\tilde\Lambda = M^2 \Lambda$")
ax_re.set_ylim(0.0, 0.45)
ax_re.set_ylabel(r"$\mathrm{Re}\,\Omega$")
hs.panel(ax_re, "(a) frequency")
ax_im.set_ylim(0.0, 0.11)
ax_im.set_ylabel(r"$-\mathrm{Im}\,\Omega$")
hs.panel(ax_im, "(b) damping rate")
ax_ratio.set_ylim(0.0, 5.2)
ax_ratio.set_ylabel(r"$\mathrm{Re}\,\Omega \,/\, (-\mathrm{Im}\,\Omega)$")
hs.panel(ax_ratio, "(c) ratio")

hs.legend_outside(fig, where="upper", ncol=3)
out = HERE / "massless-sds-modes.png"
hs.save(fig, out)
print("written", out)
