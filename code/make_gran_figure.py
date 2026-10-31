# -*- coding: utf-8 -*-
"""Figure 2: granularity of the source and protection of the past in voids.   python make_gran_figure.py   (needs granularity.json)"""
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.abspath(os.path.join(HERE, "..", "figures"))
G = json.load(open(os.path.join(HERE, "granularity.json"), encoding="utf-8"))
plt.rcParams.update({"font.size": 8.3, "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.6, "mathtext.default": "regular"})
C1, C2, C3, C4, GR = "#3E6B8A", "#B5544E", "#5C8A5C", "#8A7A3E", "#777777"

fig, ax = plt.subplots(1, 2, figsize=(7.6, 3.1))
a = ax[0]
t = G["trigger"]
lo, hi = -4.0, 15.0
bars = [("point sources averaged:  n >> $m_\\varphi^3$", lo, np.log10(t["particle_max_mass_meV"]), C2),
        ("particle-like:  occupation << 1", np.log10(t["classical_max_mass_meV"]), hi, C2),
        ("wave averaged in time:  $m_{dm}$ >> $m_\\varphi$", np.log10(t["m_phi_meV"]), hi, C1),
        ("wave classical:  occupation >> 1", lo, np.log10(t["classical_max_mass_meV"]), C1)]
for i, (lab, x0, x1, c) in enumerate(bars):
    a.barh(i, x1 - x0, left=x0, height=0.55, color=c, alpha=0.55)
    a.text(lo + 0.3, i + 0.36, lab, fontsize=6.4, va="bottom")
w0, w1 = np.log10(G["wave_window_meV"][0]), np.log10(G["wave_window_meV"][1])
a.axvspan(w0, w1, color=C3, alpha=0.25)
a.text(0.5 * (w0 + w1), 3.95, "mean field\nvalid", ha="center", fontsize=6.4, color=C3)
for m, lab in ((9.0, "1 MeV"), (12.0, "1 GeV")):
    a.axvline(m, color=GR, lw=0.6, ls=":")
    a.text(m + 0.15, -0.55, lab, fontsize=6.0, color=GR)
a.set_yticks([]); a.set_xlim(lo, hi); a.set_ylim(-0.7, 4.6)
a.set_xlabel(r"log$_{10}$($m_{dm}$ / meV)")
a.set_title("(a) When the field feels a smooth density", fontsize=8.5, fontweight="bold")
a = ax[1]
g = np.array(G["void_curves"]["one_plus_delta"])
for (k, v), c in zip(G["void_curves"]["lnI"].items(), (C2, C4, C1)):
    a.plot(g, v, color=c, label=r"mean $\theta_0/\theta_{deg}$ = %s" % k)
a.axhline(0.0, color="k", lw=0.6)
a.text(0.06, 12, "decayed in the past", fontsize=6.5); a.text(0.06, -40, "past light cone clear", fontsize=6.5)
a.set_ylim(-150, 150); a.set_xlim(0.05, 1.0)
a.set_xlabel(r"local density / mean, $1+\delta$"); a.set_ylabel("ln(expected bubbles in the past light cone)")
a.legend(fontsize=6.3, loc="upper right")
a.set_title("(b) Protection of the past in voids", fontsize=8.5, fontweight="bold")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "fig2_granularity.png"), dpi=300, bbox_inches="tight", facecolor="white")
print("fig2_granularity.png")
