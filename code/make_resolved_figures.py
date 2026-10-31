# -*- coding: utf-8 -*-
"""Figures for the spatially resolved tests: fig3_resolved_granularity.png (granular_field.json) and fig6_wave_wall.png (wave_wall_1d.json).
    python make_resolved_figures.py"""
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.abspath(os.path.join(HERE, "..", "figures"))
J = lambda f: json.load(open(os.path.join(HERE, f), encoding="utf-8"))
plt.rcParams.update({"font.size": 8.3, "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.6, "mathtext.default": "regular"})
C = ["#3E6B8A", "#B5544E", "#5C8A5C", "#8A7A3E", "#777777", "#7A5C8A", "#2E2E2E"]

# ---------------------------------------------------------------- granularity
GF, GR = J("granular_field.json"), J("granularity.json")
fig, ax = plt.subplots(1, 2, figsize=(7.6, 3.1))
a = ax[0]
for c, col in zip([c for c in GF["cells"] if c["Rc"] in (0.3, 3.0, 10.0, 30.0, 100.0)], C):
    rr = np.array(c["profile_r"]); a.plot(rr / c["Rc"], c["profile_x"], color=col, label=r"$R_c m_\varphi$ = %.1f" % c["Rc_times_mphi"])
c0 = GF["cells"][0]
a.axhline(c0["x_h_zero_density"], color="k", ls=":", lw=0.8); a.axhline(c0["x_h_mean_field"], color="k", ls="--", lw=0.8)
a.text(0.55, c0["x_h_zero_density"] - 0.006, "zero-density vacuum", fontsize=6.3); a.text(0.55, c0["x_h_mean_field"] + 0.002, "mean-field vacuum", fontsize=6.3)
a.set_ylim(0.775, 0.84); a.set_xlabel(r"r / $R_c$ (cell radius)"); a.set_ylabel(r"field x = $\varphi$/f")
a.legend(fontsize=6.0, loc="upper right")
a.set_title("(a) Field in a cell around one source", fontsize=8.5, fontweight="bold")
a = ax[1]
rows = GF["rate_vs_mass"]
m = np.array([r["m_dm_meV"] for r in rows]); fr = np.array([r["fraction_of_first_order_mean_field_suppression"] for r in rows])
a.loglog(m, np.maximum(fr, 1e-16), color=C[0])
w0, w1 = GR["wave_window_meV"]
a.axvspan(5 * GR["trigger"]["m_phi_meV"], w1, color=C[2], alpha=0.2, label="wave dark matter, mean field valid")
a.axvspan(GR["trigger"]["classical_max_mass_meV"], 1e15, color=C[1], alpha=0.12, label="particle dark matter")
a.set_xlim(1e-4, 1e15); a.set_ylim(1e-14, 3)
a.set_xlabel(r"$m_{dm}$ (meV)  [as point particles]"); a.set_ylabel("suppression / mean-field suppression")
a.legend(fontsize=6.0, loc="lower left")
a.set_title("(b) Decay rate in a granular background", fontsize=8.5, fontweight="bold")
fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig3_resolved_granularity.png"), dpi=300, bbox_inches="tight", facecolor="white"); plt.close(fig)

# ---------------------------------------------------------------- wave wall
WW = J("wave_wall_1d.json")
P = WW["params"]
runs = {(r["mode"], r["m_ratio"], r["dz"]): r for r in WW["runs"]}
fig, ax = plt.subplots(1, 2, figsize=(7.6, 3.1))
a = ax[0]
for key, lab, col in ((("none", 10.0, 0.005), "no dark matter", C[4]), (("mean field", 10.0, 0.005), "mean-field term", C[6]),
                      (("wave", 10.0, 0.005), r"wave dark matter, $m=10m_\varphi$", C[0]), (("wave", 3.0, 0.005), r"wave dark matter, $m=3m_\varphi$", C[1])):
    r = runs[key]; a.plot(r["t"][:-1], r["gamma"][:-1], color=col, label=lab)
a.set_xlabel(r"time ($1/\sqrt{\lambda}f$)"); a.set_ylabel(r"wall Lorentz factor $\gamma$"); a.legend(fontsize=6.2)
a.set_title("(a) Wall with backreaction", fontsize=8.5, fontweight="bold")
a = ax[1]
r = runs[("wave", 10.0, 0.005)]
rho1 = P["rho_star"] * (1 + P["kappa"] * (P["x_h"] - 1) ** 2 / 2)
zw = np.array(r["z_wall_E"]); E = np.array(r["E_behind"])
a.plot(zw, E / rho1, color=C[0], label="simulation (10 $m_\\varphi$)")
z0 = zw[len(zw) // 2]; E0 = np.interp(z0, zw, E / rho1)
for k, lab, col in ((P["adiabatic"], "adiabatic, r", C[2]), (P["sudden"], "sudden, (1+r$^2$)/2", C[1])):
    a.plot(zw, E0 + k * (zw - z0), color=col, ls="--", lw=0.9, label=lab)
a.set_xlabel("wall position"); a.set_ylabel(r"dark-matter energy behind the wall / $\rho_{dm}$")
a.legend(fontsize=6.2, loc="upper left")
a.set_title("(b) Energy given to the dark matter", fontsize=8.5, fontweight="bold")
fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig6_wave_wall.png"), dpi=300, bbox_inches="tight", facecolor="white"); plt.close(fig)
print("figures written")
