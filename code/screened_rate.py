# -*- coding: utf-8 -*-
"""Granular decay rate with self-screened sources.

nonpert_bounce.py shows (i) that the nonperturbative increase of the bounce action from a source at the bubble centre is below the first-order
estimate for every coupling, and (ii) that the static coupling of a source of radius r_s saturates at about 4 pi r_s (thin-shell screening),
q_eff = q / (1 + q/(4 pi r_s)). A nonrelativistic particle is spread over its de Broglie length, r_s = 1/(m v), which is the largest physical size
it can have and so the least screened case. This script recomputes the Poisson rate of granular_field.py with q_eff in place of q, for the particle
masses of Table S7, and compares the effective cell-averaged coupling n q_eff with the mean-field coupling. Writes screened_rate.json.

    python screened_rate.py
"""
import json
import os

import numpy as np

import granular_field as GFm

HERE = os.path.dirname(os.path.abspath(__file__))
J = lambda f: json.load(open(os.path.join(HERE, f), encoding="utf-8"))
NB, GR = J("nonpert_bounce.json"), J("granularity.json")
V_DM = 1e-3
UNIT = GFm.UNIT_MEV
rho_t = GR["trigger"]["rho_dm_meV4"]
# check of the interpolation against the computed saturation
chk = []
for rs, rows in NB["screening"].items():
    rs = float(rs)
    for r in rows:
        chk.append(dict(rs=rs, q=r["q"], q_eff=r["q_eff"], interp=r["q"] / (1 + r["q"] / (4 * np.pi * rs))))
dev = max(abs(c["q_eff"] / c["interp"] - 1) for c in chk if c["rs"] <= 0.1)
rows = []
for m in np.logspace(np.log10(GR["trigger"]["classical_max_mass_meV"]), 15, 30):
    n = rho_t / m / UNIT ** 3
    q = GFm.TH_T / n
    rs = UNIT / (m * V_DM)                                     # de Broglie length in units of 1/(sqrt(lambda) f)
    q_eff = q / (1 + q / (4 * np.pi * rs))
    rows.append(dict(m_dm_meV=float(m), q=float(q), r_s=float(rs), q_eff=float(q_eff), screening=float(q_eff / q),
                     theta_eff_over_mean=float(q_eff / q), ln_rate_ratio_screened=GFm.ln_rate_ratio(n, q_eff),
                     ln_rate_ratio_unscreened=GFm.ln_rate_ratio(n, q)))
out = dict(interpolation_max_dev_small_rs=float(dev), check=chk, rows=rows,
           saturated_db_np={("%.2f" % b["rs"]): b["db_np"] for b in NB["bounces"] if b["q"] >= 1e4},
           mean_field_first_order=GFm.ln_rate_ratio(1e12, GFm.TH_T / 1e12))
json.dump(out, open(os.path.join(HERE, "screened_rate.json"), "w"), indent=1)
print("interpolation max deviation (r_s <= 0.1): %.3f" % dev)
for r in rows[::4]:
    print("m %.2e meV  q %.2e  r_s %.2e  q_eff/q %.2e  ln(G/G0) screened %.3e  unscreened %.3e" % (
        r["m_dm_meV"], r["q"], r["r_s"], r["screening"], r["ln_rate_ratio_screened"], r["ln_rate_ratio_unscreened"]))
