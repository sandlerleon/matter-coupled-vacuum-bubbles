# -*- coding: utf-8 -*-
"""Validity of the mean-field (smooth-density) effective potential at the scale of the vacuum wall.

The field responds to dark matter within its Compton length 1/m_phi. The averaged potential V0 + rho*A applies only if many dark-matter quanta lie in
that volume (n >> m_phi^3) and each is heavier than m_phi (so the field cannot follow the oscillation of a dark-matter wave). Together these require
rho_dm >> m_phi^4. Writes meanfield_validity.json."""
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
J = lambda f: json.load(open(os.path.join(HERE, f), encoding="utf-8"))
WJ = J("vacuum_wall.json")
PS, TW = WJ["physical"], WJ["trigger_wall"]
OMC, OML = 0.2607, 0.6889                      # Planck 2018 cold dark matter and dark energy
HBARC_M_MEV = 1.97327e-4                       # 1/meV in metres
lam, f = PS["lam"], PS["f_meV"]
q = WJ["base"]["q"]
Vf = lam * f ** 4 * PS["dtil0"] / (1.0 + q)    # meV^4
a_t = WJ["a_trigger"]


def ppp(x, th):
    return (3 * x * x - 1) / 2.0 + th


def ratio(eta_rows_dtil, b_trig, x_h, th, a, qq, B_trig=TW["B"]):
    lam_ = b_trig / B_trig
    rho_over_lf4 = (OMC / OML) * a ** -3 * eta_rows_dtil / (1.0 + qq)
    return rho_over_lf4 / (lam_ * ppp(x_h, th) ** 2), lam_


rho_t = OMC / OML * Vf * a_t ** -3
m_phi = np.sqrt(lam * f * f * ppp(TW["x_f"], TW["theta"]))
out = dict(rho_dm_trigger_meV4=rho_t, m_phi_meV=m_phi, compton_length_mm=1e3 * HBARC_M_MEV / m_phi, rho_over_mphi4=rho_t / m_phi ** 4,
           m_dm_max_for_averaging_meV=rho_t / m_phi ** 3)
out["spacing_over_compton"] = {}
for name, m in (("1 meV", 1.0), ("1 MeV", 1e9), ("1 GeV", 1e12), ("100 GeV", 1e14)):
    sp = (m / rho_t) ** (1 / 3.0) * HBARC_M_MEV
    out["spacing_over_compton"][name] = dict(spacing_m=sp, over_compton=sp / (HBARC_M_MEV / m_phi))
# zero-density rate seen in the voids between dark-matter particles
B0 = WJ["base"]["B0"]
H0_meV = 1.4378e-30                             # H0 = 67.4 km/s/Mpc in meV
out["ln_gamma_over_H4_zero_density"] = 4 * np.log(m_phi / H0_meV) + 2 * np.log(B0 / (2 * np.pi)) - B0
# scan over tilt and q at the trigger with lambda fixed by B = B_trig
scan = []
for eta, w in WJ["walls"].items():
    for r in w["rows"]:
        for qq in (0.1, 1.0, 3.0):
            R, l_ = ratio(w["rows"][0]["dV"], r["b"], r["x_f"], r["theta"], a_t, qq)
            scan.append(dict(eta=float(eta), theta_over_deg=r["theta_over_deg"], q=qq, lam=l_, ratio=R))
out["scan_max_ratio"] = max(scan, key=lambda s: s["ratio"])
out["scan_n"] = len(scan)
out["scan_baseline_eta"] = [s for s in scan if s["eta"] == 0.30 and s["q"] == 1.0]
json.dump(out, open(os.path.join(HERE, "meanfield_validity.json"), "w"), indent=1)
for k, v in out.items():
    if k != "scan_baseline_eta":
        print(k, v)
