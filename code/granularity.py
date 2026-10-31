# -*- coding: utf-8 -*-
"""When does the scalar feel a smooth dark-matter density? Granularity of the source and the environment of the transition.

The effective potential V0 + rho_* A(phi) is a mean-field statement. The field responds to dark matter within its Compton length 1/m_phi, and on
the time scale 1/m_phi. Three kinds of coupled dark matter are distinguished:

  particle dark matter (occupation per de Broglie cell << 1): the source is a set of points; the mean field applies only if many lie inside a Compton
      volume, n >> m_phi^3, i.e. m_dm << rho/m_phi^3. Being particle-like needs m_dm >> m_cl = (rho / v^3)^(1/4). Both together need rho >> m_phi^4 v^-3/4 ...
      in practice m_cl > rho/m_phi^3 whenever rho < m_phi^4, so the window is empty. Between the particles the field sees zero density and the held
      vacuum decays at the zero-density rate.
  wave dark matter lighter than the field (m_dm << m_phi): the field follows the oscillation of the dark-matter field, whose instantaneous source
      passes through zero twice per period; tunneling proceeds in those intervals at the zero-density rate.
  wave dark matter heavier than the field (m_phi << m_dm << m_cl): the source oscillates faster than the field can respond, is coherent over the de
      Broglie length 1/(m_dm v) >> 1/m_phi, and the field feels its time average, rho_* A(phi). The mean field applies.

For the third kind the local density matters: theta_loc = theta_mean (1 + delta). The held vacuum is the true vacuum where 1 + delta > theta_deg /
theta_mean, so overdense regions keep it while the voids decay, and the past light cone must be protected in the emptiest regions. Writes
granularity.json.

    python granularity.py
"""
import json
import os

import numpy as np

import bubble_to_cycle as BC
import thermal as TH
import vacuum_wall as W

HERE = os.path.dirname(os.path.abspath(__file__))
J = lambda f: json.load(open(os.path.join(HERE, f), encoding="utf-8"))
WJ = J("vacuum_wall.json")
PS, TW, BASE = WJ["physical"], WJ["trigger_wall"], WJ["base"]
OMC, OML = 0.2607, 0.6889                       # Planck 2018 cold dark matter and dark energy
HBARC_M_MEV = 1.97327e-4                        # 1/meV in metres
H0_MEV = 1.4378e-30                             # H0 (h = 0.674) in meV
V_DM = 1e-3                                     # velocity dispersion of dark matter in haloes (c = 1); smaller in voids
lam, f = PS["lam"], PS["f_meV"]
q = BASE["q"]
Vf = lam * f ** 4 * PS["dtil0"] / (1.0 + q)     # meV^4, the observed dark-energy density
a_t = WJ["a_trigger"]
TH_DEG = WJ["walls"]["0.30"]["theta_deg"]


def ppp(x, th):
    return (3 * x * x - 1) / 2.0 + th


def rho_dm(a):
    return OMC / OML * Vf * a ** -3


def m_phi(x, th):
    return np.sqrt(lam * f * f * ppp(x, th))


def windows(a, x_h, th):
    rho, m = rho_dm(a), m_phi(x_h, th)
    return dict(a=a, rho_dm_meV4=rho, m_phi_meV=m, compton_length_mm=1e3 * HBARC_M_MEV / m, rho_over_mphi4=rho / m ** 4,
                particle_max_mass_meV=rho / m ** 3,                       # n >> m_phi^3
                classical_max_mass_meV=(rho / V_DM ** 3) ** 0.25,         # occupation per de Broglie cell >> 1
                coherence_max_mass_meV=m / V_DM)                          # de Broglie length >> Compton length


out = {}
# ------------------------------------------------------------------------------------------- (1) granularity at the trigger and today
w_t = windows(a_t, TW["x_f"], TW["theta"])
r0 = WJ["walls"]["0.30"]["rows"]
th_today = BASE["f0"] * TH_DEG
x_today = float(np.interp(BASE["f0"], [r["theta_over_deg"] for r in r0], [r["x_f"] for r in r0]))
w_0 = windows(1.0, x_today, th_today)
out["trigger"], out["today"] = w_t, w_0
out["particle_window_empty"] = bool(w_t["particle_max_mass_meV"] < w_t["classical_max_mass_meV"])
out["wave_window_meV"] = [w_t["m_phi_meV"], min(w_t["classical_max_mass_meV"], w_t["coherence_max_mass_meV"])]
out["spacing_over_compton"] = {}
for name, m in (("1 meV", 1.0), ("1 eV", 1e3), ("1 MeV", 1e9), ("1 GeV", 1e12), ("100 GeV", 1e14)):
    sp = (m / w_t["rho_dm_meV4"]) ** (1 / 3.0) * HBARC_M_MEV
    out["spacing_over_compton"][name] = dict(spacing_m=sp, over_compton=sp / (HBARC_M_MEV / w_t["m_phi_meV"]))

# structural estimate: rho/m_phi^4 = (rho/lambda f^4) / (lambda p''^2), with lambda = b/B fixed by completion (B ~ B_trigger)
scan = []
for eta, wl in WJ["walls"].items():
    d0 = wl["rows"][0]["dV"]
    for r in wl["rows"]:
        for qq in (0.1, 1.0, 3.0):
            l_ = r["b"] / TW["B"]
            ratio = (OMC / OML) * a_t ** -3 * d0 / (1.0 + qq) / (l_ * ppp(r["x_f"], r["theta"]) ** 2)
            scan.append(dict(eta=float(eta), theta_over_deg=r["theta_over_deg"], q=qq, lam=l_, ratio=ratio))
out["structural_scan_max"] = max(scan, key=lambda s: s["ratio"])
out["structural_scan_n"] = len(scan)
# approach to eta_c (the held minimum becomes shallow and the bounce action small, so lambda = b/B and p'' both fall): is there a corner with
# rho >> m_phi^4? Particle dark matter needs rho/m_phi^4 > 1/v (n >> m_phi^3 together with occupation << 1)
import vacuum_phase_map as VP
near = []
for eta in (0.36, 0.37, 0.375, 0.38, 0.383, 0.384):
    d0 = VP.classify_tb(eta, 0.0)[4]
    for fr in (0.0, 0.1, 0.3):
        th = fr * eta / 2.0
        bp = W.bounce_profile(eta, th)
        if bp is None:
            continue
        l_ = bp["b"] / TW["B"]
        for qq in (0.1, 1.0):
            ratio = (OMC / OML) * a_t ** -3 * d0 / (1.0 + qq) / (l_ * ppp(bp["x_f"], th) ** 2)
            near.append(dict(eta=eta, theta_over_deg=fr, q=qq, b=bp["b"], barrier=bp["barrier"], lam=l_, ratio=ratio))
out["near_eta_c"] = near
out["near_eta_c_max"] = max(near, key=lambda s: s["ratio"])
out["particle_needs_ratio_above"] = 1.0 / V_DM

# ------------------------------------------------------------------------------------------- (2) what particle dark matter implies
z_lc, w_lc = TH.light_cone_grid()
tabs = BC.rate_tables(WJ, "0.30")
code0, lnI0, _ = BC.classify_point(BASE["B0"], 1e-9, q, tabs, z_lc, w_lc)
lam0 = tabs[4] / BASE["B0"]
out["particle_dm"] = dict(ln_gamma_over_H4_zero_density=float(W.ln_gamma_h4(lam0, np.exp(tabs[2][0]), np.exp(tabs[3][0]), q)),
                          ln_expected_bubbles_past_light_cone=lnI0, note="the voids between particles see theta = 0 at all times")

# ------------------------------------------------------------------------------------------- (3) environment: wave dark matter heavier than the field
th_high = WJ["walls"]["0.30"]["rows"][-1]["theta"]                     # largest theta of the table (the window edge is in vacuum_phase_map.json)
PM = J("vacuum_phase_map.json")
th_high = PM["tunneling_profiles"]["0.30"]["window"][1]
out["environment"] = dict(theta_high_over_deg=th_high / TH_DEG,
                          held_true_above_1_plus_delta_at_trigger=1.0 / (BASE["f0"] * a_t ** -3),
                          only_held_minimum_above_1_plus_delta_at_trigger=(th_high / TH_DEG) / (BASE["f0"] * a_t ** -3),
                          held_true_above_1_plus_delta_today=1.0 / BASE["f0"])
# protection in the voids: the baseline evaluated with the local theta = theta_mean (1 + delta_v) along the whole past light cone (conservative:
# voids were less empty in the past), and the smallest theta_0 that keeps the past clear in such voids
voids = []
for d in (0.0, -0.5, -0.7, -0.8, -0.9):
    code, lnI, a_tr = BC.classify_point(BASE["B0"], BASE["f0"] * (1 + d), q, tabs, z_lc, w_lc)
    f_lo, f_hi = 1e-3, 30.0
    for _ in range(60):                      # smallest mean theta_0 / theta_deg with a clear past light cone in voids of this depth
        fm = np.sqrt(f_lo * f_hi)
        if BC.classify_point(BASE["B0"], fm * (1 + d), q, tabs, z_lc, w_lc)[1] < 0:
            f_hi = fm
        else:
            f_lo = fm
    c2, lnI2, a2 = BC.classify_point(BASE["B0"], f_hi * (1 + d), q, tabs, z_lc, w_lc)
    voids.append(dict(one_plus_delta=1 + d, baseline_lnI=lnI, baseline_past_clear=bool(lnI < 0), baseline_a_trigger_local=a_tr,
                      min_f0_mean=f_hi, a_trigger_local_at_min_f0=a2))
out["voids"] = voids
# curves for the figure: expected number of bubbles in the past light cone of a region of local density (1 + delta) theta_mean, for three means
curves = {}
grid = np.linspace(0.05, 1.0, 39)
for fmean in (0.3, 1.0, 3.0):
    curves["%.1f" % fmean] = [float(BC.classify_point(BASE["B0"], fmean * g, q, tabs, z_lc, w_lc)[1]) for g in grid]
out["void_curves"] = dict(one_plus_delta=grid.tolist(), lnI=curves)
# completion void by void: each region of size R must nucleate its own bubble; the rate grows at beta = d ln Gamma/dt, so a void of comoving radius R
# (R H << 1) completes after the global criterion Gamma/H^4 = 9/4pi once ln Gamma/H^4 has grown by a further ln[(9/4pi) / ((4pi/3)(R H)^3)]
HJ = J("bubble_handoff.json")
beta = HJ["percolation"]["beta_over_H"]
out["void_completion"] = {}
for R_mpc in (5.0, 20.0, 50.0):
    x = 4 * np.pi / 3 * (R_mpc / 4431.0) ** 3
    out["void_completion"]["%.0f Mpc" % R_mpc] = dict(extra_efolds_of_rate=float(np.log(9 / (4 * np.pi)) - np.log(x)),
                                                      delay_hubble_times=float((np.log(9 / (4 * np.pi)) - np.log(x)) / beta))
json.dump(out, open(os.path.join(HERE, "granularity.json"), "w"), indent=1)
for k, v in out.items():
    print(k, v)
