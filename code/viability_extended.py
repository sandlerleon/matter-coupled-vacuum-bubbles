# -*- coding: utf-8 -*-
"""Extended parameter survey and the environmental prediction.

Survey: each draw specifies eta, B0, the local coupling today in the emptiest regions theta_0/theta_deg, q, the depth of those regions (1 + delta_v),
and the dark-matter mass. Conditions (all evaluated in the emptiest regions, where the decay begins):
  T1 past light cone clear          T2 transition completes          T3 recollapse (q > 0)
  T4 mean field valid               m_dm >= 5 m_phi (time averaging, residual modulation <= 3%) and m_dm <= m_cl = (rho_loc / v^3)^(1/4), v = 1e-3
  T5 wall runaway with backreaction the energy given to the dark matter on a sudden crossing is below the driving pressure
Two sets of sampling distributions (a broad one and one concentrated on the region of interest) show how the fractions depend on the priors.

Prediction: the scale factor at which a region of density (1 + delta) times the mean completes its transition, for three mean couplings, with the
completion threshold shifted by +-ln((B/2pi)^2) to represent the uncertainty of the one-loop prefactor. Writes viability_extended.json.

    python viability_extended.py
"""
import json
import os
from multiprocessing import Pool

import numpy as np

import bubble_to_cycle as BC
import thermal as TH
import vacuum_phase_map as VP
import vacuum_wall as W

HERE = os.path.dirname(os.path.abspath(__file__))
WJ = json.load(open(os.path.join(HERE, "vacuum_wall.json"), encoding="utf-8"))
OMC, OML, VF = 0.2607, 0.6889, None
V_DM = 1e-3
PS = WJ["physical"]
VF = PS["lam"] * PS["f_meV"] ** 4 * PS["dtil0"] / (1.0 + WJ["base"]["q"])        # observed dark energy in meV^4
ETAS = ("0.20", "0.30", "0.36")
_T = {}


def init():
    global Z_LC, W_LC
    Z_LC, W_LC = TH.light_cone_grid()
    for e in ETAS:
        _T[e] = BC.rate_tables(WJ, e)


def draw(rng, kind, n):
    eta = rng.choice(ETAS, n)
    lu = lambda lo, hi: np.exp(rng.uniform(np.log(lo), np.log(hi), n))
    if kind == "broad":
        return dict(eta=eta, B0=lu(50.0, 1000.0), f0=lu(1e-3, 10.0), q=rng.uniform(-1.0, 3.0, n), dv=lu(0.05, 0.5), m=lu(1e-3, 1e12))
    return dict(eta=eta, B0=rng.uniform(150.0, 400.0, n), f0=lu(1e-2, 3.0), q=lu(0.1, 3.0), dv=rng.uniform(0.1, 0.3, n), m=lu(1.0, 1e3))


def evaluate(p):
    eta, B0, f0, q, dv, m = p
    tabs = _T[eta]
    code, lnI, a_t = BC.classify_point(B0, f0, max(q, -0.99), tabs, Z_LC, W_LC)
    r = dict(T1=bool(lnI < 0), T2=a_t is not None, T3=bool(q > 0), T4=False, T5=False, eta=float(eta), B0=float(B0), f0=float(f0), q=float(q),
             dv=float(dv), m=float(m), a_t=a_t)
    if a_t is None or q <= -0.99:
        return r
    e = float(eta)
    th_deg, fr, lb, ld, b0 = tabs
    th = f0 * th_deg / a_t ** 3
    n_min, x_t, x_top, x_h, dV, bar = VP.classify_tb(e, th)
    if n_min < 2:
        return r
    d0 = VP.classify_tb(e, 0.0)[4]
    lam = b0 / B0
    lf4 = (1.0 + q) * VF / d0                                  # lambda f^4 in meV^4
    f = (lf4 / lam) ** 0.25
    m_phi = np.sqrt(lam * f * f * VP.d2p_tb(x_h, e, th))
    rho_loc = OMC / OML * VF * dv * a_t ** -3                 # local dark-matter density at the trigger
    m_cl = (rho_loc / V_DM ** 3) ** 0.25
    r["T4"] = bool(5.0 * m_phi <= m <= m_cl)
    kappa = th / (rho_loc / lf4)                               # theta = kappa rho_*/(lambda f^4), with rho_* ~ rho_loc / A_h ~ rho_loc
    A = lambda x: 1.0 + kappa * (x - 1.0) ** 2 / 2.0
    rr = A(x_t) / A(x_h)
    deposit = (rho_loc / lf4) * ((1.0 + rr * rr) / 2.0 - 1.0)
    r["T5"] = bool(deposit < dV)
    r.update(m_phi=float(m_phi), m_cl=float(m_cl), deposit_over_dV=float(deposit / dV))
    return r


def survey(kind, n, seed):
    rng = np.random.default_rng(seed)
    d = draw(rng, kind, n)
    args = list(zip(d["eta"], d["B0"], d["f0"], d["q"], d["dv"], d["m"]))
    with Pool(8, initializer=init) as pool:
        res = pool.map(evaluate, args, chunksize=20)
    M = {k: np.array([x[k] for x in res]) for k in ("T1", "T2", "T3", "T4", "T5")}
    bg = M["T1"] & M["T2"] & M["T3"]
    out = dict(kind=kind, n=n, single={k: float(v.mean()) for k, v in M.items()}, joint_background=float(bg.mean()),
               joint_background_meanfield=float((bg & M["T4"]).mean()), joint_all=float((bg & M["T4"] & M["T5"]).mean()),
               runaway_given_background=float(M["T5"][bg].mean()) if bg.any() else None)
    dvs = d["dv"]
    edges = np.quantile(dvs, [0, 1 / 3, 2 / 3, 1])
    out["joint_all_by_void_depth"] = [dict(range=[float(edges[i]), float(edges[i + 1])],
                                           joint_all=float(((bg & M["T4"] & M["T5"])[(dvs >= edges[i]) & (dvs <= edges[i + 1])]).mean()))
                                      for i in range(3)]
    ok = [x for x, a, b, c in zip(res, bg, M["T4"], M["T5"]) if a and b and c]
    ph = {}
    for k, lab in (("B0", "B0"), ("f0", "local theta_0/theta_deg"), ("q", "q"), ("dv", "1+delta_v"), ("m", "m_dm (meV)"), ("m_phi", "m_phi (meV)"),
                   ("a_t", "a at completion (emptiest regions)")):
        v = np.array([x[k] for x in ok if x.get(k) is not None], dtype=float)
        if v.size:
            ph[lab] = dict(min=float(v.min()), q10=float(np.quantile(v, 0.1)), median=float(np.median(v)), q90=float(np.quantile(v, 0.9)),
                           max=float(v.max()))
    out["viable_physical_ranges"] = ph
    out["viable_by_eta"] = {e: int(sum(1 for x in ok if abs(x["eta"] - float(e)) < 1e-9)) for e in ETAS}
    out["n_viable"] = len(ok)
    dep = [x["deposit_over_dV"] for x in res if "deposit_over_dV" in x]
    out["deposit_over_dV_max"] = float(max(dep)) if dep else None
    return out


def prediction():
    init()
    tabs = _T["0.30"]
    B0, q = WJ["base"]["B0"], WJ["base"]["q"]
    lnB = 2.0 * np.log(WJ["trigger_wall"]["B"] / (2 * np.pi))
    base_thr = W.LN_COMPLETE
    rows = []
    for fmean in (1.0, 2.0, 3.0):
        for g in (0.1, 0.2, 0.3, 0.5, 1.0, 2.0):
            row = dict(mean_f0=fmean, one_plus_delta=g)
            for lab, sh in (("central", 0.0), ("prefactor_up", -lnB), ("prefactor_down", +lnB)):
                W.LN_COMPLETE = base_thr + sh
                code, lnI, a_t = BC.classify_point(B0, fmean * g, q, tabs, Z_LC, W_LC)
                row[lab] = a_t
                if lab == "central":
                    row["past_clear"] = bool(lnI < 0)
            W.LN_COMPLETE = base_thr
            rows.append(row)
    return rows


if __name__ == "__main__":
    out = {"surveys": [survey("broad", 3000, 20261031), survey("focused", 3000, 20261101)], "prediction": prediction()}
    json.dump(out, open(os.path.join(HERE, "viability_extended.json"), "w"), indent=1)
    for s in out["surveys"]:
        print(s["kind"], s["single"], "background %.3f +meanfield %.4f all %.4f runaway|bg %s max deposit/dV %s" % (
            s["joint_background"], s["joint_background_meanfield"], s["joint_all"], s["runaway_given_background"], s["deposit_over_dV_max"]))
        print("   by void depth", s["joint_all_by_void_depth"])
        print("   viable", s["n_viable"], s["viable_by_eta"])
        for k, v in s["viable_physical_ranges"].items():
            print("     %-36s %s" % (k, {kk: round(vv, 4) for kk, vv in v.items()}))
    for r in out["prediction"]:
        print(r)
