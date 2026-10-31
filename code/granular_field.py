# -*- coding: utf-8 -*-
"""Spatially resolved tests of the granularity bound (Section 3): the field around discrete dark-matter sources, and the decay rate in a granular
background.

(1) Static field in a Wigner-Seitz cell. Each dark-matter particle is a source of the matter term; in the units of the paper the static field obeys
        x'' + (2/r) x' = p0'(x) + theta_s(r) (x - 1),     x'(0) = x'(R_c) = 0,
    with theta_s a uniform sphere of radius r_s carrying the charge q = theta_mean (4 pi/3) R_c^3 of one cell (so that the cell average of the source is
    the mean coupling). Solved by collocation for cells from much smaller to much larger than the Compton length; the result is the fraction of the cell
    in which the field is within 1% of the zero-density vacuum and the cell-averaged field against the mean-field minimum.

(2) Decay rate in a granular background. A static particle at distance r_p from the centre of the O(4) bounce adds to the Euclidean action, at first
    order in its coupling, dB = (q/lambda) I(r_p),  I(r_p) = int dtau [(x_b(rho) - 1)^2 - (x_h - 1)^2]/2,  rho = sqrt(r_p^2 + tau^2), with x_b the
    zero-density bounce. For Poisson-distributed particles of dimensionless density n the rate is
        Gamma / Gamma_0 = < exp(-sum dB) > = exp( - n int d^3r [1 - exp(-(q/lambda) I(r))] ),
    which reduces to the first-order mean-field result exp(-theta_mean db/dtheta / lambda) for many weak sources (checked against the numerical
    derivative of the bounce action) and to Gamma_0 for few strong ones. The first-order treatment overstates the suppression from a strong source (the
    field near the particle saturates and the bounce can deform around it), so for heavy particles the true rate is at least this close to Gamma_0.
Writes granular_field.json.

    python granular_field.py
"""
import json
import os

import numpy as np
from scipy.integrate import solve_bvp, quad
from scipy.interpolate import interp1d

import thermal as TH
import vacuum_phase_map as VP
import vacuum_wall as W

HERE = os.path.dirname(os.path.abspath(__file__))
J = lambda f: json.load(open(os.path.join(HERE, f), encoding="utf-8"))
WJ, GR = J("vacuum_wall.json"), J("granularity.json")
ETA = WJ["base"]["eta"]
TW, PS = WJ["trigger_wall"], WJ["physical"]
LAM = PS["lam"]
TH_T = TW["theta"]                                   # mean coupling at the trigger
UNIT_MEV = np.sqrt(LAM) * PS["f_meV"]                # 1 unit of mass = sqrt(lambda) f in meV
dp0 = lambda x: VP.dp_tb(x, ETA, 0.0)


def x_min(theta, which="held"):
    st = VP.stationary_tb(ETA, theta)
    mins = sorted(x for x in st if VP.d2p_tb(x, ETA, theta) > 0)
    return mins[-1] if which == "held" else mins[0]


# ------------------------------------------------------------------------------------------- (1) static field in a cell
def cell(Rc, theta_mean=TH_T, rs_frac=0.05, rs_max=0.3):
    rs = min(rs_frac * Rc, rs_max)
    q = theta_mean * 4.0 / 3.0 * np.pi * Rc ** 3
    th_in = q / (4.0 / 3.0 * np.pi * rs ** 3)
    xh0, xhm = x_min(0.0), x_min(theta_mean)
    # grid refined near the source
    r = np.unique(np.concatenate([np.linspace(1e-6, rs * 1.5, 300), np.geomspace(rs * 1.5, Rc, 700)]))
    src = lambda rr: th_in * 0.5 * (1.0 - np.tanh((rr - rs) / (0.05 * rs)))

    def f(rr, y):
        return np.vstack([y[1], dp0(y[0]) + src(rr) * (y[0] - 1.0) - 2.0 * y[1] / rr])

    bc = lambda ya, yb: np.array([ya[1], yb[1]])
    guess = np.vstack([np.where(r < rs, 1.0, xh0), np.zeros_like(r)])
    sol = solve_bvp(f, bc, r, guess, tol=1e-7, max_nodes=200000)
    rr = np.linspace(1e-6, Rc, 20001)
    x = sol.sol(rr)[0]
    w = 4 * np.pi * rr ** 2
    vol = np.trapezoid(w, rr)
    near0 = np.trapezoid(w * (np.abs(x - xh0) < 0.01 * abs(1 - xh0)), rr) / vol
    nearm = np.trapezoid(w * (np.abs(x - xhm) < 0.01 * abs(1 - xh0)), rr) / vol
    xavg = np.trapezoid(w * x, rr) / vol
    return dict(Rc=Rc, Rc_times_mphi=Rc * np.sqrt(VP.d2p_tb(xh0, ETA, 0.0)), rs=rs, success=bool(sol.status == 0), x_center=float(x[0]),
                x_edge=float(x[-1]), x_cell_average=float(xavg), x_h_zero_density=xh0, x_h_mean_field=xhm,
                fraction_at_zero_density_vacuum=float(near0), fraction_at_mean_field_vacuum=float(nearm),
                profile_r=rr[::200].tolist(), profile_x=x[::200].tolist())


# ------------------------------------------------------------------------------------------- (2) rate in a granular background
def bounce0():
    """Zero-density O(4) bounce of eta = 0.30 at full resolution: x_b(rho), b(0) and the first-order density derivative db/dtheta."""
    n_min, x_t, x_top, x_f, dV, bar = VP.classify_tb(ETA, 0.0)
    v = lambda x: VP.p_tb(x, ETA, 0.0)
    dv = lambda x: VP.dp_tb(x, ETA, 0.0)
    b, x0 = TH.bounce_shoot(v, dv, x_t, x_top, x_f, 4)
    _, sol = TH._shoot(x0, dv, x_f, 4)
    rho = np.linspace(0.0, sol.t[-1], 40001)
    xb = sol.sol(np.maximum(rho, 1e-9))[0]
    # tail beyond the integration range: exponential approach to x_f with the false-vacuum mass
    return dict(b=float(b), x_f=float(x_f), rho=rho, x=xb, m_f=float(np.sqrt(VP.d2p_tb(x_f, ETA, 0.0))))


BB = bounce0()
XB = interp1d(BB["rho"], BB["x"], bounds_error=False, fill_value=BB["x_f"])
G = lambda rho: 0.5 * ((XB(rho) - 1.0) ** 2 - (BB["x_f"] - 1.0) ** 2)     # source-term density of the bounce relative to the false vacuum


def I_of(rp):
    Tm = np.sqrt(max(BB["rho"][-1] ** 2 - rp ** 2, 0.0))
    if Tm <= 0:
        return 0.0
    t = np.linspace(0.0, Tm, 4001)
    return float(2.0 * np.trapezoid(G(np.sqrt(rp ** 2 + t ** 2)), t))


RP = np.linspace(0.0, BB["rho"][-1] * 0.999, 1500)
IP = np.array([I_of(r) for r in RP])
DBDTH_FIRST = float(4 * np.pi * np.trapezoid(RP ** 2 * IP, RP))            # = 2 pi^2 int rho^3 G drho = db/dtheta at theta = 0
DBDTH_DIRECT = float(2 * np.pi ** 2 * np.trapezoid(BB["rho"] ** 3 * G(BB["rho"]), BB["rho"]))


def db_numeric(h=2e-4):
    b1 = W.bounce_profile(ETA, h)["b"]
    return (b1 - BB["b"]) / h


def ln_rate_ratio(n, q):
    """ln(Gamma/Gamma_0) for Poisson sources of dimensionless density n and charge q (first order per source)."""
    with np.errstate(over="ignore"):
        e = -np.expm1(-(q / LAM) * IP)
    return float(-n * 4 * np.pi * np.trapezoid(RP ** 2 * e, RP))


if __name__ == "__main__":
    out = {}
    out["cells"] = [cell(Rc) for Rc in (0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0)]
    out["cells_rs_check"] = [cell(30.0, rs_frac=0.05, rs_max=r) for r in (0.1, 0.3, 1.0)]
    num = db_numeric()
    out["bounce"] = dict(b0=BB["b"], db_dtheta_first_order=DBDTH_FIRST, db_dtheta_direct=DBDTH_DIRECT, db_dtheta_numerical=num,
                         first_order_vs_numerical=DBDTH_FIRST / num - 1.0)
    rho_t = GR["trigger"]["rho_dm_meV4"]
    rows = []
    for m_meV in np.logspace(-4, 15, 39):
        n = rho_t / m_meV / UNIT_MEV ** 3                                  # particles per unit (sqrt(lambda) f)^-3
        q = TH_T / n
        lr = ln_rate_ratio(n, q)
        rows.append(dict(m_dm_meV=float(m_meV), n_units=float(n), q=float(q), ln_rate_ratio=lr,
                         spacing_over_compton=float(n ** (-1 / 3.0) * GR["trigger"]["m_phi_meV"] / UNIT_MEV)))
    mf_first = -TH_T * num / LAM
    mf_full = -(TW["b"] - BB["b"]) / LAM
    for r in rows:
        r["fraction_of_first_order_mean_field_suppression"] = r["ln_rate_ratio"] / mf_first
    out["rate_vs_mass"] = rows
    out["mean_field"] = dict(ln_rate_ratio_first_order=mf_first, ln_rate_ratio_full=mf_full)
    half = [r for r in rows if r["fraction_of_first_order_mean_field_suppression"] < 0.5]
    out["m_half_meV"] = half[0]["m_dm_meV"] if half else None
    pick = lambda m: min(rows, key=lambda r: abs(np.log(r["m_dm_meV"] / m)))
    out["particle_examples"] = {k: pick(m) for k, m in (("1 eV", 1e3), ("1 keV", 1e6), ("1 MeV", 1e9), ("1 GeV", 1e12))}
    json.dump(out, open(os.path.join(HERE, "granular_field.json"), "w"), indent=1)
    for c in out["cells"]:
        print("cell R_c m_phi = %8.2f  ok %s  zero-density fraction %.3f  mean-field fraction %.3f  <x> %.4f (x_h0 %.4f, x_hm %.4f)" % (
            c["Rc_times_mphi"], c["success"], c["fraction_at_zero_density_vacuum"], c["fraction_at_mean_field_vacuum"], c["x_cell_average"],
            c["x_h_zero_density"], c["x_h_mean_field"]))
    print("r_s check:", [round(c["fraction_at_zero_density_vacuum"], 4) for c in out["cells_rs_check"]])
    print("bounce:", out["bounce"])
    print("mean field: first order %.2f, full %.2f" % (mf_first, mf_full))
    for r in rows[::3]:
        print("m %.1e meV  spacing/Compton %.2e  ln(G/G0) %.4e  fraction of MF %.3e" % (r["m_dm_meV"], r["spacing_over_compton"], r["ln_rate_ratio"],
                                                                                    r["fraction_of_first_order_mean_field_suppression"]))
    print("m_half", out["m_half_meV"])
