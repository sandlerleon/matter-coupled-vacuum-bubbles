# -*- coding: utf-8 -*-
"""A planar vacuum wall running through wave dark matter, with backreaction, in 1+1 dimensions.

Fields (units of the paper: x = phi/f, lengths and times in 1/(sqrt(lambda) f), energies in lambda f^4; chi in units of f):
    x_tt = x_zz - p0'(x) - m^2 A A' chi^2,      chi_tt = chi_zz - m^2 A^2 chi,      A = 1 + kappa (x - 1)^2 / 2,
i.e. a dark-matter field of mass m A(x). Its time average gives the mean-field term theta (x - 1) with theta = kappa rho_*, so the mean-field run
replaces the dark matter by that static term; a third run has no dark matter. The wall starts at rest as a kink between the true and the held
vacuum and accelerates through a homogeneous oscillating dark-matter background; nothing is prescribed about the crossing.

Measured: the wall trajectory and Lorentz factor in the three runs (backreaction of the energy transfer on the wall); the dark-matter energy
deposited per unit length swept by the wall, against the adiabatic (r rho) and sudden ((1 + r^2)/2 rho) predictions, r = A(x_t)/A(x_h); the mean
velocity of the dark matter behind the wall (momentum over energy), against (r^2 - 1)/(r^2 + 1) for quanta crossing a relativistic wall; and the
energy density at probe positions just behind the wall (a trailing shell compressed by 1/(1 - v_dm)). Two dark-matter masses (3 and 10 m_phi) and
two resolutions. Writes wave_wall_1d.json.

    python wave_wall_1d.py
"""
import json
import os
import time
from multiprocessing import Pool

import numpy as np

import vacuum_phase_map as VP

HERE = os.path.dirname(os.path.abspath(__file__))
WJ = json.load(open(os.path.join(HERE, "vacuum_wall.json"), encoding="utf-8"))
TJ = json.load(open(os.path.join(HERE, "interface_transport.json"), encoding="utf-8"))
ETA = WJ["base"]["eta"]
TH = WJ["trigger_wall"]["theta"]
XH0 = TJ["wall"]["x_exterior"]
A_F = TJ["dark_matter"]["A_false"]
KAPPA = 2.0 * (A_F - 1.0) / (XH0 - 1.0) ** 2         # A(x_h) of the baseline
RHO_S = TH / KAPPA                                    # rho_* in units of lambda f^4
dp0 = lambda x: x * (x * x - 1.0) / 2.0 + ETA / 2.0


def mins(theta):
    st = VP.stationary_tb(ETA, theta)
    m = sorted(x for x in st if VP.d2p_tb(x, ETA, theta) > 0)
    return m[0], m[-1]


X_T, X_H = mins(TH)
A_ = lambda x: 1.0 + KAPPA * (x - 1.0) ** 2 / 2.0
R_RATIO = A_(X_T) / A_(X_H)
M_PHI = np.sqrt(VP.d2p_tb(X_H, ETA, TH))
SIGMA = WJ["walls"]["0.30"]["degenerate"]["sigma"]


def run(mode, m_ratio=10.0, dz=0.005, L=200.0, T=190.0, z0=8.0, probes=(40, 70, 100, 130, 160)):
    t0 = time.time()
    m = m_ratio * M_PHI
    z = np.arange(0.0, L + dz / 2, dz)
    dt = 0.5 * dz
    x = X_T + (X_H - X_T) * 0.5 * (1.0 + np.tanh((z - z0) / 2.0))
    xv = np.zeros_like(z)
    dm = mode == "wave"
    th_static = TH if mode == "mean field" else 0.0
    if dm:
        Az = A_(x)
        chi = np.sqrt(2.0 * RHO_S / Az) / m                  # amplitude for number density rho_*/m, at rest
        chiv = np.zeros_like(z)

    def lap(u):
        up = np.empty(u.size + 2)
        up[1:-1] = u
        up[0], up[-1] = u[1], u[-2]                          # reflecting ends
        return (up[2:] - 2.0 * u + up[:-2]) / dz ** 2

    def acc(x, chi=None):
        ax = lap(x) - dp0(x) - th_static * (x - 1.0)
        if chi is None:
            return ax, None
        A = A_(x)
        ax -= m * m * A * KAPPA * (x - 1.0) * chi * chi
        return ax, lap(chi) - m * m * A * A * chi

    pidx = [int(round(p / dz)) for p in probes]
    rec_t, rec_zw, rec_E, rec_P = [], [], [], []
    pe = {p: [] for p in probes}
    ax, ac = acc(x, chi if dm else None)
    nsteps = int(T / dt)
    xmid = 0.5 * (X_T + X_H)
    for k in range(nsteps):
        xv += 0.5 * dt * ax
        x += dt * xv
        if dm:
            chiv += 0.5 * dt * ac
            chi += dt * chiv
        ax, ac = acc(x, chi if dm else None)
        xv += 0.5 * dt * ax
        if dm:
            chiv += 0.5 * dt * ac
        if k % 20 == 0:
            i = int(np.argmax(x[int(z0 / dz):] > xmid)) + int(z0 / dz)
            if 0 < i < z.size:
                zw = z[i - 1] + (xmid - x[i - 1]) / (x[i] - x[i - 1]) * dz
                rec_t.append((k + 1) * dt); rec_zw.append(float(zw))
                if dm:
                    sel = slice(int(1.0 / dz), max(int((zw - 1.0) / dz), int(1.0 / dz) + 1))
                    A = A_(x[sel]); cz = np.gradient(chi, dz)[sel]
                    rec_E.append(float(np.sum(0.5 * chiv[sel] ** 2 + 0.5 * cz ** 2 + 0.5 * m * m * A * A * chi[sel] ** 2) * dz))
                    rec_P.append(float(np.sum(-chiv[sel] * cz) * dz))
            if dm:
                for p, ip in zip(probes, pidx):
                    A = A_(x[ip])
                    cz = (chi[ip + 1] - chi[ip - 1]) / (2 * dz)
                    pe[p].append(0.5 * chiv[ip] ** 2 + 0.5 * cz ** 2 + 0.5 * m * m * A * A * chi[ip] ** 2)
    t_ = np.array(rec_t); zw = np.array(rec_zw)
    # the wall follows hyperbolic motion, z = z0 + L (sqrt(1 + ((t - t0)/L)^2) - 1), with L = sigma / P_eff; the Lorentz factor is taken from this fit
    # (finite differences of the tracked position are too noisy once v -> 1)
    from scipy.optimize import least_squares
    fres = lambda q: zw - (q[0] + q[1] * (np.sqrt(1.0 + ((t_ - q[2]) / q[1]) ** 2) - 1.0))
    fit = least_squares(fres, [z0, 2.0, 0.0], loss="soft_l1", f_scale=0.05)
    Lf, tf = fit.x[1], fit.x[2]
    gam = np.sqrt(1.0 + ((t_ - tf) / Lf) ** 2)
    out = dict(mode=mode, m_ratio=m_ratio, dz=dz, t=t_[::10].tolist(), z_wall=zw[::10].tolist(), gamma=gam[::10].tolist(), L_fit=float(Lf),
               gamma_end_fit=float(gam[-1]), fit_rms=float(np.sqrt(np.mean(fres(fit.x) ** 2))), seconds=time.time() - t0)
    if dm:
        E, Pm = np.array(rec_E), np.array(rec_P)
        late = t_ > 0.5 * t_[-1]
        slope = np.polyfit(zw[late], E[late], 1)[0]
        out["deposit_per_length_over_rho"] = float(slope / (RHO_S * A_(X_H)))
        out["v_dm_behind_late"] = float(np.mean(Pm[late] / E[late]))
        out["E_behind"] = E[::10].tolist()
        out["z_wall_E"] = zw[::10].tolist()
        res = []
        e_before = RHO_S * A_(X_H)
        for p in probes:
            e = np.array(pe[p]); tt = np.arange(len(e)) * 20 * dt
            tc = float(np.interp(p, zw, t_))                    # wall crossing time
            gc = float(np.interp(p, zw, gam))
            sel_b = (tt > tc - 25.0) & (tt < tc - 3.0)
            sel_a = (tt > tc + 3.0) & (tt < min(tc + 25.0, tt[-1]))
            if sel_a.sum() < 20:
                continue
            sel_s = (tt > tc + 0.5) & (tt < tc + 2.5)                # just behind the wall: the trailing shell
            res.append(dict(z=p, gamma_at_crossing=gc, e_before_over_expected=float(e[sel_b].mean() / e_before),
                            ratio_after_before=float(e[sel_a].mean() / e[sel_b].mean()),
                            shell_over_before=float(e[sel_s].mean() / e[sel_b].mean()),
                            crossing_time_over_period=float(4.6 / gc * m * A_(X_H) / (2 * np.pi))))
        out["probes"] = res
    return out


def job(a):
    return run(*a[0], **a[1])


if __name__ == "__main__":
    t0 = time.time()
    jobs = [(("wave",), dict(m_ratio=10.0)), (("wave",), dict(m_ratio=3.0)), (("mean field",), {}), (("none",), {}),
            (("wave",), dict(m_ratio=10.0, dz=0.01)), (("mean field",), dict(dz=0.01))]
    with Pool(len(jobs)) as pool:
        runs = pool.map(job, jobs)
    out = dict(params=dict(eta=ETA, theta=TH, kappa=KAPPA, rho_star=RHO_S, x_t=X_T, x_h=X_H, r=R_RATIO, adiabatic=R_RATIO,
                           sudden=(1 + R_RATIO ** 2) / 2.0, v_quanta=(R_RATIO ** 2 - 1) / (R_RATIO ** 2 + 1),
                           shell_compression=(1 + R_RATIO ** 2) / 2.0 / (1 - (R_RATIO ** 2 - 1) / (R_RATIO ** 2 + 1)), m_phi=M_PHI,
                           sigma_degenerate=SIGMA), runs=runs)
    base = {(r["mode"], r["dz"], r["m_ratio"]): r for r in runs}
    mf = base[("mean field", 0.005, 10.0)]
    for r in runs:
        if r["mode"] == "wave" and r["dz"] == 0.005:
            r["gamma_rel_diff_vs_mean_field_end"] = float(r["gamma_end_fit"] / mf["gamma_end_fit"] - 1.0)
            r["P_eff_rel_diff_vs_mean_field"] = float(mf["L_fit"] / r["L_fit"] - 1.0)
    out["convergence"] = dict(
        wave_z_end_rel=float(base[("wave", 0.01, 10.0)]["z_wall"][-1] / base[("wave", 0.005, 10.0)]["z_wall"][-1] - 1.0),
        mean_field_z_end_rel=float(base[("mean field", 0.01, 10.0)]["z_wall"][-1] / mf["z_wall"][-1] - 1.0),
        wave_L_rel=float(base[("wave", 0.01, 10.0)]["L_fit"] / base[("wave", 0.005, 10.0)]["L_fit"] - 1.0),
        mean_field_L_rel=float(base[("mean field", 0.01, 10.0)]["L_fit"] / mf["L_fit"] - 1.0))
    json.dump(out, open(os.path.join(HERE, "wave_wall_1d.json"), "w"), indent=1)
    print("r = %.3f: adiabatic %.3f, sudden %.3f; m_phi %.3f" % (R_RATIO, R_RATIO, (1 + R_RATIO ** 2) / 2, M_PHI))
    for r in runs:
        print(r["mode"], r["m_ratio"], r["dz"], "z_end %.2f gamma_end(fit) %.1f L %.3f rms %.3f [%.0fs]" % (r["z_wall"][-1], r["gamma_end_fit"],
                                                                                                         r["L_fit"], r["fit_rms"], r["seconds"]),
              r.get("gamma_rel_diff_vs_mean_field_end", ""), "deposit/rho %s v_dm %s" % (r.get("deposit_per_length_over_rho"), r.get("v_dm_behind_late")))
        for p in r.get("probes", []):
            print("   z %5.0f gamma %6.1f  e_before/expected %.3f  after/before %.3f  shell %.2f  crossing/period %.3f" % (
                p["z"], p["gamma_at_crossing"], p["e_before_over_expected"], p["ratio_after_before"], p["shell_over_before"],
                p["crossing_time_over_period"]))
    print("convergence", out["convergence"], "[%.0fs]" % (time.time() - t0))
