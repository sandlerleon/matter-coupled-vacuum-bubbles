# -*- coding: utf-8 -*-
"""Nonperturbative bounce with many discrete sources: a periodic row of sources along a line, against the same coupling smeared along the line.

Coordinates (tau, x, rho): sources are static (lines in Euclidean time) at x_k = (k + 1/2) d on the axis rho = 0, each a sphere of radius r_s with
coupling q, so the field depends on (tau, x, rho) only (axial symmetry about the row). The comparison configuration spreads the same coupling per
unit length, q/d, uniformly along a cylinder of radius r_s on the axis. The bounce is found as in nonpert_bounce.py: the maximum over the four-volume
beyond the barrier top of the constrained minimum of the action (L-BFGS, quadratic penalty, continuation), starting from the zero-density bounce
centred between two sources; the background is the static field of the sources. For d large the bubble fits between the sources and the action
returns to the zero-density value, while the smeared line, which carries the same mean coupling, keeps suppressing the decay: the discreteness of the
sources, not their mean density, controls the rate. Writes multisource_bounce.json.

    python multisource_bounce.py
"""
import json
import os
import time
from multiprocessing import Pool

import numpy as np
from scipy.optimize import minimize

import nonpert_bounce as NB

HERE = os.path.dirname(os.path.abspath(__file__))
X_H0, X_TOP = NB.X_H0, NB.X_TOP
p0, dp0 = NB.p0, NB.dp0
PROF = NB.WJ["profiles"]["0.00"]


class Grid3:
    def __init__(self, h=0.5, Lt=16.0, Lx=16.0, Lr=16.0, sources="none", d=None, q=0.0, rs=1.0):
        self.h = h
        self.nt, self.nx, self.nr = int(round(2 * Lt / h)), int(round(2 * Lx / h)), int(round(Lr / h))
        self.tau = -Lt + (np.arange(self.nt) + 0.5) * h
        self.x = -Lx + (np.arange(self.nx) + 0.5) * h
        self.r = (np.arange(self.nr) + 0.5) * h
        rf = np.arange(self.nr + 1) * h
        self.Vxr = h * np.pi * (rf[1:] ** 2 - rf[:-1] ** 2)[None, :] * np.ones((self.nx, 1))      # (x, rho) cell volume times dx
        self.Ar = 2 * np.pi * rf[1:] * h                                                          # outer rho face (per unit tau), times dx
        # source density by subsampling each (x, rho) cell
        th = np.zeros((self.nx, self.nr))
        if sources != "none":
            sub = (np.arange(6) + 0.5) / 6.0
            for i, xc in enumerate(self.x):
                xs = xc - h / 2 + sub * h
                for j, rc in enumerate(self.r):
                    rr = rf[j] + sub * h
                    X, R = np.meshgrid(xs, rr, indexing="ij")
                    w = R / R.sum()
                    if sources == "discrete":
                        k = np.round(X / d - 0.5)
                        dist2 = (X - (k + 0.5) * d) ** 2 + R ** 2
                        th[i, j] = q / (4.0 / 3.0 * np.pi * rs ** 3) * np.sum(w * (dist2 < rs ** 2))
                    else:
                        th[i, j] = (q / d) / (np.pi * rs ** 2) * np.sum(w * (R < rs))
        self.th = th
        self.ubg = self.background()

    def static(self, u):
        h = self.h
        U = u.reshape(self.nx, self.nr)
        Px = np.vstack([np.full((1, self.nr), X_H0), U, np.full((1, self.nr), X_H0)])
        dx = np.diff(Px, axis=0)
        Pr = np.hstack([U, np.full((self.nx, 1), X_H0)])
        dr = np.diff(Pr, axis=1)
        Vx = self.Vxr
        e = np.sum(0.5 * np.vstack([Vx, Vx[-1:]]) * dx ** 2) / h ** 2 + np.sum(self.Ar[None, :] * 0.5 * dr ** 2) / h \
            + np.sum(Vx * (p0(U) + 0.5 * self.th * (U - 1) ** 2))
        g = Vx * (dp0(U) + self.th * (U - 1))
        Vxe = np.vstack([Vx, Vx[-1:]])
        g += (Vxe[:-1] * dx[:-1] - Vxe[1:] * dx[1:]) / h ** 2
        g[:, :-1] += -self.Ar[None, :-1] * dr[:, :-1] / h
        g[:, 1:] += self.Ar[None, :-1] * dr[:, :-1] / h
        g[:, -1] += -self.Ar[-1] * dr[:, -1] / h
        return e, g.ravel()

    def background(self):
        res = minimize(self.static, np.full(self.nx * self.nr, X_H0), jac=True, method="L-BFGS-B",
                       options=dict(maxiter=50000, maxcor=30, gtol=1e-11, ftol=1e-16))
        return res.x.reshape(self.nx, self.nr)

    def action(self, U):
        """S[U] - S[u_bg] (Dirichlet u_bg on all outer boundaries) and its gradient; U has shape (nt, nx, nr)."""
        h, ub = self.h, self.ubg
        Vx = self.Vxr[None]
        Pt = np.concatenate([ub[None], U, ub[None]], axis=0)
        dt = np.diff(Pt, axis=0)
        Px = np.concatenate([np.full((self.nt, 1, self.nr), X_H0), U, np.full((self.nt, 1, self.nr), X_H0)], axis=1)
        dx = np.diff(Px, axis=1)
        Pr = np.concatenate([U, np.full((self.nt, self.nx, 1), X_H0)], axis=2)
        dr = np.diff(Pr, axis=2)
        Vxe = np.concatenate([self.Vxr, self.Vxr[-1:]], axis=0)[None]
        es, gs = self.static(ub.ravel())
        S = (np.sum(0.5 * Vx * dt ** 2) / h + np.sum(0.5 * Vxe * dx ** 2) / h + np.sum(self.Ar[None, None, :] * 0.5 * dr ** 2)
             + h * np.sum(Vx * (p0(U) + 0.5 * self.th[None] * (U - 1) ** 2)) - self.nt * h * es)
        G = h * Vx * (dp0(U) + self.th[None] * (U - 1))
        G += Vx * (dt[:-1] - dt[1:]) / h
        G += (Vxe[:, :-1] * dx[:, :-1] - Vxe[:, 1:] * dx[:, 1:]) / h
        G[:, :, :-1] += -self.Ar[None, None, :-1] * dr[:, :, :-1]
        G[:, :, 1:] += self.Ar[None, None, :-1] * dr[:, :, :-1]
        G[:, :, -1] += -self.Ar[-1] * dr[:, :, -1]
        return S, G

    def C(self, U, w=0.15):
        t = np.tanh((U - X_TOP) / w)
        return self.h * np.sum(self.Vxr[None] * 0.5 * (1 - t)), -self.h * self.Vxr[None] * 0.5 * (1 - t * t) / w


def bounce3(cfg, kappa=1e-3):
    t0 = time.time()
    g = Grid3(**cfg)
    T, X, R = np.meshgrid(g.tau, g.x, g.r, indexing="ij")
    rho = np.sqrt(T ** 2 + X ** 2 + R ** 2)
    xb = np.interp(rho, PROF["r"], PROF["x"], right=X_H0)
    U0 = np.minimum(g.ubg[None], xb)
    c_b = g.C(g.ubg[None] + (xb - X_H0))[0]
    out, Us = [], {}

    def run_c(c, U):
        def F(v):
            Uu = v.reshape(g.nt, g.nx, g.nr)
            S, GS = g.action(Uu)
            Cv, GC = g.C(Uu)
            return S + kappa * (Cv - c) ** 2, (GS + 2 * kappa * (Cv - c) * GC).ravel()
        res = minimize(F, U.ravel(), jac=True, method="L-BFGS-B", options=dict(maxiter=15000, maxcor=30, gtol=1e-7, ftol=1e-15))
        U = res.x.reshape(g.nt, g.nx, g.nr)
        out.append((float(c / c_b), float(g.action(U)[0]), float(g.C(U)[0] / c_b), int(res.nit)))
        return U

    Uc = U0
    for f in np.arange(1.0, 1.61, 0.15):
        Uc = run_c(f * c_b, Uc); Us[round(f, 3)] = Uc
    Uc = U0
    for f in np.arange(0.85, 0.39, -0.15):
        Uc = run_c(f * c_b, Uc); Us[round(f, 3)] = Uc
    kk = int(np.argmax([o[1] for o in out]))
    f0 = out[kk][0]
    start = Us[min(Us, key=lambda k: abs(k - f0))]
    for f in np.linspace(f0 - 0.1, f0 + 0.1, 7):
        run_c(f * c_b, start)
    out.sort()
    Ss = np.array([o[1] for o in out]); cc = np.array([o[2] for o in out])
    k = int(np.argmax(Ss))
    lo, hi = max(k - 2, 0), min(k + 3, len(Ss))
    pf = np.polyfit(cc[lo:hi], Ss[lo:hi], 2)
    b = float(np.polyval(pf, -pf[1] / (2 * pf[0]))) if pf[0] < 0 else float(Ss[k])
    return dict(cfg={k_: v for k_, v in cfg.items()}, b=b, b_grid_max=float(Ss[k]), scan=out, seconds=time.time() - t0,
                ubg_min_on_axis=float(g.ubg[:, 0].min()), ubg_max_on_axis=float(g.ubg[:, 0].max()))


if __name__ == "__main__":
    t0 = time.time()
    q, rs = 1e3, 1.0
    cfgs = [dict(sources="none")]
    for d in (8.0, 12.0, 20.0):
        cfgs += [dict(sources="discrete", d=d, q=q, rs=rs), dict(sources="smeared", d=d, q=q, rs=rs)]
    with Pool(7) as pool:
        res = pool.map(bounce3, cfgs)
    b0 = res[0]["b"]
    for r in res:
        r["db"] = r["b"] - b0
    two_d = NB.bounce(q, rs)                                               # cross-check: one source at the centre, axisymmetric code
    two_d0 = NB.bounce(0.0, rs)
    out = dict(q=q, rs=rs, b0_3d=b0, b0_shooting=NB.B0_SHOOT, runs=res, single_source_2d_db=two_d["b"] - two_d0["b"])
    json.dump(out, open(os.path.join(HERE, "multisource_bounce.json"), "w"), indent=1)
    print("b0 3D %.2f (shooting %.2f)" % (b0, NB.B0_SHOOT))
    for r in res:
        c = r["cfg"]
        print("%-9s d %-5s  b %.2f  db %.2f  axis bg %.3f-%.3f  [%.0fs]" % (c["sources"], c.get("d"), r["b"], r["db"], r["ubg_min_on_axis"],
                                                                       r["ubg_max_on_axis"], r["seconds"]))
    print("single source at the centre (2D code): db %.2f" % out["single_source_2d_db"], "[%.0fs]" % (time.time() - t0))
