# -*- coding: utf-8 -*-
"""Optional neutrino-scalar coupling (a new hypothesis, not part of the dark-matter-only model): what it would imply and what it would not.

Illustrative Dirac benchmark L = -y phi nu_bar nu (flavour-universal; a gauge-invariant completion, e.g. through a higher-dimensional operator with
the Higgs and right-handed neutrinos, is assumed and not specified). The neutrino masses shift by dm = y f x_loc, with x_loc the local minimum of the
held vacuum, which depends on the local dark-matter coupling theta (Section 4.4).

  (1) Terrestrial oscillation experiments (DUNE, Hyper-Kamiokande): along a terrestrial baseline the field sits at the minimum for the local halo
      density, which does not vary along the beam (the scalar couples to dark matter only, not to the Earth), so the shift is a constant addition to
      the mass matrix, absorbed into the measured masses: Delta P_ab(E) = 0 identically. The residual variation of x along a 1300 km baseline from the
      halo density gradient is evaluated.
  (2) Environment and redshift: x_loc differs between voids, the mean, and haloes, and between recombination and today, so the neutrino mass would
      differ by y f Delta x (a mass-varying-neutrino signature).
  (3) Backreaction: relic neutrinos add n_nu m_nu(phi) to the effective potential, a tilt delta eta = 2 n_nu y / (lambda f^3); the trigger is
      unaffected while delta eta << eta.
  (4) After a transition: neutrinos inside the true vacuum have their mass shifted by y f (x_t - x_h).
Writes neutrino_coupling.json.

    python neutrino_coupling.py
"""
import json
import os

import numpy as np

import vacuum_phase_map as VP

HERE = os.path.dirname(os.path.abspath(__file__))
J = lambda f: json.load(open(os.path.join(HERE, f), encoding="utf-8"))
WJ, GR = J("vacuum_wall.json"), J("granularity.json")
ETA = WJ["base"]["eta"]
PS = WJ["physical"]
TH_DEG = ETA / 2.0
F_MEV = PS["f_meV"]
LAM = PS["lam"]
CM3_TO_MEV3 = (1.97327e-5 * 1e3) ** 3                 # 1 cm^-3 in meV^3 (hbar c = 1.97327e-5 eV cm)
N_NU = 336.0                                           # relic neutrinos and antineutrinos per cm^3, all species
F0_MEAN = 3.0                                          # mean coupling today / theta_deg (Section 4.4 baseline: voids at 10% of the mean)
RHO_LOCAL_OVER_MEAN = 0.3 / 1.26e-6                   # local dark matter 0.3 GeV/cm^3 against the cosmic mean


def x_held(theta):
    st = VP.stationary_tb(ETA, theta)
    return max(x for x in st if VP.d2p_tb(x, ETA, theta) > 0)


def x_true(theta):
    st = sorted(x for x in VP.stationary_tb(ETA, theta) if VP.d2p_tb(x, ETA, theta) > 0)
    return st[0] if len(st) == 2 else None


out = {}
th0 = F0_MEAN * TH_DEG
env = {"deepest voids (10% of mean)": 0.1, "mean density": 1.0, "galactic halo (200x mean)": 200.0, "solar neighbourhood": RHO_LOCAL_OVER_MEAN}
xs = {k: x_held(th0 * v) for k, v in env.items()}
ref = xs["solar neighbourhood"]
out["environments"] = {k: dict(one_plus_delta=env[k], theta_over_deg=th0 * env[k] / TH_DEG, x_held=xs[k],
                               dm_nu_over_y_meV=F_MEV * (xs[k] - ref)) for k in env}
zs = (0.0, 0.5, 1.0, 3.0, 1100.0)
out["redshift_mean_density"] = [dict(z=z, x_held=x_held(th0 * (1 + z) ** 3), dm_nu_over_y_meV=F_MEV * (x_held(th0 * (1 + z) ** 3) - x_held(th0)))
                                for z in zs]
# (1) terrestrial: variation of the local minimum across a 1300 km baseline from the halo density gradient (rho ~ r^-2 at 8 kpc)
grad = 2.0 * 1300.0 / (8.0 * 3.086e16)                 # relative density change over 1300 km
th_loc = th0 * RHO_LOCAL_OVER_MEAN
dx = abs(x_held(th_loc * (1 + grad)) - x_held(th_loc))
out["terrestrial"] = dict(theta_local_over_deg=th_loc / TH_DEG, x_local=x_held(th_loc), one_minus_x_local=1.0 - x_held(th_loc),
                          relative_density_change_over_baseline=grad, delta_x_over_baseline=dx, dm_nu_over_y_meV_over_baseline=F_MEV * dx,
                          delta_P="identically zero for a constant shift (absorbed into the measured masses)")
# (3) backreaction on the trigger
n_nu = N_NU * CM3_TO_MEV3
out["backreaction"] = dict(n_nu_meV3=n_nu, delta_eta_per_y=2.0 * n_nu / (LAM * F_MEV ** 3), y_for_delta_eta_1pct_of_eta=0.01 * ETA / (2.0 * n_nu / (LAM * F_MEV ** 3)))
# (4) after a transition
TW = WJ["trigger_wall"]
out["after_transition"] = dict(x_true=TW["x_t"], x_held=TW["x_f"], dm_nu_over_y_meV=F_MEV * (TW["x_t"] - TW["x_f"]))
json.dump(out, open(os.path.join(HERE, "neutrino_coupling.json"), "w"), indent=1)
for k, v in out.items():
    print(k, v)
