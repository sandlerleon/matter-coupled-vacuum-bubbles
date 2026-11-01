# -*- coding: utf-8 -*-
"""Paper 2 (Physics of the Dark Universe): density-dependent vacuum bubble nucleation in matter-coupled dark-sector cosmology.

    python paper2_build.py     -> manuscript/PDU_Vacuum_Bubbles.{docx,tex}, supplement, cover letter
                                  (after interface_validation.py, granularity.py, make_int_figures.py, make_gran_figure.py)
"""
import datetime
import json
import os
import platform
import re
import sys

import numpy as np
import scipy
from docx import Document
from docx.shared import Pt, Inches

import copyedit as CE
import int_blocks as IB
import render as RN

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, ".."))
OUT = os.path.join(REPO, "manuscript")
FIGDIR = os.path.join(REPO, "figures")
J = lambda f: json.load(open(os.path.join(HERE, f), encoding="utf-8"))
if not IB.VJ["all_passed"]:
    sys.exit("interface checks have not all passed; run interface_validation.py")
DOIS = J("dois.json")
COMP = J("companion.json") if os.path.exists(os.path.join(HERE, "companion.json")) else {}
GITHUB = "https://github.com/sandlerleon/matter-coupled-vacuum-bubbles"
REFS = dict(J("refs_all.json"))
REFS.update({
    "Companion": "L. Sandler, Matter-dilution-triggered recollapse and conditional nonsingular cycles in a metastable dark-energy model. Preprint, "
                 "Zenodo (2026). https://doi.org/%s" % COMP.get("paper", "(to be assigned)"),
    "HOTW17": "L. Hui, J.P. Ostriker, S. Tremaine, E. Witten, Ultralight scalars as cosmological dark matter. Phys. Rev. D 95, 043541 (2017). "
              "https://doi.org/10.1103/PhysRevD.95.043541",
    "PWW83": "J. Preskill, M.B. Wise, F. Wilczek, Cosmology of the invisible axion. Phys. Lett. B 120, 127--132 (1983). "
             "https://doi.org/10.1016/0370-2693(83)90637-8",
    "TG79": "S. Tremaine, J.E. Gunn, Dynamical role of light neutral leptons in cosmology. Phys. Rev. Lett. 42, 407--410 (1979). "
            "https://doi.org/10.1103/PhysRevLett.42.407",
    "GMW14": "R. Gregory, I.G. Moss, B. Withers, Black holes as bubble nucleation sites. J. High Energy Phys. 03, 081 (2014). "
             "https://doi.org/10.1007/JHEP03(2014)081",
    "BGM15": "P. Burda, R. Gregory, I.G. Moss, Vacuum metastability with black holes. J. High Energy Phys. 08, 114 (2015). "
             "https://doi.org/10.1007/JHEP08(2015)114",
    "FNW04": "R. Fardon, A.E. Nelson, N. Weiner, Dark energy from mass varying neutrinos. J. Cosmol. Astropart. Phys. 10, 005 (2004). "
             "https://doi.org/10.1088/1475-7516/2004/10/005",
    "DUNE20": "DUNE Collaboration, B. Abi et al., Deep Underground Neutrino Experiment (DUNE), Far Detector Technical Design Report, Volume II: DUNE "
              "Physics. arXiv:2002.03005 (2020).",
    "HK18": "Hyper-Kamiokande Proto-Collaboration, K. Abe et al., Hyper-Kamiokande Design Report. arXiv:1805.04163 (2018).",
    "CRS14": "F.-Y. Cyr-Racine, K. Sigurdson, Limits on neutrino-neutrino scattering in the early Universe. Phys. Rev. D 90, 123533 (2014). "
             "https://doi.org/10.1103/PhysRevD.90.123533",
    "AH14": "M. Archidiacono, S. Hannestad, Updated constraints on non-standard neutrino interactions from Planck. J. Cosmol. Astropart. Phys. 07, 046 "
            "(2014). https://doi.org/10.1088/1475-7516/2014/07/046",
    "SvdW04": "R.K. Sheth, R. van de Weygaert, A hierarchy of voids: much ado about nothing. Mon. Not. R. Astron. Soc. 350, 517--538 (2004). "
              "https://doi.org/10.1111/j.1365-2966.2004.07661.x",
})
SV, DBJ, GRJ = J("survival.json"), J("density_bounce.json"), J("granularity.json")
TITLE = ("Density-Dependent Vacuum Bubble Nucleation in Matter-Coupled Dark-Sector Cosmology: Dark-Matter Granularity, Wall Dynamics and "
         "Transmission")
GF, WWJ, RSJ, VEJ, NUJ = J("granular_field.json"), J("wave_wall_1d.json"), J("region_sensitivity.json"), J("viability_extended.json"), J("neutrino_coupling.json")
NBJ, SRJ, EEV, MSJ = J("nonpert_bounce.json"), J("screened_rate.json"), J("env_evolution.json"), J("multisource_bounce.json")
NBC = {(c["h"], c["L"]): c for c in NBJ["convergence"]}
MSR = {(r["cfg"]["sources"], r["cfg"].get("d")): r for r in MSJ["runs"]}
NBQ = {(b["q"], b["rs"]): b for b in NBJ["bounces"]}
SR1 = min(SRJ["rows"], key=lambda r: abs(np.log(r["m_dm_meV"] / 1e9)))
CMB = NUJ["cmb"]
CONV = WWJ["convergence_table"]
CONV_F = [c for c in CONV if c["L_box"] == 200.0]
texsci, sci = IB.texsci, IB.sci
PM, WJ, TJ, GJ, CJ, EJ, VJ, HJ, VM = IB.PM, IB.WJ, IB.TJ, IB.GJ, IB.CJ, IB.EJ, IB.VJ, IB.HJ, IB.VM
COH, PS, TW, DM, LNG, NOREG, DW = IB.COH, IB.PS, IB.TW, IB.DM, IB.LNG, IB.NOREG, IB.DW
PC, WB, IS, RP = HJ["percolation"], HJ["wall_bounds"], HJ["isolated"], HJ["reproduction"]
VAR = {v["variant"]: v for v in HJ["variants"]}
W0 = VAR["wall energy as massive scalar quanta until the bounce (w = 0)"]
W13 = VAR["wall energy as dark radiation until the bounce (w = 1/3)"]
TM = VM["tunneling"]
B_P, B_I = IB.B_P, IB.B_I
ev = IB.ev
TB = WJ["walls"]["0.30"]["rows"]
r0 = TB[0]
r6 = [r for r in TB if abs(r["theta_over_deg"] - 0.6) < 1e-9][0]
deg = WJ["walls"]["0.30"]["degenerate"]
VC = IB.VC
GT, G0 = GRJ["trigger"], GRJ["today"]
VOID = {round(v["one_plus_delta"], 2): v for v in GRJ["voids"]}
DV_ = 0.1                                       # density of the emptiest regions relative to the mean (deepest voids)
F0_LOC = WJ["base"]["f0"]                        # the baseline theta_0/theta_deg, read as the local value in those regions
F0_MEAN = F0_LOC / DV_
THR_T = 1.0 / (F0_MEAN * WJ["a_trigger"] ** -3)  # 1 + delta above which the held vacuum is still the true vacuum at the trigger
NEARC = GRJ["near_eta_c_max"]
NEAR01 = max((r for r in GRJ["near_eta_c"] if r["theta_over_deg"] >= 0.1), key=lambda r: r["ratio"])
VCOMP = GRJ["void_completion"]
CELL = {c["Rc"]: c for c in GF["cells"]}
PEX = GF["particle_examples"]
WWR = {(r["mode"], r["m_ratio"], r["dz"]): r for r in WWJ["runs"] if r.get("L_box", 200.0) == 200.0}
W10, W3 = WWR[("wave", 10.0, 0.005)], WWR[("wave", 3.0, 0.005)]
SHELL10 = float(np.median([q["shell_over_before"] for q in W10["probes"]]))
VEB, VEF = VEJ["surveys"]
PRED = {(r["mean_f0"], r["one_plus_delta"]): r for r in VEJ["prediction"]}
PREF = max(max(abs(r["prefactor_up"] / r["central"] - 1), abs(r["prefactor_down"] / r["central"] - 1)) for r in VEJ["prediction"] if r["central"] > 1.0)
RSI, RSR = RSJ["isolated"], RSJ["regions"]
NUE = NUJ["environments"]
HIGHLIGHTS = ["Dilution makes a matter-held vacuum metastable below a degeneracy density",
              "Discrete sources suppress the decay far less than the mean field predicts",
              "Wave dark matter between 3 meV and 0.25 eV restores the smooth mean field",
              "Regions decay in order of their density; galactic halos keep the old vacuum",
              "Dark matter crossing the runaway wall is kicked forward at three quarters of c"]

B = []
add = B.append
ABSTRACT = (
    "Can a declining dark-matter density switch on the decay of a vacuum? For a scalar whose coupling to dark matter is centered on a false vacuum, the "
    "mean-field effective potential becomes metastable below a degeneracy density, and Coleman--De Luccia tunneling switches on. We test the assumption "
    "behind this picture, that the field feels a smooth density. For a vacuum at the dark-energy scale the field's Compton length is about %.2f mm, and "
    "the dark-matter energy inside a Compton volume is only %.2f$m_\\varphi^4$. A spatially resolved treatment, the static field around discrete sources "
    "and the decay rate with randomly placed sources at first order in each source's coupling, together with nonperturbative bounces that deform "
    "around a single source and around a row of sources, indicates that discrete dark-matter sources suppress the decay far less than a homogeneous "
    "mean-field treatment predicts: for 1 MeV particles, a fraction $%s$ of the mean-field suppression survives at first order, and less once the "
    "particles screen themselves. Within this model, particle dark matter therefore leaves the vacuum close to its zero-density rate; general "
    "many-body source configurations remain open. The mean field holds for wave "
    "dark matter heavier than the field and still classical, about %.0f meV to %.2f eV, which a 1+1-dimensional simulation of a wall running through "
    "such dark matter confirms. The decay is then environmental: a region completes its transition at a scale factor proportional to $(1+\\delta)^{1/3}$, "
    "voids first, while galactic halos keep the held vacuum. The walls run away; the dark matter crossing them gains the sudden-crossing energy and trails "
    "the wall at %.2fc. Transitioned regions recollapse, and within the thin-wall junction model, an isolated Hubble-scale bubble develops a trapped "
    "surface before its interior bounces. An optional neutrino coupling is limited by an approximate early-universe free-streaming criterion to "
    "$y\\lesssim10^{-4}$, at which it has no observable consequence."
    % (GT["compton_length_mm"], GT["rho_over_mphi4"], texsci(PEX["1 MeV"]["fraction_of_first_order_mean_field_suppression"], 0), GT["m_phi_meV"],
       GT["classical_max_mass_meV"] / 1e3, W10["v_dm_behind_late"]))
add(("title", TITLE))
add(("author", "Leon Sandler\nIndependent Researcher, Northbrook, Illinois, USA\nsandler.leon@gmail.com"))
add(("abstract", ABSTRACT))
add(("keywords", "vacuum decay; bubble nucleation; coupled dark matter; wave dark matter; dark energy; cosmic voids; mass-varying neutrinos"))

# ============================================================================================ 1
add(("h1", "1. Introduction"))
add(("p", "False-vacuum decay proceeds by the nucleation of bubbles of the true vacuum, whose rate follows from the Euclidean bounce \\cite{Coleman,CC} and "
          "whose gravitational corrections were derived by Coleman and De Luccia \\cite{CDL}. The theory is usually applied to a vacuum whose stability does not "
          "depend on its surroundings. In a universe filled with dark matter, however, a scalar that couples to the dark matter feels an effective potential "
          "that depends on the dark-matter density, as in chameleon and symmetron screening \\cite{KW,HK,HKLM11} and in coupled dark energy \\cite{Amendola}. "
          "If the dark energy is the energy of such a scalar held in a vacuum by the matter term, the dilution of the dark matter by the expansion changes "
          "the vacuum structure itself, and the decay of the vacuum is placed under the control of the cosmic density. This paper asks three questions. Can a "
          "declining dark-matter density turn a stable vacuum into a metastable one? When is the density that the field feels the smooth mean density "
          "that such an argument assumes? What do the resulting vacuum walls do: what crosses them, how do they move, and what becomes of the region they "
          "enclose?"))
add(("p", "Environment-dependent vacua are familiar from screened modified gravity, where the field's minimum and mass track the local density "
          "\\cite{KW,HK,HKLM11}; in those applications, the field's Compton length is long compared with the separation of the matter particles, and the "
          "density can be averaged. Vacuum decay out of a metastable state has been studied for vacua whose stability is fixed, from the Standard Model Higgs "
          "\\cite{IRS,EGR} to landscapes of many vacua \\cite{BP,FKMS}, including the signatures of bubble collisions \\cite{Kleban,AJ}. Combining the two "
          "raises the question that neither literature had to face: a vacuum at the dark-energy scale has a Compton length below a millimeter, and whether the "
          "dark matter looks smooth on that scale depends on what the dark matter is. The interfaces add a further question. A wall between vacua in which the "
          "dark-matter mass differs is a mass step; slow walls reflect particles whose momentum is below the step, which has been proposed as a way to filter "
          "dark matter at first-order transitions \\cite{BKL20}, while fast walls transmit them and the pressure they exert is bounded \\cite{BM09,BM17}. "
          "Bubbles of a negative vacuum enclose regions that recollapse, and bubbles in an expanding exterior can be cut off behind horizons "
          "\\cite{BGG,FG,FGG}."))
add(("p", "The contributions of this study are as follows. (i) The phase structure of a matter-coupled potential is mapped as a function of density and tilt "
          "(Section 2). (ii) The conditions under which the field feels a smooth density are derived and tested with spatially resolved calculations; within "
          "them, particle dark matter does not provide a smooth source, and a window remains for wave dark matter heavier than the field, a bound that "
          "follows from the dark-energy scale itself (Section 3). (iii) The bounce, the "
          "wall and the density-dependent action are computed and validated, and the nucleation history is followed in the emptiest regions, where it begins, "
          "with the consequences of the environment for the protection of the past and for the regions that keep the held vacuum (Section 4). (iv) The "
          "transmission of fields and of the dark matter through the wall is computed (Section 5), (v) the gravitational dynamics of the wall and the "
          "completion of the transition are computed with explicit bounds on the approximations (Section 6), and (vi) the region enclosed by the wall is "
          "evolved as a cosmology (Section 7). (vii) Controls that remove the coupling, a viability map and parameter surveys identify where the conditions "
          "hold (Section 8), and (viii) an optional coupling of the scalar to neutrinos is examined for its experimental consequences (Section 9)."))
add(("p", "Units are $\\hbar=c=k_B=1$, $M_{\\rm Pl}$ is the reduced Planck mass, and every number quoted is read from the output of the script named in "
          "parentheses in the archived code. The paper does not require a cosmological bounce. A companion paper treats the symmetron version of the trigger, "
          "whose release is a classical roll on cosmological scales, and its continuation into cycles \\cite{Companion}; the two papers share the class of "
          "models but not their results."))

# ============================================================================================ 2
add(("h1", "2. Matter-coupled vacuum metastability"))
add(("p", "A scalar coupled conformally to dark matter has the action"))
add(("eq", r"S=\int d^4x\sqrt{-g}\left[\frac{M_{\rm Pl}^2R}{2}-\frac{1}{2}(\nabla\varphi)^2-V_0(\varphi)\right]+S_{\rm dm}\left[A^2(\varphi)g_{\mu\nu},\psi_{\rm dm}\right]"))
add(("p", "and, where the dark matter can be described by a smooth nonrelativistic density, the conserved Einstein-frame density $\\rho_*$ enters the field "
          "equation through $V_{\\rm eff}(\\varphi,\\rho_*)=V_0(\\varphi)+\\rho_*A(\\varphi)$ \\cite{Amendola,HK}. When this mean-field description holds is "
          "the subject of Section 3; here it is assumed. We take the tilted quartic $V_0$ in the dimensionless form"))
add(("eq", r"p(x;\eta,\theta)=\frac{(x^2-1)^2}{8}+\frac{\eta x}{2}+\frac{\theta}{2}(x-1)^2,\qquad x=\frac{\varphi}{f},\qquad \theta=\frac{\rho_*}{\lambda f^2M^2},\qquad V_{\rm eff}=\lambda f^4\,p+{\rm const}"))
add(("p", "where $\\eta$ tilts the potential so that $x\\approx-1$ is the zero-density true vacuum and $x\\approx1$ the false one, and the coupling "
          "$A(\\varphi)=1+(\\varphi-f)^2/2M^2$ is centered on the false vacuum, so that the matter term $\\theta(x-1)^2/2$ vanishes there and deepens it. "
          "Because $\\rho_*\\propto a^{-3}$, the mean $\\theta=\\theta_0a^{-3}$ decreases monotonically as the universe expands (vacuum_phase_map.py)."))
add(("p", "The stationary points solve $x^3+(2\\theta-1)x+(\\eta-2\\theta)=0$. Two minima coexist for $\\theta_{\\rm low}<\\theta<\\theta_{\\rm high}$, "
          "with $\\theta_{\\rm low}=0$ for $\\eta<\\eta_c=2/(3\\sqrt3)$, and the window edges are roots of the cubic discriminant (Fig. 1b). The matter-held "
          "minimum $x_h\\approx1$ is the true vacuum above the degeneracy density, which is exactly"))
add(("eq", r"\theta_{\rm deg}=\frac{\eta}{2},\qquad p(x;\eta,\eta/2)=\frac{(x^2-1)^2}{8}+\frac{\eta}{4}\left(x^2+1\right)"))
add(("p", "because $p$ is then even in $x$. Below $\\theta_{\\rm deg}$ the held vacuum is metastable: dilution turns a stable vacuum into a decaying one, "
          "and tunneling is meaningful for $0<\\theta<\\eta/2$. Above $\\theta_{\\rm high}$ (%.3f for $\\eta=0.30$) the other minimum disappears; thus, at high "
          "density, the matter term alone determines the vacuum. For $\\eta>\\eta_c$ the held minimum itself disappears below $\\theta_{\\rm low}$ (for example, "
          "%.3f at $\\eta=0.45$), and the release is a classical roll rather than a tunneling event (Fig. 1a)."
          % (PM["tunneling_profiles"]["0.30"]["window"][1], PM["tunneling_profiles"]["0.45"]["window"][0])))
add(("p", "The same structure is absent in the symmetron form of the coupling, $A=1+\\varphi^2/2M^2$ with a symmetric double well: there the held point "
          "$\\varphi=0$ loses stability exactly when the broken minima appear, the occupied minimum is never metastable, and Coleman--De Luccia tunneling "
          "does not apply. Its interfaces are domain walls between $\\pm v$, which form if the roll is not seeded coherently; they would carry up to %.0f%% of "
          "the critical density per Hubble volume, and a tilt $\\epsilon\\gtrsim%s\\,\\mu^2v$ suffices to prevent them. The coupling energy of the dark matter "
          "between the vacua is $\\rho_*\\Delta A$ exactly; the first-order estimate $\\delta(\\Delta V)=\\rho_m\\Delta\\ln m$ agrees with it to %.2f%% when the "
          "coupling changes the rate by a factor $e$ and to %.0f%% when it moves the action across the survival window. The tunneling branch is the one in "
          "which density-dependent bubble nucleation can occur, and the rest of the paper focuses on this branch."
          % (100 * COH["domain_walls"][-1]["wall_energy_per_Hubble_volume_over_rho_crit"], texsci(COH["eps_hat_min"], 0),
             100 * (min(c["exact_over_log"] for c in PM["coupling_check"]) - 1), 100 * (max(c["exact_over_log"] for c in PM["coupling_check"]) - 1))))
add(("fig", "fig1_vacuum_landscape.png", "Fig. 1. (a) Effective potential at five densities ($\\eta=0.30$): dilution turns the matter-held minimum near "
            "$x=1$ from the true into the false vacuum at $\\theta_{\\rm deg}=\\eta/2$. (b) Phase diagram: number of minima and which is lower, with the "
            "degeneracy line and the lines where a minimum disappears; the gray line indicates $\\eta_c$.", 6.3))

# ============================================================================================ 3 granularity
add(("h1", "3. When does the field feel a smooth dark-matter density?"))
add(("h2", "3.1 The scales of the problem"))
add(("p", "The effective potential $V_0+\\rho_*A$ is a mean-field statement. The field responds to sources within its Compton length "
          "$m_\\varphi^{-1}$, with $m_\\varphi^2=\\lambda f^2\\,\\partial_x^2p(x_h)$ at the held vacuum, and on the time scale $m_\\varphi^{-1}$; the bounce "
          "itself has a radius of a few Compton lengths. Normalized to the observed dark energy (Section 4.3), $m_\\varphi=%.2f$ meV at the trigger and %.2f meV "
          "today, a Compton length of %.3f mm (granularity.py). The mean dark-matter energy density at the trigger is %.2f meV$^4$, so"
          % (GT["m_phi_meV"], G0["m_phi_meV"], GT["compton_length_mm"], GT["rho_dm_meV4"])))
add(("eq", r"\frac{\rho_{\rm dm}}{m_\varphi^4}=%.3f\ \ (\text{trigger}),\qquad %.3f\ \ (\text{today})" % (GT["rho_over_mphi4"], G0["rho_over_mphi4"])))
add(("p", "Less than one field-mass worth of dark matter lies in a Compton volume. Whether the source is nevertheless smooth depends on what the dark "
          "matter is, and three cases must be distinguished (Table 1, Fig. 2a)."))
add(("h2", "3.2 Particle dark matter"))
sp = GRJ["spacing_over_compton"]
add(("p", "If the dark matter consists of particles whose occupation of phase space is small, the source of the field is a set of points. The mean field "
          "applies only if many of them lie inside a Compton volume, $n_{\\rm dm}\\gg m_\\varphi^3$, which requires $m_{\\rm dm}\\ll\\rho_{\\rm dm}/m_\\varphi^3"
          "=%.2f$ meV. A particle description, on the other hand, needs an occupation per de Broglie cell below one, $m_{\\rm dm}\\gg"
          "(\\rho_{\\rm dm}/v^3)^{1/4}=%.0f$ meV for a velocity dispersion $v=10^{-3}$. The two conditions exclude each other, and fermions lighter than a few "
          "hundred eV are excluded in any case by phase-space bounds \\cite{TG79}. For thermal relics, the mismatch is large: the mean separation of 1 MeV "
          "particles is %.0f times the Compton length and that of 1 GeV particles %.0f times. Each particle is a screened point source around which the field "
          "is pinned within a Compton length, as for the thin-shell objects of chameleon screening \\cite{KW}; between them, the field sees zero density. If so, "
          "the held vacuum decays in the space between the particles at the zero-density rate, $\\ln(\\Gamma/H^4)=%.0f$ in the baseline, and the expected "
          "number of bubbles in our past light cone is $e^{%.0f}$. Under this local mean-field criterion, particle dark matter does not provide a smooth "
          "source on the response scale of the field; Section 3.5 tests the consequence with a spatially resolved calculation. The bound is least severe where it does not matter: in halos the density is highest, but there the "
          "held vacuum is the true vacuum and nothing decays (Section 4.4); in the voids, where the decay must begin, the density is lower than the mean and "
          "the bound is stronger."
          % (GT["particle_max_mass_meV"], GT["classical_max_mass_meV"], sp["1 MeV"]["over_compton"], sp["1 GeV"]["over_compton"],
             GRJ["particle_dm"]["ln_gamma_over_H4_zero_density"], GRJ["particle_dm"]["ln_expected_bubbles_past_light_cone"])))
add(("h2", "3.3 Wave dark matter"))
add(("p", "Bosonic dark matter with a large occupation number is a classical field oscillating at its mass \\cite{HOTW17}, and its coupling through "
          "$A^2(\\varphi)g_{\\mu\\nu}$ gives the field a source proportional to the instantaneous trace of its stress-energy. Two regimes follow. If the dark "
          "matter is lighter than the field, $m_{\\rm dm}\\ll m_\\varphi$, the field follows the oscillation; the source passes through zero twice per period, "
          "and tunneling, which takes a time of order $m_\\varphi^{-1}$, proceeds in those intervals at the zero-density rate, with the same consequence as for "
          "particles. If it is heavier, $m_{\\rm dm}\\gg m_\\varphi$, the source oscillates faster than the field can respond, and the field feels its time "
          "average, which for a nonrelativistic field is $\\rho_*A(\\varphi)$. The source is also smooth in space, because the field is coherent over its de "
          "Broglie length $1/(m_{\\rm dm}v)$, which exceeds the Compton length for $m_{\\rm dm}\\ll m_\\varphi/v=%.1f$ eV. Together with classicality, the mean "
          "field of Section 2 applies for"
          % (GT["coherence_max_mass_meV"] / 1e3)))
add(("eq", r"m_\varphi=%.1f\ {\rm meV}\ \ll\ m_{\rm dm}\ \ll\ \left(\frac{\rho_{\rm dm}}{v^3}\right)^{1/4}=%.0f\ {\rm meV}" % (GRJ["wave_window_meV"][0], GRJ["wave_window_meV"][1])))
add(("p", "at the trigger, a window of about two decades. Its edges have definite physical origins. The lower edge is the mass of the field, a "
          "property of the vacuum fixed by the observed dark energy and by the requirement $B\\approx%.0f$ (Section 3.4); below it the field resolves the "
          "oscillation of the dark matter. The upper edge is the classicality of the dark matter, where its occupation per de Broglie cell falls to one; it "
          "scales as $(\\rho_{\\rm dm}/v^3)^{1/4}$ and is therefore higher in the colder, emptier voids where the decay begins (by a factor of about 3 for a "
          "tenth of the mean density and $v=10^{-4}$). The coherence condition $m_{\\rm dm}\\ll m_\\varphi/v$ is weaker. Dark matter in this mass range is a "
          "nonthermal condensate produced, for instance, by misalignment \\cite{PWW83}, which behaves as cold dark matter on all scales that have been "
          "observed \\cite{HOTW17} and does not need to couple to the Standard Model. The trigger therefore specifies the dark matter, a classical bosonic field of meV "
          "to sub-eV mass, and everything that follows assumes it." % TW["B"]))
add(("p", "The averaging can be made explicit. We take the dark matter to be a real scalar $\\chi$ whose mass depends on the field, $m_{\\rm dm}A(\\varphi)$, "
          "the coupling used in the simulation of Section 5. Its source in the field equation is $m_{\\rm dm}^2AA'\\chi^2$, which for a homogeneous "
          "nonrelativistic oscillation is $\\rho_*A'[1+\\cos(2m_{\\rm dm}t)]$: never negative, with the time average of a pressureless fluid. The response of a "
          "field of mass $m_\\varphi$ to the oscillating part is smaller than its response to a static source by $m_\\varphi^2/(4m_{\\rm dm}^2-m_\\varphi^2)$, so "
          "the residual modulation of the density the field feels is $(m_\\varphi/2m_{\\rm dm})^2$: 1% for $m_{\\rm dm}=5m_\\varphi$, the lower edge used in "
          "Section 8, growing to order one as $m_{\\rm dm}$ approaches $m_\\varphi$, where the two regimes meet. A coupling through the metric, "
          "$A^2(\\varphi)g_{\\mu\\nu}$, would give instead the source $-(A'/A)T_{\\rm dm}$ with $T_{\\rm dm}=-\\rho[1+3\\cos(2m_{\\rm dm}t)]$, which changes sign "
          "during each period; at early times, when the field is far heavier than the dark matter (Section 9), the field would follow that sign change, so "
          "the mass coupling is the one assumed. Fluctuations in the dark-matter density on scales between the de Broglie length and the size of a void are "
          "part of the local density of Section 4.4, not of its granularity. Section 5 tests the averaging dynamically, with a wall running through wave "
          "dark matter."))
rowsT1 = [["particle, any mass", "point sources; $n\\gg m_\\varphi^3$ needs $m<%.2f$ meV, particle regime needs $m>%.0f$ meV" % (GT["particle_max_mass_meV"],
                                                                                                                         GT["classical_max_mass_meV"]),
           "zero-density rate between particles", "excluded"],
          ["fermion, light", "phase-space (Tremaine--Gunn) bound", "--", "excluded"],
          ["wave, $m_{\\rm dm}\\ll m_\\varphi$", "field follows the oscillation", "zero-density rate at the nodes", "excluded"],
          ["wave, $m_\\varphi\\ll m_{\\rm dm}\\ll%.0f$ meV" % GT["classical_max_mass_meV"], "time average, coherent over the de Broglie length",
           "mean field $\\rho_*A(\\varphi)$, local density", "viable"]]
add(("table", "Table 1. Coupled dark matter and the density the field feels at the trigger ($m_\\varphi=%.1f$ meV, $\\rho_{\\rm dm}=%.2f$ meV$^4$, $v=10^{-3}$)."
              % (GT["m_phi_meV"], GT["rho_dm_meV4"]), ["Dark matter", "Condition", "What the field feels", "Trigger"], rowsT1, None))
add(("h2", "3.4 Why the bound is generic"))
add(("p", "The smallness of $\\rho_{\\rm dm}/m_\\varphi^4$ is not a property of the baseline. With $\\lambda f^4=(1+q)V_f/\\Delta\\tilde V_0$ fixed by the "
          "observed dark energy and $m_\\varphi^4=\\lambda^2f^4\\,(\\partial_x^2p)^2$,"))
add(("eq", r"\frac{\rho_{\rm dm}}{m_\varphi^4}=\frac{\rho_{\rm dm}}{V_f}\,\frac{\Delta\tilde V_0}{1+q}\,\frac{1}{\lambda\,(\partial_x^2p)^2},\qquad \lambda=\frac{b(\eta,\theta)}{B}"))
add(("p", "The first factor is of order one because the dark matter and the dark energy have comparable densities, which is what makes a dilution "
          "trigger possible in the first place; and $\\lambda$ is fixed by the requirement that the transition happen within a Hubble time, $B\\approx%.0f$. "
          "Over all %d combinations of tilt, density and $q$ in the bounce tables, the ratio remains below %.1f. Toward $\\eta_c$, where the held minimum becomes "
          "shallow, both $\\lambda$ and $\\partial_x^2p$ fall and the ratio increases, but only at vanishing density: it reaches %.0f at $\\eta=%.3f$ and "
          "$\\theta=0$, and stays below %.1f once $\\theta\\ge0.1\\,\\theta_{\\rm deg}$, the densities at which the decay is triggered. Particle dark matter "
          "would need a ratio above $1/v\\approx10^3$. The coincidence of the dark-matter and dark-energy densities, which motivates a dilution trigger, is "
          "thus what makes particle dark matter fail the criterion."
          % (TW["B"], GRJ["structural_scan_n"], GRJ["structural_scan_max"]["ratio"], NEARC["ratio"], NEARC["eta"], NEAR01["ratio"])))
add(("fig", "fig2_granularity.png", "Fig. 2. (a) Conditions for the field to feel a smooth dark-matter density at the trigger: particle dark matter "
            "(red) needs two incompatible conditions; wave dark matter (blue) is averaged where it is heavier than the field and still classical (green "
            "band). (b) Expected number of bubbles in the past light cone of a region whose density is a fraction $1+\\delta$ of the mean, for three mean "
            "couplings today ($\\eta=0.30$, $q=1$, $B_0=%.0f$): the emptiest regions decide whether the past is protected." % WJ["base"]["B0"], 6.6))
add(("h2", "3.5 Spatially resolved tests"))
add(("p", "Two calculations replace the criterion with a resolved treatment (granular_field.py). The first solves the static field equation "
          "$\\nabla^2x=\\partial_xp_0(x)+\\theta_s(\\mathbf r)(x-1)$ in a spherical Wigner--Seitz cell around one source, with the source spread over a "
          "sphere small compared with the cell and normalized so that the cell average of $\\theta_s$ is the mean coupling at the trigger. For cells smaller "
          "than the Compton length the field is uniform at the mean-field minimum; for a cell radius $R_c=%.1f/m_\\varphi$ the field lies within 1%% of the "
          "zero-density vacuum in %.0f%% of the cell volume, and in %.1f%% for $R_c=%.0f/m_\\varphi$ (Fig. 10a); the fraction changes by less than 1%% "
          "when the radius of the source is varied tenfold. The crossover lies at a cell radius of a few Compton lengths."
          % (CELL[10.0]["Rc_times_mphi"], 100 * CELL[10.0]["fraction_at_zero_density_vacuum"], 100 * CELL[30.0]["fraction_at_zero_density_vacuum"],
             CELL[30.0]["Rc_times_mphi"])))
add(("p", "The second computes the decay rate in a granular background. A static source at distance $r_p$ from the centre of the zero-density bounce "
          "adds to the Euclidean action, at first order in its coupling $q$ (the coupling carried by one particle), $\\Delta B=(q/\\lambda)I(r_p)$ with "
          "$I(r_p)=\\int d\\tau\\,\\frac12[(x_b-1)^2-(x_h-1)^2]$ along its world line. For sources placed at random with number density $n$,"))
add(("eq", r"\frac{\Gamma}{\Gamma_0}=\Big\langle e^{-\sum_i\Delta B_i}\Big\rangle=\exp\!\Big[-n\int d^3r\,\big(1-e^{-qI(r)/\lambda}\big)\Big]"))
add(("p", "For many weak sources, this reduces to the first-order mean-field suppression $\\exp(-\\theta\\,\\partial_\\theta b/\\lambda)$, and the "
          "derivative obtained in this way agrees with the numerical derivative of the bounce action to %.2f%%. For few strong sources, it reduces to "
          "$\\Gamma_0$ times the probability that the bounce region contains no source. At the trigger density, the suppression decreases to half of its "
          "mean-field value for point sources of mass %.0f meV, and to fractions $%s$ and $%s$ of it for 1 MeV and 1 GeV particles (Fig. 10b). These are "
          "first-order results; the particles are treated as static over the time of the bounce, which holds for $v\\ll1$."
          % (100 * abs(GF["bounce"]["first_order_vs_numerical"]), GF["m_half_meV"],
             texsci(PEX["1 MeV"]["fraction_of_first_order_mean_field_suppression"], 0),
             texsci(PEX["1 GeV"]["fraction_of_first_order_mean_field_suppression"], 0))))
add(("p", "For heavy particles, the coupling of a single source is large, and the first-order treatment must be tested beyond first order. This was "
          "done with a nonperturbative bounce (nonpert_bounce.py). With a source on the axis of the bubble, the worst case, the field depends only on "
          "Euclidean time and on the distance from the source, and the bounce is found as the maximum, over a collective coordinate (the four-volume on the "
          "true-vacuum side of the barrier), of the constrained minimum of the action, so that the bubble is free to deform around the source. Without a "
          "source, the method reproduces the shooting action to %.1f%% (%.2f%% on the finest grid). The numbers below are on a grid with a spacing 0.4. The "
          "convergence tests (Table S12b) show that enlarging the box does not change the value and that the saturation in $q$ is present at every resolution: "
          "between $q=10^3$ and $10^4$ the increase changes by %.1f%% on a grid with a spacing 0.2 and by %.1f%% at 0.25. The saturation value itself depends "
          "on how the source, one to three grid cells across, is represented: %.0f--%.0f for grid spacings 0.2--0.4 at $q=10^4$, and %.0f--%.0f at $q=10$. "
          "It does not converge and is quoted only to show that a single source increases the action by a bounded amount—at most about %.0f%% of the "
          "zero-density value and three to four orders of magnitude below the first-order estimate. The nonperturbative increase in the action is below the first-order estimate for "
          "every coupling tested: their ratio is %.2f, %.2f, %.3f and %.4f for $q=1$, 10, 100 and $10^4$, because the field near the source saturates and "
          "the bubble deforms around it (Fig. 10c). The increase saturates at a value set by the size of the source, $\\Delta b=%.0f$, %.0f and %.0f for "
          "source radii 0.25, 0.5 and 1, and the static coupling of a source saturates in the same way, at $q_{\\rm eff}\\approx4\\pi r_s$ (within %.0f%% for "
          "$r_s\\le0.1$; every case converges at two radial resolutions, Table S13), the thin-shell screening known from chameleon fields \\cite{KW}. A nonrelativistic particle is spread at most over its de Broglie "
          "length $1/(m_{\\rm dm}v)$, far below the Compton length of the field for particle dark matter, so its effective coupling is smaller than $q$ by a "
          "factor $%s$ for 1 MeV particles, and the granular suppression of the rate falls accordingly (screened_rate.py)."
          % (100 * abs(NBJ["grid_error"]), 100 * abs(NBJ["saturation_fine"]["h02"]["b0"] / NBJ["b0_shooting"] - 1),
             100 * abs(NBJ["saturation_fine"]["h02"]["db_q1e4"] / NBJ["saturation_fine"]["h02"]["db_q1e3"] - 1),
             100 * abs(NBJ["saturation_fine"]["h025"]["db_q1e4"] / NBJ["saturation_fine"]["h025"]["db_q1e3"] - 1),
             min([c["db_q1e4"] for c in NBJ["convergence"]] + [NBJ["saturation_fine"]["h02"]["db_q1e4"]]),
             max([c["db_q1e4"] for c in NBJ["convergence"]] + [NBJ["saturation_fine"]["h02"]["db_q1e4"]]),
             min(c["db_q10"] for c in NBJ["convergence"]), max(c["db_q10"] for c in NBJ["convergence"]),
             100 * max([c["db_q1e4"] for c in NBJ["convergence"]] + [NBJ["saturation_fine"]["h02"]["db_q1e4"]]) / NBJ["b0_shooting"],
             NBQ[(1.0, 0.5)]["ratio_np_over_first"], NBQ[(10.0, 0.5)]["ratio_np_over_first"],
             NBQ[(100.0, 0.5)]["ratio_np_over_first"], NBQ[(1e4, 0.5)]["ratio_np_over_first"], NBQ[(1e4, 0.25)]["db_np"], NBQ[(1e4, 0.5)]["db_np"],
             NBQ[(1e4, 1.0)]["db_np"], 100 * SRJ["interpolation_max_dev_small_rs"], texsci(SR1["screening"], 0))))
md = lambda k: MSR[k]["db"]
add(("p", "Several sources were treated with a three-dimensional bounce (multisource_bounce.py): a periodic row of discrete sources along a line, "
          "spacing $d$, coupling $q=%.0f$ and radius %.0f each, so that the field depends on Euclidean time, the position along the row and the distance from "
          "it, against the same coupling per unit length smeared uniformly along the line. The bubble is free to deform and to sit between the sources. "
          "Without sources, the three-dimensional code gives $b=%.1f$ against %.1f from shooting. For $d=8$, 12 and 20 the discrete row increases the action by "
          "%.1f, %.1f and %.1f, the smeared line by %.1f, %.1f and %.1f: the same mean coupling suppresses the decay far less when it is concentrated in "
          "discrete sources, and once the spacing exceeds the size of the bubble the discrete row hardly matters while the smeared line still does. "
          "This finding supports the first-order conclusion beyond first order: discrete sources suppress the decay much less than a homogeneous mean-field "
          "treatment predicts. What is established is a constraint within this coupling and vacuum model, for a single source and for a row of sources; a "
          "fully three-dimensional bounce among many sources at random positions is not computed, and general many-body configurations remain open."
          % (MSJ["q"], MSJ["rs"], MSJ["b0_3d"], MSJ["b0_shooting"], md(("discrete", 8.0)), md(("discrete", 12.0)), md(("discrete", 20.0)),
             md(("smeared", 8.0)), md(("smeared", 12.0)), md(("smeared", 20.0)))))
add(("fig", "fig3_resolved_granularity.png", "Fig. 10. (a) Static field in a Wigner--Seitz cell around one dark-matter source at the trigger coupling, "
            "for cells from 0.2 to 65 Compton lengths: small cells reproduce the mean-field vacuum, large cells the zero-density vacuum away from the "
            "source. (b) Suppression of the decay rate by randomly placed point sources, relative to the first-order mean-field suppression, against the "
            "particle mass at the trigger density (first order); shaded, the wave-dark-matter window of Section 3.3 and the particle regime. (c) Increase in the "
            "bounce action from a source at the bubble centre: nonperturbative bounce against first order, source radius 0.5 (triangles: radii 0.25 and 1).", 6.9))

# ============================================================================================ 4 nucleation
add(("h1", "4. Density-dependent bubble nucleation"))
add(("h2", "4.1 The bounce and its action"))
add(("p", "Where the mean field applies, the O(4)-symmetric bounce solves, in units of $f$ for the field and $(\\sqrt\\lambda f)^{-1}$ for lengths,"))
add(("eq", r"x''+\frac{3}{r}x'=\partial_xp(x;\eta,\theta),\qquad x'(0)=0,\qquad x(\infty)=x_h,\qquad b(\eta,\theta)=2\pi^2\int_0^\infty r^3\left[\frac{x'^2}{2}+p(x)-p(x_h)\right]dr"))
add(("p", "and the action is $B=b/\\lambda$ \\cite{Coleman,CC}. Gravity enters only through $g=f^2/M_{\\rm Pl}^2$, of order $10^{-60}$ at the dark-energy "
          "scale, so the flat-space bounce is exact for all practical purposes \\cite{CDL,Parke}; Section 6 confirms that the wall's own gravity is negligible. "
          "The bounce is found by overshoot and undershoot on the central value with an eighth-order Runge--Kutta integrator, where the central value lies "
          "within $10^{-5}$ of the true minimum (the thin-wall regime) shooting becomes exponentially sensitive, and collocation on the boundary-value problem "
          "is used instead (vacuum_wall.py). At $\\theta=0$ the actions reproduce an independent table to $%s$, and where both methods apply, they agree to the "
          "tolerance of the collocation."
          % texsci(max(abs(v["b"] / v["table"] - 1) for v in DBJ["check_theta0"].values()) + 1e-16, 0)))
add(("h2", "4.2 The wall"))
add(("p", "The interface follows from the field equations rather than being imposed. At the degeneracy density the static planar wall, $x''=\\partial_xp$, has "
          "the first integral $z(x)=\\int dx/\\sqrt{2(p-p_{\\rm min})}$, tension $\\sigma=\\int\\sqrt{2(p-p_{\\rm min})}\\,dx=%.3f\\,\\sqrt\\lambda f^3$ and 10--90%% "
          "thickness $%.1f/(\\sqrt\\lambda f)$ for $\\eta=0.30$; the symmetric case $\\eta=\\theta=0$ reproduces $\\sigma=2/3$ and $x=\\tanh(z/2)$. Below the "
          "degeneracy density a static planar wall does not exist, since the energy difference accelerates it, and the wall is read from the bounce: its tension "
          "$\\sigma_b=\\int x'^2dr$ and its thickness change little with density, whereas the radius increases from $%.1f$ to $%.1f$ and the action from $b=%.0f$ to "
          "$%.0f$ between $\\theta=0$ and $0.6\\,\\theta_{\\rm deg}$ (Table 2, Fig. 3a,b). Thin-wall action"
          % (deg["sigma"], deg["thickness"], r0["R_w"], r6["R_w"], r0["b"], r6["b"])))
add(("eq", r"B_{\rm tw}=\frac{27\pi^2\sigma^4}{2\,\Delta V^3},\qquad \sigma=\int_{x_t}^{x_h}\sqrt{2\,[p(x)-p_{\rm lin}(x)]}\,dx"))
add(("p", "with $p_{\\rm lin}$ the straight line through the two minima, which removes the energy difference from the tension, is recovered within 1%% from "
          "$0.6\\,\\theta_{\\rm deg}$ upward for $\\eta=0.20$ and 0.30, and overestimates the thick-wall bounce near zero density ($b/B_{\\rm tw}=%.2f$ at "
          "$\\theta=0$ for $\\eta=0.30$). The density dependence of the action is therefore dominated by $\\Delta V(\\theta)$, which decreases linearly to zero at "
          "$\\theta_{\\rm deg}$, while the tension is nearly constant." % r0["b_over_B_tw"]))
rows2 = [["%.2f" % r["theta_over_deg"], "%.3f" % r["dV"], "%.0f" % r["b"], "%.2f" % r["b_over_B_tw"], "%.1f" % r["R_w"], "%.2f" % r["thickness"],
          "%.3f" % r["sigma_bounce"], r["method"]] for r in TB if abs(r["theta_over_deg"] * 10 - round(r["theta_over_deg"] * 10)) < 1e-6]
add(("table", "Table 2. Bounce and wall against density ($\\eta=0.30$, units $\\lambda f^4$, $(\\sqrt\\lambda f)^{-1}$ and $\\sqrt\\lambda f^3$).",
     ["$\\theta/\\theta_{\\rm deg}$", "$\\Delta V$", "$b$", "$b/B_{\\rm tw}$", "$R_w$", "thickness", "$\\sigma_b$", "method"], rows2, None))
add(("h2", "4.3 Physical scale, survival and the trigger"))
add(("p", "We normalize to the observed universe: the false vacuum carries the observed dark-energy density $V_f$, the true vacuum $V_t=-qV_f$ with $q=1$, "
          "so that $\\lambda f^4=(1+q)V_f/\\Delta\\tilde V_0$, and the action at zero density is $B_0=b(0)/\\lambda=%.0f$ for $\\eta=0.30$. Then $f=%.2f$ meV and the "
          "wall scale $1/\\sqrt\\lambda f$ is %.3f mm. The rate per unit four-volume is $\\Gamma=\\mathcal A\\,e^{-B(\\theta)}$ with the one-loop prefactor "
          "estimated as $(B/2\\pi)^2$ times the fourth power of the false-vacuum mass scale \\cite{CC}, and the relevant comparison is with the Hubble "
          "four-volume, $\\Gamma/H^4$. The probability that no bubble has nucleated in the past light cone of a present observer is $\\exp(-\\int\\Gamma\\,dV_4)$; "
          "in the observed $\\Lambda$CDM history the four-volume of the past light cone is $V_4=%.4f\\,H_\\Lambda^{-4}$, %.0f times smaller than the Euclidean "
          "four-sphere often used as a shortcut, so that survival requires only $\\Gamma/H_\\Lambda^4\\lesssim%.0f$ for a constant rate \\cite{GW,TWW}. With the "
          "coupling the rate in the past was suppressed by the larger density, and the survival integral is evaluated with the density history "
          "$\\theta=\\theta_0(1+z)^3$. In the future, the rate increases as the matter dilutes (Fig. 3c), and the transition completes, in the sense that bubbles "
          "percolate, once $\\Gamma/H^4>9/4\\pi$ \\cite{GW}. The prefactor is the least certain input: a full one-loop determinant \\cite{DM} could "
          "differ from the estimate by a factor comparable to $(B/2\\pi)^2$ itself, $e^{%.1f}$ at the trigger. Because the rate increases at $\\beta=%.0f\\,H$ "
          "(Section 6), such a factor moves the completion by %.2f Hubble times, and it moves the boundary of the region with a protected past by the same "
          "%.1f units in $B_0$; the qualitative results do not depend on it."
          % (WJ["base"]["B0"], PS["f_meV"], PS["length_unit_mm"], SV["V4_HL4"], SV["ratio_euclid_to_lorentzian"], SV["gamma_over_hl4_max_95"],
             2 * np.log(TW["B"] / (2 * np.pi)), PC["beta_over_H"], 2 * np.log(TW["B"] / (2 * np.pi)) / PC["beta_over_H"], 2 * np.log(TW["B"] / (2 * np.pi)))))
add(("p", "Because the field responds to the local density (Section 3.3), these statements hold region by region, with $\\theta_0$ the local value. For "
          "a region with $\\theta_0=%.1f\\,\\theta_{\\rm deg}$, the baseline of this paper, the past light cone is clear, and the transition is complete at "
          "$a=%.2f$, when $\\theta=%.3f\\,\\theta_{\\rm deg}$ and $B=%.0f$; the nucleated bubble has radius %.2f mm, wall thickness %.2f mm and tension "
          "%.0f meV$^3$."
          % (F0_LOC, WJ["a_trigger"], TW["theta_over_deg"], TW["B"], TW["R_w_mm"], TW["thickness_mm"], TW["sigma_meV3"])))
add(("h2", "4.4 The environment: voids first, halos last"))
add(("p", "The rate depends exponentially on the local density, so the emptiest regions decide both whether the past is protected and when the decay "
          "begins. The deepest cosmic voids have central densities of 10--20%% of the mean \\cite{SvdW04}. If a region with $1+\\delta$ times the mean "
          "density is followed through the whole past light cone with that contrast (conservative, since voids were less empty in the past), the past light "
          "cone of a region at the mean density is clear for any mean coupling above $\\theta_0=%.2f\\,\\theta_{\\rm deg}$, but a region at %.0f%% or %.0f%% "
          "of the mean requires $\\theta_0\\ge%.2f$ or $%.2f\\,\\theta_{\\rm deg}$ on average (Fig. 2b; granularity.py). With the mean coupling of Section 4.3 "
          "taken as the mean rather than the void value, regions at %.0f%% of the mean or less would already have decayed. A protected past therefore requires a "
          "mean coupling today close to the degeneracy density, and the baseline of Section 4.3 is to be read as the local value in voids at %.0f%% of the "
          "mean, $\\theta_0=%.1f\\,\\theta_{\\rm deg}$ on average."
          % (VOID[1.0]["min_f0_mean"], 100 * 0.2, 100 * 0.1, VOID[0.2]["min_f0_mean"], VOID[0.1]["min_f0_mean"],
             100 * max(k for k, v in VOID.items() if not v["baseline_past_clear"]), 100 * DV_, F0_MEAN)))
add(("p", "With that mean, the held vacuum is still the true vacuum wherever $1+\\delta>\\theta_{\\rm deg}/\\bar\\theta(a)$: at the trigger this is any region "
          "denser than %.2f times the mean, and the threshold increases as $a^3$. A bubble nucleated in a void expands through all connected regions in which the "
          "held vacuum is metastable, and its wall stalls where the local density reaches the degeneracy value, because the driving pressure, proportional to "
          "$\\Delta V_{\\rm eff}$, vanishes there and changes sign beyond. The transition therefore advances from the voids into sheets and filaments as the "
          "universe expands, following the contour $\\theta_{\\rm loc}=\\theta_{\\rm deg}$, while virialized halos, whose density does not dilute, keep the "
          "held vacuum indefinitely: a galactic halo at about two hundred times the mean density has $\\theta\\approx%.0f\\,\\theta_{\\rm deg}$, far above "
          "$\\theta_{\\rm high}=%.2f\\,\\theta_{\\rm deg}$, where only the held minimum exists. The dark energy measured inside galaxies is therefore unchanged "
          "by the transition, and the halos become islands of the held vacuum bounded by walls whose tension energy is negligible (%s kg per square meter "
          "for the baseline tension). Separate voids are separated by such walls at the trigger and must each nucleate a bubble; since the rate grows at "
          "$\\beta=%.0f\\,H$ (Section 6), a void of comoving radius 5--50 Mpc does so within %.3f--%.3f Hubble times of the global completion criterion."
          % (THR_T, 200 * F0_MEAN, GRJ["environment"]["theta_high_over_deg"],
             texsci(TW["sigma_meV3"] * 1e-9 / (1.97327e-7) ** 2 * 1.602e-19 / (2.998e8) ** 2, 0).join("$$"),
             PC["beta_over_H"], VCOMP["50 Mpc"]["delay_hubble_times"], VCOMP["5 Mpc"]["delay_hubble_times"])))
fr = lambda a: a ** 3 / F0_MEAN
add(("p", "The advance of the transition can be read off the threshold $1+\\delta_*(a)=a^3\\theta_{\\rm deg}/\\bar\\theta_0$: %.2f at the trigger, %.1f at "
          "$a=2$, %.1f at $a=3$ and %.0f at $a=5$. Sheets and the outskirts of filaments, a few times denser than the mean, are swept within a Hubble time of the "
          "trigger; the cores of filaments, ten or more times the mean, only later; virialized structures never, because their density does not fall. Two "
          "features distinguish this from vacuum decay seeded by inhomogeneities elsewhere in the literature. Compact objects such as black holes catalyze "
          "the decay of a metastable vacuum by lowering the action in their vicinity \\cite{GMW14,BGM15}; here the seed is the absence of matter, and the "
          "decay is suppressed, not catalyzed, where matter has collapsed. The wall does not need to overcome the denser regions: it stops where the "
          "held vacuum becomes the true vacuum, so that the outcome of the transition is a map of the cosmic density field at the time of the trigger."
          % (fr(WJ["a_trigger"]), fr(2.0), fr(3.0), fr(5.0))))
add(("p", "The environment yields a quantitative relation. Because the field responds to the local density and the trigger depends on the coupling only "
          "through $\\theta_0a^{-3}$, a region whose density is a fraction $1+\\delta$ of the mean completes its transition at "
          "$a_{\\rm tr}(\\delta)=a_{\\rm tr}(0)\\,(1+\\delta)^{1/3}$ while its contrast is frozen. For a mean coupling $\\bar\\theta_0=%.0f\\,\\theta_{\\rm deg}$ the "
          "computed values are $a_{\\rm tr}=%.2f$, %.2f, %.2f and %.2f for $1+\\delta=0.1$, 0.3, 1 and 2 (viability_extended.py), and a change in the one-loop "
          "prefactor by $(B/2\\pi)^2$ moves them by at most %.0f%%. The order of the transition across the cosmic web, voids first and in the order of their "
          "emptiness, is a prediction of the mechanism. It is not observable today, because a protected past means that no transition has occurred inside "
          "our past light cone; what is observable now is listed in Section 10."
          % (F0_MEAN, PRED[(3.0, 0.1)]["central"], PRED[(3.0, 0.3)]["central"], PRED[(3.0, 1.0)]["central"], PRED[(3.0, 2.0)]["central"], 100 * PREF)))
EV = {round(r["one_plus_delta_today"], 2): r for r in EEV["rows"]}
add(("p", "The relation assumes a frozen contrast, but voids keep emptying. With the contrast evolved by the $\\Lambda$CDM growth factor and the "
          "nonlinear evolution of a spherical void, $1+\\delta=(1-\\delta_L/\\nu)^{-\\nu}$ with $\\nu=%.1f$ (env_evolution.py), regions at 10--50%% of the mean "
          "density today complete their transition %.0f--%.0f%% earlier in scale factor than with a frozen contrast, and a region at twice the mean %.0f%% later; "
          "the order is preserved and the relation steepens slightly. The growth of structure after the trigger is modest because the expansion is then "
          "dominated by the vacuum energy (the growth factor increases by %.0f%% between $a=1$ and 2)."
          % (EEV["nu"], 100 * min(abs(EV[g]["shift"]) for g in (0.1, 0.2, 0.3, 0.5)), 100 * max(abs(EV[g]["shift"]) for g in (0.1, 0.2, 0.3, 0.5)),
             100 * EV[2.0]["shift"], 100 * ([d["D"] for d in EEV["D_at"] if d["a"] == 2.0][0] - 1))))
add(("fig", "fig3_interface_properties.png", "Fig. 3. (a) Wall tension (solid) and thickness/10 (dashed) of the O(4) bounce against density for three tilts. "
            "(b) Bounce action $b=\\lambda B$, numerical (points) and thin-wall (dotted). (c) $\\ln(\\Gamma/H^4)$ along the dilution history of a region with "
            "$\\theta_0=%.1f\\,\\theta_{\\rm deg}$, with the completion threshold and the trigger." % F0_LOC, 6.6))

# ============================================================================================ 5
add(("h1", "5. Interface tension and field transmission"))
sp_ = TJ["species"]
g1 = sp_["1.00"]
add(("p", "A wall makes the mass of every field coupled to $\\varphi$ depend on position. For a test scalar with $\\mathcal L_{\\rm int}=-\\frac12g^2\\varphi^2\\psi^2$, "
          "$m_\\psi^2(z)=m_0^2+g^2f^2x(z)^2$, and the modes of the planar wall in its rest frame obey"))
add(("eq", r"\left[-\frac{d^2}{dz^2}+m_\psi^2(z)\right]u=(\omega^2-k_\parallel^2)\,u,\qquad T=\frac{k_+}{k_-}\frac{1}{|A|^2},\qquad R=\left|\frac{B}{A}\right|^2"))
add(("p", "with a transmitted wave $e^{ik_+z}$ on the false-vacuum side, integrated back across the wall and decomposed into $Ae^{ik_-z}+Be^{-ik_-z}$ on the "
          "true-vacuum side (interface_transport.py). The same equation describes the wave dark matter of Section 3.3, whose mass changes by the factor "
          "$A(x_t)/A(x_h)$ across the wall. The solver is validated on the symmetric kink, for which $m^2=g^2f^2-g^2f^2\\,{\\rm sech}^2(z/2)$ is a "
          "Poeschl--Teller well with the exact transmission $T=\\sinh^2(2\\pi q)/[\\sinh^2(2\\pi q)+\\cos^2(\\frac\\pi2\\sqrt{1+16g^2/\\lambda})]$; the numerical "
          "values agree with it to %s, including the reflectionless coupling $g^2/\\lambda=1/2$, and $T+R=1$ holds to %s (Fig. 4b)."
          % (sci(TJ["max_PT_error"], 0), sci(TJ["max_unitarity_error"], 0))))
add(("p", "Through the wall at the trigger the interior is heavier for every coupling ($m_{\\rm in}/m_{\\rm out}=%.2f$, since $|x_t|>|x_h|$). Modes below the "
          "interior threshold are reflected, and those above it $T$ exceeds 0.5 within $0.1\\sqrt\\lambda f$ of the threshold for every coupling tested (Fig. 4a). "
          "The interface is selective, but kinematically: it filters by mass, as any step does, with little additional reflection from its shape. Photons, "
          "neutrinos and baryons do not couple to $\\varphi$ and cross unchanged. A thermal bath of test scalars exerts on a static wall a reflected momentum "
          "flux that is at most %s of the driving pressure at temperatures up to the wall scale, and in the ultrarelativistic limit at most "
          "$\\Delta m^2T^2/24$ per degree of freedom \\cite{BM09}; at the CMB temperature of the trigger epoch, this value is below %s of the driving pressure."
          % (g1["m_in"] / g1["m_out"], sci(max(r["P_static_over_dV"] for r in TJ["pressure"] if r["T"] <= 0.2), 0), sci(0.31 * TJ["T_cmb_at_trigger_units"] ** 2, 0))))
add(("p", "For the coupled dark matter, the mass ratio across the wall is $A(x_t)/A(x_h)=%.2f$, far greater than its kinetic energy; thus, a slow wall would "
          "reflect all of it, as at the filtering transitions of \\cite{BKL20}. The wall is not slow. Its driving pressure exceeds the maximum pressure the "
          "dark matter can exert, $\\rho_*\\Delta A$ (the Boedeker--Moore limit for quanta that are transmitted at high Lorentz factor \\cite{BM09}), by a "
          "factor %.0f, so it runs away; dark matter at rest crosses once the wall's Lorentz factor exceeds $A(x_t)/A(x_h)$, which occurs within %.2f mm of the "
          "nucleation point. The fraction of the dark matter swept by the wall before it becomes transparent is $(\\gamma_{\\rm th}R_0/R_{\\rm coll})^3=%s$. "
          "The interface is therefore transparent to the dark matter for all practical purposes, and the transitioned region inherits the dark matter of its "
          "surroundings with the heavier mass of the true vacuum."
          % (DM["A_true"] / DM["A_false"], 1.0 / DM["P_dm_max_over_dV_vac"], DM["R_threshold_mm"], texsci(DM["filtered_fraction"], 0))))
rA = DM["A_true"] / DM["A_false"]
dep_ad, dep_su = DM["A_true"] - DM["A_false"], DM["A_false"] * (1 + rA * rA) / 2 - DM["A_false"]
add(("p", "For wave dark matter, the crossing was simulated directly (wave_wall_1d.py). A planar wall, starting at rest as a kink, accelerates through a "
          "homogeneous oscillating dark-matter field whose mass is $mA(x)$, with the full backreaction of the dark matter on the scalar and nothing "
          "prescribed about the crossing (Fig. 11). The wall follows hyperbolic motion, and a fit of it gives the Lorentz factor; at the end of the run "
          "this is %.1f%% and %.1f%% below that of a run with the static mean-field term for $m_{\\rm dm}=10m_\\varphi$ and $3m_\\varphi$, a difference that "
          "changes by %.1f%% when the grid spacing is doubled, while the wall position changes by %.1f%%. The energy given to the dark matter per unit length swept is %.2f times its initial density (%.2f for "
          "$3m_\\varphi$), close to the sudden-crossing value $(1+r^2)/2=%.2f$, with $r=A(x_t)/A(x_h)=%.2f$, and well above the adiabatic value $r$: once the "
          "wall is relativistic, the crossing is sudden, and the field maintains its amplitude and velocity rather than its number of quanta. The dark matter is "
          "not left at rest. Quanta that cross a relativistic wall into the heavier vacuum leave it with momentum $m(r^2-1)/2$ in the frame of the "
          "background and velocity $(r^2-1)/(r^2+1)=%.2f$; the simulated mean velocity of the dark matter behind the wall is %.2f, and the dark matter forms "
          "a shell that trails the wall with its energy density raised by $1/(1-v)$ (%.0f against %.1f predicted). The energy taken from the wall is %.2f of "
          "its driving pressure, so within the friction mechanisms considered here (Table 3) the wall still runs away. In the transitioned region, the dark matter carries %.0f%% more quanta than the adiabatic "
          "estimate and moves outward until the shells of neighboring bubbles collide, after which its momentum redshifts away. Table S8 lists the "
          "diagnostics at four resolutions and two box sizes. The arrival time of the wall at a fixed position changes by %s between the two finest "
          "resolutions, the Lorentz factor relative to the mean-field run by %.1f%%, the deposited energy by %.1f%%, the velocity of the dark matter by "
          "%.2f%% and the shell compression by %.0f%%; the absolute fitted Lorentz factor changes by %.0f%%, because a wall contracted to a few grid cells is "
          "underresolved, which is why only ratios to the mean-field run are used. The larger box gives identical results, as it must, since nothing from "
          "its far end can reach the wall during the run. The simulation is planar and in flat space; it does not establish the same behavior for a "
          "spherical wall in an expanding universe, where the wall's Lorentz factor is far larger and the crossing correspondingly more sudden."
          % (100 * abs(W10["gamma_rel_diff_vs_mean_field_end"]), 100 * abs(W3["gamma_rel_diff_vs_mean_field_end"]),
             100 * abs(WWR[("wave", 10.0, 0.01)]["gamma_end_fit"] / WWR[("mean field", 10.0, 0.01)]["gamma_end_fit"] - 1.0
                       - W10["gamma_rel_diff_vs_mean_field_end"]),
             100 * abs(WWJ["convergence"]["wave_z_end_rel"]), W10["deposit_per_length_over_rho"], W3["deposit_per_length_over_rho"],
             WWJ["params"]["sudden"], rA, WWJ["params"]["v_quanta"], W10["v_dm_behind_late"], SHELL10, WWJ["params"]["shell_compression"],
             DM["P_dm_max_over_dV_vac"] * dep_su / dep_ad, 100 * ((1 + rA * rA) / (2 * rA) - 1),
             texsci(abs(CONV_F[-1]["t_at_z150"] / CONV_F[-2]["t_at_z150"] - 1), 0).join("$$"),
             100 * abs(CONV_F[-1]["gamma_rel_mean_field"] - CONV_F[-2]["gamma_rel_mean_field"]),
             100 * abs(CONV_F[-1]["deposit"] / CONV_F[-2]["deposit"] - 1), 100 * abs(CONV_F[-1]["v_dm"] / CONV_F[-2]["v_dm"] - 1),
             100 * abs(CONV_F[-1]["shell"] / CONV_F[-2]["shell"] - 1), 100 * abs(CONV_F[-1]["gamma_fit_end"] / CONV_F[-2]["gamma_fit_end"] - 1))))
add(("fig", "fig4_selective_transmission.png", "Fig. 4. (a) Transmission of a test scalar through the wall at the trigger against momentum above the "
            "interior threshold, for four couplings. (b) Validation of the symmetric kink: numerical transmission (dots) against the exact Poeschl--Teller "
            "result (lines), including the reflectionless coupling $g^2/\\lambda=1/2$.", 6.3))
add(("fig", "fig6_wave_wall.png", "Fig. 11. A planar wall running through wave dark matter (1+1 dimensions, backreaction included). (a) Lorentz factor of "
            "the wall without dark matter, with the static mean-field term, and with wave dark matter of mass $10m_\\varphi$ and $3m_\\varphi$. (b) "
            "Dark-matter energy behind the wall against the wall position, with the adiabatic and sudden slopes.", 6.6))

# ============================================================================================ 6
add(("h1", "6. Gravitational wall and collision dynamics"))
add(("p", "For a thin wall of tension $\\sigma$ and areal radius $R$ between an interior of density $\\rho_-$ and an exterior of density $\\rho_+$, the Israel "
          "junction condition can be written with the Misner--Sharp mass on each side \\cite{Israel66,BGG}:"))
add(("eq", r"\Gamma_--\Gamma_+=\frac{\sigma R}{2M_{\rm Pl}^2},\qquad \Gamma_\pm=\sqrt{1+\dot R^2-\frac{\rho_\pm R^2}{3M_{\rm Pl}^2}}\ \Rightarrow\ \Gamma_-=\frac{(\rho_+-\rho_-)R}{3\sigma}+\frac{\sigma R}{4M_{\rm Pl}^2}"))
add(("p", "(bubble_gravity.py). From $R_0=3\\sigma/\\Delta V$ (%.2f mm) to the Hubble radius the wall follows the flat-space runaway $\\Gamma=R/R_0$ to "
          "round-off and reaches $\\Gamma=$%s; its self-gravity, $\\sigma R/2\\Gamma$, is %s; and the Misner--Sharp mass it carries equals the vacuum-energy "
          "deficit of the region it encloses (Fig. 5a). The released vacuum energy is therefore carried by the wall, not deposited in the interior, until walls "
          "collide or stall at the boundary of a denser region."
          % (GJ["wall"]["R0_mm"], sci(GJ["wall"]["Gamma_at_RH"], 1), sci(GJ["wall"]["max_grav"], 0))))
add(("p", "The completion of the transition in a connected region is computed from its rate history (bubble_handoff.py): the false-vacuum fraction "
          "$f=e^{-I}$, $I=\\frac{4\\pi}{3}\\int\\Gamma a^3r^3dt$ with $r$ the comoving radius of a bubble nucleated at $t$ \\cite{GW,TWW}, decreases to one half at "
          "$a=%.3f$ and to 1%% at $a=%.3f$ (Fig. 6a). The rate grows at $\\beta=d\\ln\\Gamma/dt=%.0f\\,H$, so the transition is completed in %.3f Hubble times and "
          "the bubbles collide at a mean radius $R_{\\rm coll}=(8\\pi)^{1/3}/\\beta=%.4f/H$, comparable to the size of large voids. During a collision, the walls "
          "have Lorentz factors on the order of %s."
          % (PC["a_half"], PC["a_one_percent"], PC["beta_over_H"], PC["completion_duration_over_hubble_time"], PC["R_coll_times_H"],
             sci(WB["gamma_at_collision"], 0))))
gmin = min(r["gamma_terminal"] for r in WB["gauge_nlo_terminal_gamma"] if r["T_meV"] < 1.0 and r["g"] <= 1.0)
add(("p", "The junction treatment omits several effects, which are bounded instead (Table 3). The Lorentz-contracted thickness of the wall is %s of its "
          "radius, so the thin-wall description is excellent once the wall has grown. The pressure of the dark matter is at most %.2f of the drive. A dark gauge "
          "boson acquiring mass from $\\varphi$ would add a pressure that increases with the Lorentz factor \\cite{BM17} and stop the runaway; for couplings "
          "$g\\le1$ and any bath colder than 1 meV the terminal Lorentz factor is at least %.0f, still above the threshold %.1f at which the dark matter "
          "crosses, so the transparency of Section 5 is unchanged. Gravitational waves from the collisions carry at most a fraction $\\sim%s$ of the wall "
          "energy, suppressed by $(H/\\beta)^2$, and since $R_{\\rm coll}H\\approx%.3f$ no region of Hubble size is overdense enough to form black holes."
          % (sci(WB["contracted_thickness_over_radius"], 0), WB["dm_pressure_over_vacuum_drive"], gmin, WB["gamma_threshold_dm"],
             texsci(WB["gw_fraction_bound"], 0), WB["pbh_scale_ratio"])))
add(("p", "None of these signals are observable today. The trigger lies in the future, at $a\\approx%.1f$ for the baseline, and a protected past means "
          "precisely that no bubble has nucleated within our past light cone; the imprints of bubble collisions on the microwave background that have been "
          "studied for transitions in the early universe \\cite{Kleban,AJ} therefore do not arise. What is observable now is indirect: the dark energy is "
          "exactly $w=-1$ while the vacuum is held, and the dark matter must be of the kind specified in Section 3."
          % WJ["a_trigger"]))
add(("table", "Table 3. Approximations of the wall dynamics and of the hand-off to the enclosed region, and their bounds (bubble_handoff.py, granularity.py).",
     ["Effect", "Bound or result", "Consequence"],
     [["mean-field source", "wave dark matter, %.0f meV $\\ll m_{\\rm dm}\\ll$ %.0f meV" % tuple(GRJ["wave_window_meV"]), "required (Section 3)"],
      ["completion of the transition", "%.3f Hubble times; $R_{\\rm coll}H=%.4f$" % (PC["completion_duration_over_hubble_time"], PC["R_coll_times_H"]),
       "void by void"],
      ["equation of state of collision products", "$w=1/3$ and $w=0$ bracketed", "turnaround moves by %.0f%%" %
       (100 * abs(W0["a_max_over_a_nuc"][1] / W13["a_max_over_a_nuc"][1] - 1))],
      ["thin wall after contraction", "thickness/radius %s" % texsci(WB["contracted_thickness_over_radius"], 0).join("$$"), "valid"],
      ["self-gravity of the wall", "$\\sigma R/2\\Gamma=%s$" % texsci(GJ["wall"]["max_grav"], 0), "negligible"],
      ["dark-matter pressure", "%.2f of the drive (adiabatic), %.2f (sudden crossing)" % (WB["dm_pressure_over_vacuum_drive"],
                                                                                      WB["dm_pressure_over_vacuum_drive"] * dep_su / dep_ad), "runaway"],
      ["dark-matter energy after crossing", "sudden: +%.0f%% quanta; simulated %.0f%% of the sudden energy" % (100 * ((1 + rA * rA) / (2 * rA) - 1),
                                                                                                             100 * W10["deposit_per_length_over_rho"] / WWJ["params"]["sudden"]),
       "turnaround +%.1f%%" % (100 * (RSR["sudden"]["a_max_over_a_nuc"] / RSR["adiabatic"]["a_max_over_a_nuc"] - 1))],
      ["wave dark matter against the mean field", "wall $\\gamma$ %.1f%% ($10m_\\varphi$) and %.1f%% ($3m_\\varphi$) lower" %
       (100 * abs(W10["gamma_rel_diff_vs_mean_field_end"]), 100 * abs(W3["gamma_rel_diff_vs_mean_field_end"])), "mean field holds"],
      ["gauge-boson (NLO) pressure", "terminal $\\Gamma\\ge%.0f$ ($g\\le1$, $T<1$ meV)" % gmin, "dark matter still crosses"],
      ["gravitational waves", "fraction $\\lesssim%s$" % texsci(WB["gw_fraction_bound"], 0), "negligible"],
      ["black holes from collisions", "$R_{\\rm coll}H=%.4f$" % WB["pbh_scale_ratio"], "none"],
      ["isolated Hubble-scale bubble", "trapped surface at $t=%.2f$ (%.2f with sudden crossing), before the interior bounce at %.2f" %
       (IS["t_trapped"], RSI["sudden, observed"]["t_trapped"], IS["t_bounce"]), "trapped region"]], None))
add(("fig", "fig5_daughter_evolution.png", "Fig. 5. (a) Lorentz factor of the wall from the Israel junction against the flat-space runaway, and its "
            "self-gravity (right axis). (b) Scale factor of a transitioned region continued formally through an effective bounce (the conditional "
            "continuation of Section 7.3, not supported by the dark matter of Section 3). (c) The same equations for the open region enclosed by an isolated "
            "Hubble-scale bubble, which becomes trapped before its restoration (Section 7.2).", 6.6))

# ============================================================================================ 7
add(("h1", "7. Energy accounting and the fate of the transitioned region"))
add(("h2", "7.1 The region enclosed by the wall"))
nB, tB, rB, bB = ev(B_P, "nucleation"), ev(B_P, "turnaround"), ev(B_P, "restoration"), ev(B_P, "bounce")
tI = ev(B_I, "turnaround")
add(("p", "A wall says nothing by itself about the cosmology inside it; the interior has to be solved. It is evolved as a Friedmann region with radiation, "
          "baryons, dark matter of field-dependent mass $A(x)$ and the scalar in one of its minima: the scalar mass, of order meV, exceeds $H$ by thirty orders "
          "of magnitude, so the field follows its minimum adiabatically, and the stress-energy is $\\rho=\\rho_r+\\rho_b+\\rho_*+U$, "
          "$\\rho+P=\\frac43\\rho_r+\\rho_b+\\rho_*A(x)$, with $U(x,\\theta)$ the vacuum energy including the coupling energy $\\rho_*(A-1)$. The junction fixes "
          "the initial data: the fluid inside maintains the Hubble rate of the parent at nucleation while its density decreases by $\\Delta V_{\\rm eff}$, so an isolated "
          "Hubble-scale bubble encloses an open region with $K/a^2=\\Delta V_{\\rm eff}/3$, and a completed transition leaves a flat region whose wall energy is "
          "released at the collisions. Events (nucleation, turnaround, restoration) are found as roots in $\\ln a$ with $H$ from the constraint, so the "
          "constraint holds exactly between events; an independent integration as an initial-value problem reproduces the turnaround and the restoration to "
          "%s in time with a constraint residual of %s."
          % (sci(abs(RP["restoration_time_ivp"] - RP["restoration_time_root"]) / RP["restoration_time_root"], 0), sci(RP["constraint_residual_ivp"], 0))))
add(("p", "In the environmental picture of Section 4.4, the transitioned volume is a network of voids and, later, sheets and filaments, which together fill "
          "most of the space, while the halos keep the held vacuum. Its volume-averaged vacuum energy is negative once it fills most of the space, and a region in "
          "which the transition is complete holds the negative true vacuum, turns around at %.2f times the scale factor at nucleation and contracts. When the "
          "dark-matter density in the contraction reaches $\\theta_{\\rm high}$ the true minimum disappears and the field returns to the matter-held minimum, "
          "releasing $%.2f\\,V_f$; the matter term restores the held vacuum, a direct consequence of the density dependence. This restoration needs "
          "dissipation: undamped, the oscillation about the held minimum is an adiabatic invariant that returns to the same amplitude when the density falls "
          "back through $\\theta_{\\rm high}$. The inhomogeneous evolution of a universe with islands of the held vacuum is not computed here; the "
          "volume-averaged statement is that it recollapses. With the dark-matter energy of a sudden crossing (Section 5; region_sensitivity.py) the "
          "turnaround moves from %.3f to %.3f times the nucleation scale factor and the restoration occurs earlier in the contraction; part of that energy "
          "is kinetic and redshifts away, so the two values bracket the evolution."
          % (tB["a_max_over_a_nuc"], rB["dE_released"] / VC, RSR["adiabatic"]["a_max_over_a_nuc"], RSR["sudden"]["a_max_over_a_nuc"])))
add(("h2", "7.2 Isolated Hubble-scale bubbles develop trapped surfaces"))
add(("p", "A bubble that remains isolated on the Hubble scale, which happens when the transition does not complete, encloses an open region that also "
          "turns around (at %.2f times the nucleation scale factor), but its continuation is not part of the parent universe. Following the wall as a null ray "
          "in the open interior and comparing the Misner--Sharp mass inside it with that of the exterior false-vacuum universe at the same areal radius, the "
          "interior becomes a trapped region ($\\rho_{\\rm in}R^2/3>1$) at $t=%.2f/H_0$ and denser than the exterior at $t=%.2f/H_0$, both before its "
          "restoration ($%.2f$) and before any bounce ($%.2f$) (Fig. 6b). Within the thin-wall junction model, the interior then contains a trapped surface, "
          "an apparent horizon, and the junction requires $\\Gamma_+$ to change sign, as for the false-vacuum bubbles of \\cite{BGG,FG,FGG}. Whether a global "
          "event horizon forms, and hence a black hole, depends on the subsequent evolution, which is not computed. The trapping time does not depend on "
          "the exterior cosmology (it is unchanged for $\\Omega_\\Lambda\\pm0.02$ and with the coupling energy of the held vacuum in the exterior) and moves "
          "to %.2f when the dark matter receives the sudden-crossing energy, still before the restoration (%.2f) and the bounce (%.2f) of that run "
          "(region_sensitivity.py). Void-sized regions, much smaller than the Hubble radius, are not trapped in this way: their trapping parameter is smaller "
          "by $(R_{\\rm void}H)^2\\sim10^{-4}$, and they follow the averaged evolution."
          % (tI["a_max_over_a_nuc"], IS["t_trapped"], IS["t_interior_denser"], IS["t_restoration"], IS["t_bounce"],
             RSI["sudden, observed"]["t_trapped"], RSI["sudden, observed"]["t_restoration"], RSI["sudden, observed"]["t_bounce"])))
add(("fig", "fig6_handoff.png", "Fig. 6. (a) Fraction of space in the false vacuum along the dilution history of a connected region, with the completion "
            "criterion $\\Gamma/H^4=9/4\\pi$. (b) Isolated Hubble-scale bubble: density inside the wall and in the exterior, and the trapping parameter "
            "$\\rho_{\\rm in}R^2/3$, with the turnaround, restoration and bounce of the interior.", 6.4))
add(("h2", "7.3 Energy and entropy closure, and why a continuation is not supported"))
add(("p", "Every species crosses the wall, so the comoving entropy of the region is inherited, and it is increased when the wall energy is thermalized at "
          "the collisions (by a factor of about %.0f) and when the vacuum is restored (bubble_entropy.py); it never decreases. The interior bookkeeping "
          "$dE_{\\rm in}=F_E\\,dV+W_{\\rm wall}$, with the swept-in content $F_E=\\rho_+$ and the work $W_{\\rm wall}=-\\Delta V_{\\rm eff}dV$ handed to the "
          "wall, closes to round-off. The decay removes the cosmological horizon of the false vacuum, whose entropy ($%s$) exceeds that of the radiation in a "
          "Hubble volume ($%s$); this is the usual situation for decays out of de Sitter space, accounted for by the action itself."
          % (EJ["entropy"][B_P]["nucleation_jumps"][1], texsci(EJ["horizon"]["S_dS"], 1), texsci(EJ["horizon"]["S_radiation_Hubble_volume"], 1))))
add(("p", "The collision products are not computed; their equation of state is bracketed: if they behave as radiation, the region evolves exactly as with "
          "thermalized wall energy, and if they are massive scalar quanta ($w=0$) the turnaround moves from %.2f to %.2f times the nucleation scale factor. "
          "Whether the recollapsing region could be continued through a nonsingular bounce into a new cycle is a separate question. With an effective "
          "loop-quantum-cosmology bounce \\cite{APS06,AS11} the equations do repeat, self-similarly with growth $\\ln G=%.3f$ per cycle ($%.3f$ for $w=0$ "
          "products), but only if the bounce returns the dark matter to the state that holds the vacuum in the next cycle. For a thermal relic, a hot bounce "
          "does that through freeze-out; for the wave dark matter that Section 3 requires it does not, because a meV-scale boson thermalized at a hot bounce is "
          "a hot relic, not a cold coherent condensate. The dark matter that allows the trigger is therefore not the dark matter that would allow a cyclic "
          "continuation, and the latter is not supported in this branch."
          % (W13["a_max_over_a_nuc"][1], W0["a_max_over_a_nuc"][1], LNG["B"], float(np.mean(W0["lnG"])))))

# ============================================================================================ 8
add(("h1", "8. Parameter-space analysis and numerical controls"))
add(("p", "Controls separate the contribution of each ingredient (bubble_to_cycle.py, granularity.py; Table 4). Without the coupling (D), there is no "
          "density dependence, and a single transition that has neither already occurred in our past light cone nor fails to complete requires $B_0$ in a "
          "window of width %.0f around %.0f. With the coupling but a granular source, particle dark matter or wave dark matter lighter than the field, the "
          "rate is the zero-density rate, and the vacuum decays in the past. With the computed transmission (C), the region is identical to that with every "
          "species transmitted (B) to %s: the selectivity of the wall changes nothing, and a hypothetical wall that keeps all the dark matter out (C*) differs "
          "only in the conditional continuation."
          % (DW[1] - DW[0], 0.5 * (DW[0] + DW[1]), sci(DM["filtered_fraction"], 0))))
rows4 = [["D", "no matter coupling", "constant rate", "window $%.0f\\le B_0\\le%.0f$" % tuple(DW), "fine tuning"],
         ["G1", "coupling, particle dark matter", "zero-density rate between particles", "$e^{%.0f}$ bubbles in the past" %
          GRJ["particle_dm"]["ln_expected_bubbles_past_light_cone"], "excluded"],
         ["G2", "coupling, wave dark matter $m_{\\rm dm}\\ll m_\\varphi$", "zero-density rate at the nodes", "as G1", "excluded"],
         ["B", "coupling, wave dark matter in the window, all species transmitted", "local mean field", "past clear in voids for $\\bar\\theta_0\\gtrsim%.1f\\,"
          "\\theta_{\\rm deg}$" % VOID[0.1]["min_f0_mean"], "viable"],
         ["C", "as B with the computed transmission", "local mean field", "identical to B", "viable"],
         ["C*", "fully selective wall (hypothetical)", "local mean field", "differs only in a continuation", "--"]]
add(("table", "Table 4. Controls ($\\eta=0.30$, $q=1$, $B_0=%.0f$; voids at %.0f%% of the mean density)." % (WJ["base"]["B0"], 100 * DV_),
     ["Model", "Description", "Release", "Past light cone", "Outcome"], rows4, None))
add(("p", "In the plane of $B_0$ and the local coupling today in the emptiest regions, the region in which the past light cone is clear and the transition "
          "completes extends from $B_0=%.0f$ to $%.0f$ for suitable couplings and covers %.0f%% of the plane scanned (Fig. 7); above the upper end, the vacuum "
          "rate itself is too low for completion. A survey over stated ranges (viability_mc.py; %d draws, log-uniform $B_0\\in[50,1000]$, local "
          "$\\theta_0/\\theta_{\\rm deg}\\in[10^{-3},10]$, uniform $q\\in[-1,3]$, $\\eta\\in\\{0.20,0.30,0.36\\}$) finds a clear past light cone in %.0f%% of draws, "
          "completion in %.0f%% and a negative true vacuum, needed for recollapse, in %.0f%%; the three hold together in %.0f%%. These fractions describe the "
          "sample space of the model, not a probability for our universe, and they presuppose the wave dark matter of Section 3; for particle dark matter, "
          "the fraction is zero."
          % (CJ["map"]["B0_range_viable"][0], CJ["map"]["B0_range_viable"][1], 100 * CJ["map"]["fraction_viable"], TM["n"],
             100 * TM["single"]["T1 past light cone clear"], 100 * TM["single"]["T2 completes"], 100 * TM["single"]["T3 recollapse"],
             100 * TM["joint_background"])))
add(("p", "These fractions depend on the sampling distributions and do not yet include the conditions of Sections 3 and 5. An extended survey "
          "(viability_extended.py; two sets of %d draws) adds the dark-matter mass, the depth of the emptiest regions and the runaway of the wall with the "
          "sudden-crossing energy, and evaluates every condition in the emptiest regions. With broad distributions (log-uniform $B_0\\in[50,1000]$, local "
          "$\\theta_0/\\theta_{\\rm deg}\\in[10^{-3},10]$, $q\\in[-1,3]$, $1+\\delta_v\\in[0.05,0.5]$, $m_{\\rm dm}$ from 1 $\\mu$eV to 1 GeV) the background "
          "conditions hold in %.0f%% of draws, adding the mean-field condition ($5m_\\varphi\\le m_{\\rm dm}\\le m_{\\rm cl}$) leaves %.1f%%, and adding the "
          "runaway leaves %.2f%%; with distributions concentrated in the region of interest ($B_0\\in[150,400]$, $q\\in[0.1,3]$, $m_{\\rm dm}$ from 1 meV "
          "to 1 eV) the fractions are %.0f%%, %.0f%% and %.0f%%. The joint fraction varies little with the depth of the emptiest regions (%.1f--%.1f%% across "
          "the terciles of the broad set). The mean-field condition is the binding condition; the aim is not a large fraction but a region that survives every "
          "condition, and it does (Table S9)."
          % (VEB["n"], 100 * VEB["joint_background"], 100 * VEB["joint_background_meanfield"], 100 * VEB["joint_all"], 100 * VEF["joint_background"],
             100 * VEF["joint_background_meanfield"], 100 * VEF["joint_all"], 100 * min(x["joint_all"] for x in VEB["joint_all_by_void_depth"]),
             100 * max(x["joint_all"] for x in VEB["joint_all_by_void_depth"]))))
PR_ = VEF["viable_physical_ranges"]
rg = lambda k, f="%.2f": (f + "--" + f) % (PR_[k]["q10"], PR_[k]["q90"])
add(("p", "Because these fractions depend on how the ranges and distributions are chosen, the viable region is better described by its physical "
          "parameters. Among the draws of the focused set that satisfy every condition, the central 80%% have $B_0=%s$, a local coupling today in the "
          "emptiest regions $\\theta_0/\\theta_{\\rm deg}=%s$, $q=%s$, voids at %s of the mean density, a scalar mass $m_\\varphi=%s$ meV, a dark-matter mass "
          "$m_{\\rm dm}=%s$ meV and completion in the emptiest regions at $a=%s$ (Table S9); all three tilts contribute (%s)."
          % (rg("B0", "%.0f"), rg("local theta_0/theta_deg"), rg("q"), rg("1+delta_v"), rg("m_phi (meV)", "%.1f"), rg("m_dm (meV)", "%.0f"),
             rg("a at completion (emptiest regions)"), ", ".join("%d for $\\eta=%s$" % (v, k) for k, v in VEF["viable_by_eta"].items()))))
add(("p", "All %d checks of interface_validation.py pass (Table S1): the analytic degeneracy and window edges, the kink tension and thickness, the bounce "
          "actions against an independent table, the thin-wall limit, shooting against collocation, the transmission against Poeschl--Teller and the sharp "
          "step, unitarity and tolerance convergence, the flat-space limit of the junction and the energy carried by the wall, the constraint and continuity "
          "of the region evolution, entropy monotonicity, the controls, the completion time and the independent integration. The granularity estimates "
          "(granularity.py) and the resolved tests (granular_field.py, wave_wall_1d.py) are tabulated with their convergence checks in Tables S4--S8." % VJ["n"]))
add(("fig", "fig7_viability_map.png", "Fig. 7. Outcome over the vacuum action $B_0$ and the local coupling today in the emptiest regions ($\\eta=0.30$, "
            "$q=1$): already decayed in the past light cone, no completed transition, or a completed transition; the red band is the window without "
            "coupling (D).", 4.6))

# ============================================================================================ 9 neutrinos
TT = NUJ["terrestrial"]
add(("h1", "9. Optional neutrino--scalar coupling and the CMB"))
add(("p", "The dark-matter-only model does not predict neutrinos. If the scalar also coupled to neutrinos, a new hypothesis, the neutrino masses "
          "would track the local field. As an illustrative benchmark, consider a flavor-universal Dirac coupling $\\mathcal L=-y\\varphi\\bar\\nu\\nu$, which "
          "requires an electroweak completion (for instance, a higher-dimensional operator with the Higgs field and right-handed neutrinos) that is not "
          "specified here; the masses shift by $\\delta m_\\nu=yfx_{\\rm loc}$, with $x_{\\rm loc}$ the held minimum at the local dark-matter coupling "
          "(neutrino_coupling.py)."))
add(("p", "For the homogeneous local benchmark and the flavor-universal coupling assumed here, terrestrial oscillation experiments would see none "
          "of it. The scalar couples to dark matter only, and at the dark-matter density of the solar "
          "neighborhood, about $%s$ times the cosmic mean, the field sits at $1-x=%s$. Along a 1300 km baseline such as that of DUNE \\cite{DUNE20}, and "
          "along the baselines of Hyper-Kamiokande \\cite{HK18}, the halo density changes by a relative $%s$ and the field by less than $10^{-15}$. The shift "
          "is then a constant addition to the mass matrix, absorbed into the measured masses and mixing, so that $\\Delta P_{\\alpha\\beta}(E)=0$ "
          "identically; an interface, across which a neutrino would accumulate a different oscillation phase, does not exist inside a protected past light "
          "cone. This null result applies to the coupling considered; flavor-dependent couplings, or couplings to a field that varies on terrestrial "
          "scales, are not covered."
          % (texsci(NUE["solar neighbourhood"]["one_plus_delta"], 1), texsci(TT["one_minus_x_local"], 1), texsci(TT["relative_density_change_over_baseline"], 0))))
add(("p", "The coupling would instead make the neutrino masses depend on environment and on time, as in mass-varying-neutrino models \\cite{FNW04}. "
          "Relative to the solar neighborhood the masses would be smaller by $%.2fy$ meV at the mean density and $%.2fy$ meV in the deepest voids today, and "
          "the mean-density value would be larger by $%.2fy$ meV at recombination than it is at present (Table 5). For a perturbative coupling, these shifts are below a "
          "meV, beneath the sensitivity of the present cosmological and laboratory mass measurements; a larger coupling would require a separate assessment of "
          "neutrino self-interaction and free-streaming constraints, which is not attempted. Relic neutrinos also act back on the trigger, adding "
          "$n_\\nu m_\\nu(\\varphi)$ to the effective potential, a tilt $\\delta\\eta=%s\\,y$, so the trigger is unaffected for $y\\lesssim%.0f$ (a shift of 1%% "
          "of $\\eta$). Inside a transitioned region, the masses would shift by $%.1fy$ meV. Neutrino production in wall collisions and a neutrino-driven "
          "trigger are separate questions and are not addressed."
          % (abs(NUE["mean density"]["dm_nu_over_y_meV"]), abs(NUE["deepest voids (10% of mean)"]["dm_nu_over_y_meV"]),
             [r for r in NUJ["redshift_mean_density"] if r["z"] == 1100.0][0]["dm_nu_over_y_meV"], texsci(NUJ["backreaction"]["delta_eta_per_y"], 1),
             NUJ["backreaction"]["y_for_delta_eta_1pct_of_eta"], NUJ["after_transition"]["dm_nu_over_y_meV"])))
wz = CMB["worst_in_window"]
add(("p", "The cosmic microwave background constrains the coupling much more strongly than any of these effects could be measured. In the early "
          "universe, the dark-matter density was greater by $(1+z)^3$, and with it, the mass of the scalar: $m_\\varphi\\approx%.0f$ eV at $z=10^3$ against a "
          "neutrino temperature of %.2f eV. The scalar is then a heavy mediator for neutrino--neutrino scattering, with a contact cross-section "
          "$\\sigma\\sim y^4s/16\\pi m_\\varphi^4$, and since $m_\\varphi\\propto(1+z)^{3/2}$ grows faster than the temperature, the ratio of the scattering rate to "
          "the expansion rate rises toward late times. The CMB requires neutrinos to stream freely while the modes it measures enter the horizon "
          "\\cite{CRS14,AH14}; demanding $\\Gamma_\\nu<H$ for $10^3<z<10^5$ (neutrino_coupling.py) gives the approximate theoretical bound $y\\lesssim%s$, "
          "set at $z\\approx10^3$. It is a free-streaming criterion, not an observational confidence limit, which would need a consistent collision term "
          "and a comparison with CMB data through a Boltzmann code. Its assumptions hold over the window: the scalar is heavier than the momentum transfer "
          "at every redshift (Table S16), so the contact form applies, and the neutrinos are relativistic; below $z\\sim10^2$, where they become "
          "nonrelativistic, the estimate does not apply. The scalar is always at least %.0f times heavier than the neutrino temperature, so its thermal "
          "production is expected to be strongly suppressed in the modeled regime, subject to initial conditions and other production channels. With "
          "$y\\lesssim%s$ the environmental and redshift-dependent mass shifts of Table 5 are below $%s$ meV, so the optional "
          "coupling, within this bound, has no observable consequence. The CMB tests the scalar here only through this optional coupling; a computation of "
          "the full spectra with a Boltzmann code is not attempted."
          % (wz["m_phi_meV"] / 1e3, wz["T_nu_meV"] / 1e3, texsci(CMB["y_max_free_streaming"], 0), CMB["min_m_over_T"],
             texsci(CMB["y_max_free_streaming"], 0), texsci(CMB["max_mass_shift_meV_at_y_max"], 0))))
rows5 = [[k.replace("neighbourhood", "neighborhood"), "%s" % (texsci(v["one_plus_delta"], 1).join("$$") if v["one_plus_delta"] > 1e3 else "%g" % v["one_plus_delta"]), "%.4f" % v["x_held"],
          "%.3f" % v["dm_nu_over_y_meV"]] for k, v in NUE.items()]
rows5 += [["mean density at $z=%g$" % r["z"], "--", "%.4f" % r["x_held"], "%+.3f" % (r["dm_nu_over_y_meV"] + NUE["mean density"]["dm_nu_over_y_meV"])]
          for r in NUJ["redshift_mean_density"] if r["z"] in (1.0, 3.0, 1100.0)]
add(("table", "Table 5. Optional neutrino coupling: the held minimum and the neutrino-mass shift per unit coupling relative to the solar neighborhood "
              "($\\bar\\theta_0=%.0f\\,\\theta_{\\rm deg}$)." % F0_MEAN, ["Environment", "$1+\\delta$", "$x_{\\rm loc}$", "$\\delta m_\\nu/y$ (meV)"], rows5, None))

# ============================================================================================ 10
add(("h1", "10. Implications, limitations and conclusions"))
add(("p", "A coupling of a scalar to dark matter centered on a false vacuum makes the stability of that vacuum a function of the dark-matter density, and "
          "dilution changes the held vacuum from the true vacuum to a metastable vacuum at $\\theta_{\\rm deg}=\\eta/2$. For a vacuum at the dark-energy scale this "
          "mean-field statement holds only if the dark matter looks smooth to the field on a scale below a tenth of a millimeter and a time below a picosecond. "
          "In a spatially resolved treatment, confirmed beyond first order by a bounce that deforms around a source, particle dark matter does not, and with "
          "it the vacuum decays at its zero-density rate; the coincidence of the "
          "dark-matter and dark-energy densities, which motivates a dilution trigger, is what makes it fail. Wave dark matter heavier than the field and still classical, approximately 3 meV to 0.25 eV, "
          "does, and then the decay is environmental. It begins in the emptiest voids, the past is protected only if the mean coupling today is close to the "
          "degeneracy density, the transition advances from voids into sheets and filaments as the universe expands, and galactic halos keep the held "
          "vacuum. Within the voids, the wall is a fraction of a millimeter thick, runs away, becomes transparent to the dark matter and carries the released "
          "energy, while the dark matter it crosses is kicked forward and trails it; the transitioned regions recollapse in their negative vacuum, and isolated "
          "Hubble-scale bubbles develop trapped surfaces."))
add(("p", "The observational implications follow from the environment. The dark energy measured inside halos is unchanged until the global recollapse; "
          "thus, local tests reveal no transition, whereas the volume-averaged expansion responds to the voids. The dark matter must be a nonthermal boson in "
          "the meV to sub-eV range with a coupling to the scalar of strength $M\\approx f$, which changes its mass by a factor %.1f between the vacua; such a "
          "coupling acts only over the Compton length of the scalar and is invisible to astrophysical fifth-force tests. The dark energy is exactly $w=-1$ "
          "until the trigger, and the order of the future transition follows the density of the cosmic web, $a_{\\rm tr}\\propto(1+\\delta)^{1/3}$. Within "
          "the resolved treatment of Section 3.5, a detection of thermal-relic dark matter, or of dark matter heavier than approximately one eV, would exclude the "
          "density-triggered decay of a dark-energy-scale vacuum in this class of models."
          % (DM["A_true"] / DM["A_false"])))
add(("p", "The limitations are specific. The granular decay rate is computed at first order in each source's coupling; the nonperturbative test, "
          "which shows that first order overstates the suppression, is performed for a single source on the bubble axis, not for many sources at arbitrary "
          "positions, and the conclusion is a constraint within this coupling and vacuum model, not a universal statement about particle dark matter. The wave-dark-matter crossing is simulated in 1+1 dimensions for a planar wall, not for an expanding "
          "bubble. The local density is treated as constant along the past light cone of a region, which is conservative for voids; "
          "the void statistics enter only through the depth of the emptiest regions; the advance of the walls into denser regions and the islands of held "
          "vacuum are described through the degeneracy contour, not simulated; and the trapped surfaces of isolated bubbles are found within the thin-wall "
          "junction model, without establishing a global event horizon. The one-loop prefactor is estimated, not computed, and enters the completion time only "
          "logarithmically; the collision products are bracketed by their equation of state; and the next-to-leading friction is estimated for a hypothetical "
          "dark gauge boson only."))
add(("p", "The results do not require a cosmological bounce, and they do not support one in this branch. They establish when a density-dependent vacuum "
          "can decay at all, where it decays first, what its walls transmit and how they move, and what the regions they enclose become. The homogeneous, "
          "classical version of the trigger, which does not face the granularity bound because its field has a cosmological Compton length, is examined in the "
          "companion paper \\cite{Companion}."))

# ============================================================================================ declarations
add(("h1", "Declarations"))
add(("pb", "Funding.", "This research received no specific grant from funding agencies in the public, commercial, or not-for-profit sectors."))
add(("pb", "Declaration of competing interest.", "The author declares no known competing financial interests or personal relationships that could have "
        "influenced the work reported in this paper."))
add(("pb", "CRediT authorship contribution statement.", "Leon Sandler: Conceptualization, Methodology, Software, Formal analysis, Investigation, Writing -- "
        "original draft, Writing -- review and editing."))
add(("pb", "Data and code availability.", "All code (Python, NumPy, SciPy, Matplotlib \\cite{NumPy,SciPy}) and the data behind every figure and table are at %s "
        "(MIT license) and archived at https://doi.org/%s; this manuscript is archived at https://doi.org/%s. Input tables shared with the companion paper "
        "are documented as such in the README." % (GITHUB, DOIS.get("code", "(to be assigned)"), DOIS.get("paper", "(to be assigned)"))))
add(("pb", "Related manuscript.", "A companion manuscript on the homogeneous version of the trigger and its continuation into cycles has been prepared for "
        "General Relativity and Gravitation%s; it shares the class of models but none of the results reported here."
        % ((", archived at https://doi.org/%s" % COMP["paper"]) if COMP.get("paper") else "")))
add(("pb", "Declaration of generative AI and AI-assisted technologies in the writing process.", "During the preparation of this work, the author used Claude "
        "(Anthropic) to assist with the numerical analysis, coding, literature checking and drafting. After using this tool, the author reviewed and edited the "
        "content as needed and takes full responsibility for the content of the publication."))
add(("h1", "References"))

_fig_order = [int(re.match(r"Fig\. (\d+)\.", b[2]).group(1)) for b in B if b[0] == "fig"]
FIGMAP = {old: i + 1 for i, old in enumerate(_fig_order)}
_fr = lambda t: re.sub(r"\bFig\. (\d+)", lambda m: "Fig. %d" % FIGMAP.get(int(m.group(1)), int(m.group(1))), t) if isinstance(t, str) else t
B = [(b[0], _fr(b[1])) + tuple(b[2:]) if b[0] in ("p", "abstract") else ((b[0], b[1], _fr(b[2])) + tuple(b[3:]) if b[0] in ("fig", "pb") else b)
     for b in B]
add = B.append
order = []
for blk in B:
    chunks = [blk[-1]] if blk[0] in ("p", "abstract", "pb") else (list(blk[1]) if blk[0] == "bullets" else
                                                                  ([str(c) for row in blk[3] for c in row] if blk[0] == "table" else []))
    for ch in chunks:
        for m in re.finditer(r"\\cite\{([^}]+)\}", ch):
            for key in m.group(1).split(","):
                if key not in order:
                    order.append(key)
undefined = [k for k in order if k not in REFS]
if undefined:
    sys.exit("undefined citation keys %s" % undefined)
numbering = {k: i + 1 for i, k in enumerate(order)}
rep = lambda s: re.sub(r"\\cite\{([^}]+)\}", lambda m: "[" + ", ".join(str(numbering[x]) for x in m.group(1).split(",")) + "]", s)
for k in order:
    add(("ref", "[%d] %s" % (numbering[k], REFS[k])))
B2 = []
for b in B:
    if b[0] in ("p", "abstract", "note"):
        B2.append((b[0], rep(b[1])))
    elif b[0] == "pb":
        B2.append(("pb", b[1], rep(b[2])))
    elif b[0] == "table":
        B2.append(("table", b[1], b[2], [[rep(str(c)) for c in row] for row in b[3]], b[4]))
    else:
        B2.append(b)
B = [CE.block(b) for b in B2]
BASENAME = "PDU_Vacuum_Bubbles"


def supplement():
    S = []
    a = S.append
    a(("title", "Supplementary material: " + TITLE))
    a(("author", "Leon Sandler\nIndependent Researcher, Northbrook, Illinois, USA"))
    a(("h1", "S1. Methods"))
    a(("p", "granularity.py evaluates the scales of Section 3 (the field mass at the held vacuum, the Compton length, the dark-matter density, the particle, "
            "classicality and coherence bounds), the zero-density rate and past-light-cone integral seen between particles, the structural ratio over the "
            "bounce tables and toward $\\eta_c$ (with new bounces computed by vacuum_wall.bounce_profile), the protection of the past in regions of reduced "
            "density, and the void-by-void completion delay. vacuum_phase_map.py maps the stationary points over density and tilt and evaluates the coherence "
            "of the symmetron release. vacuum_wall.py solves the O(4) bounce with its profile, the thin-wall tension with the energy difference removed, the "
            "static wall at the degeneracy density, and the nucleation history with the actual Hubble rate. interface_transport.py integrates the mode equation "
            "across the planar wall with complex amplitudes (DOP853, rtol $10^{-10}$) and computes the reflected momentum flux of a Bose bath with the "
            "transverse integral in closed form. bubble_gravity.py evaluates the Israel junction and evolves the enclosed regions by root finding in $\\ln a$; "
            "bubble_to_cycle.py runs the controls and the viability map; bubble_entropy.py the energy and entropy closure; bubble_handoff.py the percolation "
            "integral, the bounds of Table 3, the isolated-bubble test, the collision-product variants and the independent integration; viability_mc.py the "
            "survey; interface_validation.py the checks of Table S1. granular_field.py solves the static field in a cell by collocation (solve_bvp, tolerance "
            "$10^{-7}$) and evaluates the first-order granular rate from the zero-density bounce at full resolution; wave_wall_1d.py integrates the two "
            "coupled fields with a leapfrog scheme ($\\Delta t=\\Delta z/2$, reflecting ends) at two resolutions; region_sensitivity.py repeats the region "
            "evolution and the trapping test with the sudden-crossing energy and other exterior cosmologies; viability_extended.py runs the extended survey "
            "and the prediction; env_evolution.py repeats the prediction with evolving contrasts; nonpert_bounce.py finds the bounce with a source on its axis "
            "by constrained minimization (L-BFGS, quadratic penalty on the four-volume beyond the barrier top, continuation in the constraint) and solves "
            "the static field of one source on a fine radial grid; screened_rate.py recomputes the granular rate with self-screened couplings; "
            "neutrino_coupling.py evaluates the optional neutrino coupling and its early-universe history."))
    a(("table", "Table S1. Validation (interface_validation.py).", ["Check", "Result"], IB.validation_rows(), None))
    rows = [["%.2f" % r["gh"], "%.2f" % r["T"], texsci(max(r["P_static_over_dV"], 1e-300), 1).join("$$"), texsci(max(r["P_BM_over_dV"], 1e-300), 1).join("$$")]
            for r in TJ["pressure"]]
    a(("table", "Table S2. Pressure of a thermal test-scalar bath on the wall at the trigger relative to the driving pressure: static wall and "
                "ultrarelativistic limit $\\Delta m^2T^2/24$. Temperatures in units of $\\sqrt\\lambda f$; the CMB at the trigger is %.3f."
                % TJ["T_cmb_at_trigger_units"], ["$g/\\sqrt\\lambda$", "T", "static", "relativistic limit"], rows, None))
    rows = [[v["variant"], ", ".join("%.4f" % x for x in v["lnG"]), "%.2f" % v["T_b_MeV"], ", ".join("%.3f" % x for x in v["a_max_over_a_nuc"])]
            for v in HJ["variants"]]
    a(("table", "Table S3. Formal continuation variants (not supported by the dark matter of Section 3; Section 7.3): growth per cycle, bounce temperature "
                "($\\rho_r=T^4$) and turnaround relative to nucleation.", ["Variant", "$\\ln G$", "$T_b$ (MeV)", "$a_{\\max}/a_{\\rm nuc}$"], rows, None))
    rows = [["%.3f" % r["eta"], "%.1f" % r["theta_over_deg"], "%.1f" % r["q"], "%.1f" % r["b"], texsci(r["barrier"], 1).join("$$"), "%.3f" % r["lam"],
             "%.2f" % r["ratio"]] for r in GRJ["near_eta_c"]]
    a(("table", "Table S4. Approach to $\\eta_c$: bounce action, barrier, $\\lambda=b/B$ and $\\rho_{\\rm dm}/m_\\varphi^4$ at the trigger epoch.",
       ["$\\eta$", "$\\theta/\\theta_{\\rm deg}$", "$q$", "$b$", "barrier", "$\\lambda$", "$\\rho_{\\rm dm}/m_\\varphi^4$"], rows, None))
    rows = [["%.1f" % v["one_plus_delta"], "%.1f" % v["baseline_lnI"], "yes" if v["baseline_past_clear"] else "no", "%.3f" % v["min_f0_mean"]]
            for v in GRJ["voids"]]
    a(("table", "Table S5. Regions of reduced density: ln of the expected number of bubbles in the past light cone for a mean coupling "
                "$\\theta_0=%.1f\\,\\theta_{\\rm deg}$, and the smallest mean coupling for which the past is clear." % F0_LOC,
       ["$1+\\delta$", "ln expected bubbles", "past clear", "minimum mean $\\theta_0/\\theta_{\\rm deg}$"], rows, None))
    rows = [["%.2f" % c["Rc_times_mphi"], "%.2f" % c["rs"], "%.4f" % c["x_center"], "%.4f" % c["x_cell_average"], "%.3f" % c["fraction_at_zero_density_vacuum"],
             "%.3f" % c["fraction_at_mean_field_vacuum"]] for c in GF["cells"]]
    rows += [["%.2f (source %.1f)" % (c["Rc_times_mphi"], c["rs"]), "%.2f" % c["rs"], "%.4f" % c["x_center"], "%.4f" % c["x_cell_average"],
              "%.3f" % c["fraction_at_zero_density_vacuum"], "%.3f" % c["fraction_at_mean_field_vacuum"]] for c in GF["cells_rs_check"]]
    a(("table", "Table S6. Static field in a Wigner--Seitz cell around one source (granular_field.py): field at the source and cell average, and the "
                "fraction of the cell within 1%% of the zero-density and of the mean-field vacuum (zero-density %.4f, mean-field %.4f). The last rows vary the "
                "radius of the source." % (GF["cells"][0]["x_h_zero_density"], GF["cells"][0]["x_h_mean_field"]),
       ["$R_cm_\\varphi$", "source radius", "$x$ at source", "$\\langle x\\rangle$", "zero-density", "mean-field"], rows, None))
    rows = [[texsci(r["m_dm_meV"], 1).join("$$"), texsci(r["spacing_over_compton"], 1).join("$$"), "%.4g" % r["ln_rate_ratio"],
             texsci(max(r["fraction_of_first_order_mean_field_suppression"], 1e-300), 1).join("$$")] for r in GF["rate_vs_mass"][::2]]
    a(("table", "Table S7. Decay rate with randomly placed point sources at the trigger density: $\\ln(\\Gamma/\\Gamma_0)$ and its fraction of the "
                "first-order mean-field suppression (%.2f; the full mean-field value is %.2f). First-order $\\partial_\\theta b$ = %.1f against %.1f from the "
                "numerical bounce." % (GF["mean_field"]["ln_rate_ratio_first_order"], GF["mean_field"]["ln_rate_ratio_full"],
                                       GF["bounce"]["db_dtheta_first_order"], GF["bounce"]["db_dtheta_numerical"]),
       ["$m_{\\rm dm}$ (meV)", "spacing / Compton length", "$\\ln(\\Gamma/\\Gamma_0)$", "fraction of mean field"], rows, None))
    rows = [["%.4f" % c["dz"], "%.0f" % c["L_box"], "%.4f" % c["t_at_z150"], "%.1f" % c["gamma_fit_end"], "%.4f" % c["gamma_rel_mean_field"],
             "%.4f" % c["deposit"], "%.4f" % c["v_dm"], "%.2f" % c["shell"]] for c in CONV]
    a(("table", "Table S8. Convergence of the 1+1-dimensional wall through wave dark matter ($m_{\\rm dm}=10m_\\varphi$): arrival time of the wall at "
                "$z=150$, fitted Lorentz factor at the end of the run and its ratio to the mean-field run, energy deposited per unit length (sudden %.3f), "
                "mean dark-matter velocity behind the wall (predicted %.4f) and shell compression (predicted %.1f), for four grid spacings and two box sizes."
                % (WWJ["params"]["sudden"], WWJ["params"]["v_quanta"], WWJ["params"]["shell_compression"]),
       ["$\\Delta z$", "box", "$t(z=150)$", "$\\gamma$", "$\\gamma/\\gamma_{\\rm MF}-1$", "deposit", "$v_{\\rm dm}$", "shell"], rows, None))
    rows = []
    for r in WWJ["runs"]:
        rows.append([r["mode"], "%.0f" % r["m_ratio"] if r["mode"] == "wave" else "--", "%.3f" % r["dz"], "%.2f" % r["z_wall"][-1], "%.1f" % r["gamma_end_fit"],
                     "%.2f" % r["deposit_per_length_over_rho"] if r.get("deposit_per_length_over_rho") else "--",
                     "%.3f" % r["v_dm_behind_late"] if r.get("v_dm_behind_late") else "--"])
    a(("table", "Table S8b. All runs of the 1+1-dimensional wall through wave dark matter (wave_wall_1d.py): final wall position and Lorentz factor, energy deposited per "
                "unit length (in units of the initial dark-matter density; adiabatic %.2f, sudden %.2f) and the mean dark-matter velocity behind the wall "
                "(predicted %.3f). Lorentz factors from a fit of hyperbolic motion to the wall position; at the coarser resolution the fitted values "
                "change by about 6%% but their ratio to the mean-field run by less than 0.5%%, and the wall position by %.2f%%."
                % (WWJ["params"]["adiabatic"], WWJ["params"]["sudden"], WWJ["params"]["v_quanta"], 100 * abs(WWJ["convergence"]["wave_z_end_rel"])),
       ["run", "$m/m_\\varphi$", "$\\Delta z$", "$z_{\\rm wall}$", "$\\gamma$", "deposit", "$v_{\\rm dm}$"], rows, None))
    rows = []
    for sv in VEJ["surveys"]:
        rows.append([sv["kind"], "%d" % sv["n"]] + ["%.3f" % sv["single"][k] for k in ("T1", "T2", "T3", "T4", "T5")] +
                    ["%.3f" % sv["joint_background"], "%.4f" % sv["joint_background_meanfield"], "%.4f" % sv["joint_all"]])
    a(("table", "Table S9. Extended survey (viability_extended.py): single fractions T1 (past clear), T2 (completion), T3 (recollapse), T4 (mean field), "
                "T5 (runaway with the sudden-crossing energy), and the joint fractions.",
       ["priors", "draws", "T1", "T2", "T3", "T4", "T5", "T1-3", "+T4", "+T5"], rows, None))
    rows = [["%.0f" % r["mean_f0"], "%.1f" % r["one_plus_delta"], "%.3f" % r["central"], "%.3f--%.3f" % (r["prefactor_up"], r["prefactor_down"]),
             "yes" if r["past_clear"] else "no"] for r in VEJ["prediction"]]
    a(("table", "Table S10. Scale factor at which a region of density $(1+\\delta)$ times the mean completes its transition, with the range for a change "
                "of the one-loop prefactor by $(B/2\\pi)^{\\pm2}$ ($\\eta=0.30$, $q=1$, $B_0=%.0f$; the grid of scale factors has a spacing of 2.3%%)."
                % WJ["base"]["B0"], ["mean $\\theta_0/\\theta_{\\rm deg}$", "$1+\\delta$", "$a_{\\rm tr}$", "prefactor range", "past clear"], rows, None))
    rows = [[k] + ["%.3f" % v[x] if v.get(x) is not None else "--" for x in ("t_turnaround", "t_trapped", "t_interior_denser", "t_restoration", "t_bounce")]
            for k, v in RSJ["isolated"].items()]
    a(("table", "Table S11. Isolated Hubble-scale bubble: turnaround, trapping, interior denser than the exterior, restoration and bounce (in $H_0^{-1}$) "
                "for the adiabatic and sudden crossings and three exterior cosmologies (region_sensitivity.py).",
       ["case", "turnaround", "trapped", "denser", "restoration", "bounce"], rows, None))
    rows = [["%g" % b["q"], "%.2f" % b["rs"], "%.2f" % b["b"], "%.2f" % b["db_np"], "%.2f" % b["db_first_order"],
             "%.4f" % b["ratio_np_over_first"] if b["ratio_np_over_first"] else "--", "%.4f" % b["ubg_center"]] for b in NBJ["bounces"]]
    a(("table", "Table S12. Nonperturbative bounce with a source on the bubble axis (nonpert_bounce.py; grid spacing 0.4, zero-density bounce %.2f "
                "against %.2f from shooting): action, its increase, the first-order estimate, their ratio and the background field at the source."
                % (NBJ["b0_grid"], NBJ["b0_shooting"]),
       ["$q$", "$r_s$", "$b$", "$\\Delta b$", "first order", "ratio", "$x_{\\rm bg}(0)$"], rows, None))
    rows = [["%s" % rs, "%g" % c["q"], "%.4f" % c["q_eff"], "%.4f" % (c["q"] / (1 + c["q"] / (4 * np.pi * float(rs)))), "%.1e" % c["resolution_change"],
             "yes" if c["success"] else "no"] for rs, v in NBJ["screening"].items() for c in v]
    a(("table", "Table S13. Static coupling of a single source of radius $r_s$ (Yukawa tail of the nonlinear solution, Newton minimization at two radial "
                "resolutions) against the thin-shell interpolation $q/(1+q/4\\pi r_s)$; a case is accepted if both resolutions converge and agree to 1%.",
       ["$r_s$", "$q$", "$q_{\\rm eff}$", "interpolation", "resolution change", "accepted"], rows, None))
    rows = [["%.2f" % c["h"], "%.0f" % c["L"], "%.2f" % c["b0"], "%.2f" % c["db_q10"], "--", "%.2f" % c["db_q1e4"]] for c in NBJ["convergence"]]
    sf = NBJ["saturation_fine"]
    rows += [["0.25", "20", "%.2f" % sf["h025"]["b0"], "--", "%.2f" % sf["h025"]["db_q1e3"], "%.2f" % sf["h025"]["db_q1e4"]],
             ["0.20", "20", "%.2f" % sf["h02"]["b0"], "--", "%.2f" % sf["h02"]["db_q1e3"], "%.2f" % sf["h02"]["db_q1e4"]]]
    a(("table", "Table S12b. Convergence of the nonperturbative bounce (source radius 0.5): zero-density action (shooting %.2f) and the increase of the "
                "action for $q=10$, $10^3$ and $10^4$ against grid spacing and box size. The saturated increase falls as the source is resolved and is "
                "an upper bound." % NBJ["b0_shooting"],
       ["grid spacing", "box", "$b_0$", "$\\Delta b(q=10)$", "$\\Delta b(q=10^3)$", "$\\Delta b(q=10^4)$"], rows, None))
    rows = [[r["cfg"]["sources"], "%g" % r["cfg"].get("d", 0) if r["cfg"].get("d") else "--", "%.2f" % r["b"], "%.2f" % r["db"],
             "%.3f--%.3f" % (r["ubg_min_on_axis"], r["ubg_max_on_axis"])] for r in MSJ["runs"]]
    a(("table", "Table S18. Three-dimensional bounce with a periodic row of discrete sources ($q=%.0f$, radius %.0f, spacing $d$) and with the same "
                "coupling per unit length smeared along the line (multisource_bounce.py; grid spacing 0.5). Last column: range of the background field on "
                "the axis. The single-source axisymmetric code gives $\\Delta b=%.1f$ for one source at the bubble centre."
                % (MSJ["q"], MSJ["rs"], MSJ["single_source_2d_db"]), ["sources", "$d$", "$b$", "$\\Delta b$", "axis background"], rows, None))
    rows = [[texsci(r["m_dm_meV"], 1).join("$$"), texsci(r["q"], 1).join("$$"), texsci(r["r_s"], 1).join("$$"), texsci(r["screening"], 1).join("$$"),
             texsci(max(abs(r["ln_rate_ratio_screened"]), 1e-300), 1).join("$$"), texsci(max(abs(r["ln_rate_ratio_unscreened"]), 1e-300), 1).join("$$")]
            for r in SRJ["rows"][::3]]
    a(("table", "Table S14. Granular decay rate with self-screened particles (screened_rate.py): coupling, de Broglie radius, screening factor and "
                "$|\\ln(\\Gamma/\\Gamma_0)|$ with and without screening.", ["$m_{\\rm dm}$ (meV)", "$q$", "$r_s$", "$q_{\\rm eff}/q$", "screened",
                                                                              "unscreened"], rows, None))
    rows = [["%.2f" % r["one_plus_delta_today"], "%.3f" % r["a_tr_frozen"], "%.3f" % r["a_tr_evolving"], "%.3f" % r["one_plus_delta_at_tr"],
             "%+.3f" % r["shift"]] for r in EEV["rows"]]
    a(("table", "Table S15. Transition order with evolving density contrasts (env_evolution.py).",
       ["$1+\\delta$ today", "$a_{\\rm tr}$ frozen", "$a_{\\rm tr}$ evolving", "$1+\\delta$ at $a_{\\rm tr}$", "shift"], rows, None))
    rows = [["%.0f" % h["z"], "%.3g" % h["m_phi_meV"], "%.3g" % h["T_nu_meV"], "%.0f" % h["m_over_T"], texsci(h["gamma_over_H_per_y4"], 1).join("$$")]
            for h in CMB["history"]]
    a(("table", "Table S16. Early-universe history for the optional neutrino coupling: scalar mass, neutrino temperature (meV) and the ratio of the "
                "neutrino--neutrino scattering rate to the expansion rate per $y^4$ (contact estimate).",
       ["$z$", "$m_\\varphi$ (meV)", "$T_\\nu$ (meV)", "$m_\\varphi/T_\\nu$", "$\\Gamma/Hy^4$"], rows, None))
    rows = [[sv["kind"], k, "%.3g" % v["min"], "%.3g" % v["q10"], "%.3g" % v["median"], "%.3g" % v["q90"], "%.3g" % v["max"]]
            for sv in VEJ["surveys"] for k, v in sv["viable_physical_ranges"].items()]
    a(("table", "Table S17. Physical parameters of the draws that satisfy every condition (viability_extended.py).",
       ["priors", "parameter", "min", "10%", "median", "90%", "max"], rows, None))
    a(("h1", "S2. Reproduction"))
    a(("p", "Environment of the published run: Python %s, NumPy %s, SciPy %s (Windows). From the code directory: vacuum_phase_map.py, vacuum_wall.py, "
            "interface_transport.py (about 25 minutes), bubble_gravity.py, bubble_to_cycle.py, bubble_handoff.py, viability_mc.py (about 15 minutes), "
            "bubble_entropy.py, interface_validation.py, granularity.py, granular_field.py, wave_wall_1d.py (about 3 minutes), region_sensitivity.py, "
            "viability_extended.py (about 5 minutes), env_evolution.py, nonpert_bounce.py (about 7 minutes; 'convergence' and 'screening' modes), "
            "multisource_bounce.py (about 20 minutes), screened_rate.py, neutrino_coupling.py, make_int_figures.py, make_gran_figure.py, make_resolved_figures.py, "
            "paper2_build.py. Input tables shared with "
            "the companion paper (bounce_tables.json, thermal_tables.json, survival.json, density_bounce.json, coupling.json, continuous_cycles.json, "
            "cycle_map.json) are included as data and documented in the README. Code: %s." % (platform.python_version(), np.__version__, scipy.__version__,
                                                                                           GITHUB)))
    S = [CE.block(b) for b in S]
    RN.render_docx(S, os.path.join(OUT, BASENAME + "_Supplement.docx"), FIGDIR)


def letter():
    d = Document()
    st = d.styles["Normal"]; st.font.name = "Times New Roman"; st.font.size = Pt(11); st.paragraph_format.space_after = Pt(8)
    for s in d.sections:
        s.left_margin = s.right_margin = Inches(1.1); s.top_margin = s.bottom_margin = Inches(1.0)

    def para(text, bold=False, after=None):
        p = d.add_paragraph(); r = p.add_run(CE.ce(text)); r.bold = bold
        if after is not None:
            p.paragraph_format.space_after = Pt(after)
    para("Leon Sandler", True, 0)
    for t in ("Independent Researcher", "Northbrook, Illinois, USA", "sandler.leon@gmail.com", "ORCID 0009-0007-4584-808X"):
        para(t, after=0)
    para("")
    para(datetime.date.today().strftime("%d %B %Y").lstrip("0"))
    para("The Editors", after=0); para("Physics of the Dark Universe", after=0); para("Elsevier")
    para("Dear Editors,")
    para("I submit the manuscript “%s” as a research article for consideration in Physics of the Dark Universe." % TITLE)
    para("The paper asks whether a declining dark-matter density can switch on the decay of a vacuum, and what the resulting vacuum walls do. For a scalar "
         "coupled to dark matter it maps where the held vacuum is stable and where it becomes metastable, and then asks a question that the literature on "
         "environment-dependent vacua has not had to face: when does the field feel the smooth mean density that such an argument assumes? For a vacuum at "
         "the dark-energy scale the field responds within about 0.07 mm, inside which the dark-matter energy is only a few percent of the field mass to the "
         "fourth power. Particle dark matter of any mass is then a set of point sources between which the vacuum decays at its zero-density rate: "
         "particle dark matter does not provide a smooth source on the response scale of the field; the coincidence of the dark-matter and dark-energy densities, which motivates such a trigger, is "
         "what makes it fail. Spatially resolved calculations, the static field around discrete sources and the decay rate with randomly placed "
         "sources, confirm this within a first-order treatment, and the mean field is restored only for wave dark matter between about 3 meV and "
         "0.25 eV, which a 1+1-dimensional simulation of a wall running through such dark matter supports.")
    para("For that dark matter the transition is environmental: it begins in the emptiest voids, the past is protected only if the mean coupling today is "
         "close to the degeneracy density, and galactic halos keep the old vacuum. The paper then computes the Coleman-De Luccia bounce and the wall, the "
         "transmission of fields and of the dark matter through the wall (which runs away within a millimeter and filters a fraction of order 1e-88), the "
         "Israel junction, the completion of the transition and the fate of the enclosed regions; an isolated Hubble-scale bubble develops a trapped "
         "surface. An optional neutrino coupling is shown to be invisible to terrestrial oscillation experiments, with its cosmological signatures given. A "
         "continuation through a bounce is shown not to be supported, because it would need a thermal relic.")
    para("Related submission. The homogeneous, classical version of the trigger and its continuation into cycles are the subject of a separate manuscript "
         "submitted to General Relativity and Gravitation%s. Both grew out of one earlier, broader draft that a review recommended dividing; the granularity "
         "bound reported here does not affect that paper, whose field has a cosmological Compton length. The two share the class of models but not their "
         "results, and I will provide the companion manuscript on request."
         % ((" (preprint https://doi.org/%s)" % COMP["paper"]) if COMP.get("paper") else ""))
    para("The analysis is fully reproducible: the code (Python, MIT license) and the data behind every figure and table are at %s and archived at "
         "https://doi.org/%s; the manuscript is archived at https://doi.org/%s, and it is built only if the %d validation checks pass."
         % (GITHUB, DOIS.get("code", "(to be assigned)"), DOIS.get("paper", "(to be assigned)"), VJ["n"]))
    para("A generative AI tool (Claude, Anthropic) assisted with the numerical analysis, code, literature checking and drafting; I reviewed all of it and take "
         "full responsibility for the content, as declared in the manuscript. The manuscript is original, has not been published elsewhere and is not under "
         "consideration at any other journal. I am the sole author, I have no competing interests, and the work received no external funding.")
    para("Thank you for considering this submission.")
    para("Sincerely,")
    para("Leon Sandler", True)
    d.save(os.path.join(OUT, "PDU_Cover_Letter.docx"))


if __name__ == "__main__":
    for h in HIGHLIGHTS:
        assert len(h) <= 85, (len(h), h)
    os.makedirs(OUT, exist_ok=True)
    Bd = B[:3] + [("h1", "Highlights"), ("bullets", HIGHLIGHTS)] + B[3:]
    RN.render_docx(Bd, os.path.join(OUT, BASENAME + ".docx"), FIGDIR)
    RN.render_tex_els(B, os.path.join(OUT, BASENAME + ".tex"), "Physics of the Dark Universe", HIGHLIGHTS, "sandler.leon@gmail.com",
                      "organization={Independent Researcher},city={Northbrook},state={Illinois},country={USA}")
    supplement()
    letter()
    words = sum(len(b[1].split()) for b in B if b[0] in ("p", "abstract")) + sum(len(x.split()) for b in B if b[0] == "bullets" for x in b[1])
    print("Paper 2: %d words, %d figures, %d tables, %d equations, %d references" % (words, sum(1 for b in B if b[0] == "fig"),
                                                                                    sum(1 for b in B if b[0] == "table"), sum(1 for b in B if b[0] == "eq"), len(order)))
