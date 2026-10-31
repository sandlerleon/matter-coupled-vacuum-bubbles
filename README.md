# Density-Dependent Vacuum Bubble Nucleation in Matter-Coupled Dark-Sector Cosmology

**Granularity of the Source, Wall Dynamics, Transmission, and Gravitational Viability**

Leon Sandler, Independent Researcher — sandler.leon@gmail.com — ORCID 0009-0007-4584-808X

Code, data and manuscript for a paper prepared for **Physics of the Dark Universe** (Elsevier). Manuscript:
`manuscript/PDU_Vacuum_Bubbles.{docx,tex,pdf}` (elsarticle), supplement `manuscript/PDU_Vacuum_Bubbles_Supplement.docx`, cover letter
`manuscript/PDU_Cover_Letter.docx`.

## What the paper shows

A scalar coupled to dark matter through A(φ) = 1 + (φ − f)²/2M², centred on a false vacuum, is held there by the matter term while the dark-matter
density is high. In the mean-field description, dilution turns the held vacuum from the true into a metastable one at θ_deg = η/2, and
Coleman–De Luccia tunneling switches on below it.

* **Granularity bound.** For a vacuum at the dark-energy scale the field's Compton length is about 0.07 mm, and the dark-matter energy inside a Compton
  volume is only ≈ 0.05 m_φ⁴. Particle dark matter of any mass is a set of point sources between which the field sees zero density; the vacuum would
  decay at its zero-density rate (e^86 expected bubbles in our past light cone). The mean field is restored only for **wave dark matter heavier than the
  field and still classical, ≈ 3 meV – 0.25 eV**. The bound follows from the coincidence of the dark-matter and dark-energy densities and holds over
  all tilts up to η_c (the ratio ρ/m_φ⁴ exceeds 6 only at vanishing density).
* **Environment.** With such dark matter the decay is local: it begins in the emptiest voids; protecting the past requires a mean coupling today close
  to the degeneracy density; walls stall where the local density reaches the degeneracy value, and galactic haloes keep the held vacuum.
* **Walls.** About 0.2 mm thick; they run away, become transparent to the dark matter (filtered fraction 2×10⁻⁸⁸), carry the released energy, and
  are gravitationally negligible (Israel junction). Completion: β/H = 351; void by void within a few hundredths of a Hubble time.
* **Fate.** Transitioned regions recollapse in their negative vacuum; an isolated Hubble-scale bubble is trapped before its restoration and ends as a
  black hole. A continuation through a bounce would need a thermal relic and is **not** supported by the dark matter the trigger requires.

## Scripts and results

| Result | Paper | Scripts |
|---|---|---|
| Phase structure, degeneracy θ = η/2, symmetron contrast, coupling check | 2 | `vacuum_phase_map.py` |
| Granularity bound, wave-dark-matter window, structural ratio and approach to η_c, voids, void-by-void completion | 3, 4.4 | `granularity.py` |
| O(4) bounce, wall tension and thickness, thin-wall comparison, nucleation history | 4 | `vacuum_wall.py`, `cdl_lib.py`, `thermal.py` |
| Transmission through the wall (Pöschl–Teller validated), thermal pressure, dark-matter runaway | 5 | `interface_transport.py` |
| Israel junction, region evolution, controls, viability map | 6, 7, 8 | `bubble_gravity.py`, `bubble_to_cycle.py` |
| Percolation, bounds on the approximations, isolated bubble, variants, independent integration | 6, 7 | `bubble_handoff.py` |
| Energy and entropy closure | 7.3 | `bubble_entropy.py` |
| Parameter survey | 8 | `viability_mc.py` |
| 40 validation checks (the manuscript is built only if they pass) | 8, S1 | `interface_validation.py` |

## Reproduction

Python 3.11 with the packages in `requirements.txt`; pdflatex with `elsarticle.cls` for the PDF. From `code/`:

```
python vacuum_phase_map.py && python vacuum_wall.py && python interface_transport.py      # ~25 min for the transport
python bubble_gravity.py && python bubble_to_cycle.py && python bubble_handoff.py && python bubble_entropy.py
python viability_mc.py && python interface_validation.py && python granularity.py
python make_int_figures.py && python make_gran_figure.py && python paper2_build.py
cd ../manuscript && pdflatex PDU_Vacuum_Bubbles.tex && pdflatex PDU_Vacuum_Bubbles.tex
```

Input tables shared with the companion paper (computed by scripts in `sandlerleon/late-time-vacuum-decay`, tag v1.8.2) are included as data:
`bounce_tables.json`, `thermal_tables.json`, `survival.json`, `density_bounce.json`, `coupling.json`, `continuous_cycles.json`, `cycle_map.json`.

## Companion paper

The homogeneous (symmetron) version of the trigger, whose field has a cosmological Compton length and is therefore not affected by the granularity
bound, and its continuation into cycles are treated in "Matter-Dilution-Triggered Recollapse and Conditional Nonsingular Cycles in a Metastable
Dark-Energy Model" (The European Physical Journal C; code `sandlerleon/late-time-vacuum-decay`). The two papers share the class of models but not
their results.

## License

MIT (code and data); the manuscript text is © the author.
