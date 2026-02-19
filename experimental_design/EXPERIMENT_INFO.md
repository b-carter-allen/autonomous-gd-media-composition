# Experiment: GD Media Optimization — Run 1

**Date**: 2026-02-12
**Organism**: *Vibrio natriegens*
**Base media**: Novel_Bio
**Plate format**: 96-well, row-wise iterations (A-H)

## Hypothesis

Novel_Bio base media may be oversalted for optimal *V. natriegens* growth. By supplementing with Glucose (carbon source), MOPS (pH buffer), and DiH2O (to dilute excess salt), the gradient descent optimizer should converge on an improved media composition that increases growth rate as measured by OD600.

We expect:
- **Glucose**: Moderate supplementation (10-20 uL of 100 mg/mL stock) will provide additional carbon and improve growth
- **MOPS**: Small amounts (5-15 uL of 1 M stock) will stabilize pH during growth
- **DiH2O**: Small amounts (5-10 uL) will dilute the salt concentration to an optimal level

## Reagents

| Reagent | Well on Reagent Plate | Volume | Concentration | Role |
|---------|----------------------|--------|---------------|------|
| Glucose | A1 | 5 mL | 100 mg/mL | Carbon source |
| MOPS | B1 | 5 mL | 1 M | pH buffer |
| DiH2O | C1 | 5 mL | Pure (deionized water) | Salt dilution |
| Novel_Bio | D1 | 9 mL | Undiluted | Base media (fills remaining volume) |
| Novel_Bio | D2 | 9 mL | Undiluted | Base media (overflow from D1) |
| NM+Cells | A2 | 5 mL | — | Seed culture for warmup (refrigerated) |

## Starting Composition

| Component | Volume (uL) | Notes |
|-----------|-------------|-------|
| Glucose 100 mg/mL | 20 | Starting volume |
| MOPS 1 M | 20 | Starting volume |
| DiH2O | 20 | Starting volume |
| Novel_Bio | 120 | Auto-calculated (180 - 60) |
| Cells | 20 | Seeded from column 1 |
| **Total** | **200** | |

## Objective Function

**Delta OD600** (endpoint - baseline): measures actual growth over the 90-minute monitoring window, normalizing for varying initial cell densities across wells.

## Constraints

| Parameter | Value | Why |
|-----------|-------|-----|
| `REAGENT_VOLUME_UL` | 180 uL | Reagent mix volume (before cells) |
| `SEED_TRANSFER_VOLUME` | 20 uL | Cells added from seed well |
| `WELL_VOLUME_UL` | 200 uL | Total volume (reagent + cells) |
| `MIN_SUPPLEMENT_UL` | 5 uL | Minimum pipettable volume (Flex accuracy) |
| `MAX_SUPPLEMENT_UL` | 90 uL | Prevents any single supplement from dominating |
| `MIN_NOVEL_BIO_UL` | 90 uL | Ensures cells have base nutrients |

## Algorithm Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `DELTA_UL` | 10 | Base perturbation size (uL) |
| `ALPHA` | 1.0 | Initial learning rate (decays to floor of 0.1) |
| `MAX_ITERATIONS` | 8 | One iteration per row (A-H) |
| `MONITORING_INTERVAL_MINUTES` | 5 | Time between OD600 readings |
| `MONITORING_READINGS` | 18 | Number of readings (18 x 5 min = 90 min) |

Step sizes narrow as alpha decays: 10 → 5 → 2 → 1 uL, allowing the optimizer to converge on the exact optimal volumes.

## Plate Layout (per row)

| Column | Role | Composition |
|--------|------|-------------|
| 1 | Seed well | NM+Cells (pre-loaded, 220 uL) |
| 2 | Negative control | 200 uL Novel_Bio, no cells |
| 3 | Positive control | 180 uL Novel_Bio + 20 uL cells |
| 4 | Center | Current best recipe + 20 uL cells |
| 5-6 | +Glucose | Center + delta Glucose (2 reps) + cells |
| 7-8 | +MOPS | Center + delta MOPS (2 reps) + cells |
| 9-10 | +DiH2O | Center + delta DiH2O (2 reps) + cells |
| 11-12 | Extra | Center duplicate (2 reps) + cells |

## Notes

- The "GD Iteration Combined" routine on the workcell handles all liquid handling in one step: reagent transfers, on-plate seeding from column 1, and NM+Cells warmup for the next row's seed well.
- Reagent plate is a 24-well deep well plate (type: "AGD Stock Plate", modeled as 96-well in the system). Novel_Bio is split across D1 and D2 (9 mL each, 18 mL total) to ensure enough volume for all 8 iterations. Supplements have 5 mL each.
- NM+Cells in A2 is refrigerated until needed for seed well warmup.
- DiH2O was chosen instead of NaCl supplementation because the working hypothesis is that Novel_Bio is already oversalted — adding water achieves the same gradient direction while diluting rather than adding more salt.
