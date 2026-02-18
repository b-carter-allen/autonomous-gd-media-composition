# Experimental Design

This folder is where you define **what** to test before running the autonomous optimizer.

Fill in the sections below (or create your own files) to document your experimental design. The gradient descent daemon in `experimental_execution/` will handle the **how**.

## What to Define

### Hypothesis

What are you testing? What do you expect to find?

Example: *"Adding glucose to Novel_Bio media will increase V. natriegens growth rate, while NaCl and MgSO4 will have minimal or negative effects at the concentrations tested."*

### Reagents

Which supplements are you optimizing? The daemon currently supports 3 supplements (rows C-H on the plate, 2 replicate wells each). To change which reagents are used, update `SUPPLEMENT_NAMES` and `REAGENT_WELLS` in `gradient_descent.py`.

| Reagent | Well on Reagent Plate | Role |
|---------|----------------------|------|
| Glucose | A1 | Supplement 1 |
| NaCl | B1 | Supplement 2 |
| MgSO4 | C1 | Supplement 3 |
| Novel_Bio | D1 | Base media (fills remaining volume) |

### Starting Composition

What volumes (in uL) should the optimizer start with? Set these in `INITIAL_COMPOSITION` in `gradient_descent.py`.

Total well volume is 180 uL. The remainder after supplements is filled with Novel_Bio base media.

| Parameter | Value | Notes |
|-----------|-------|-------|
| Glucose | 20 uL | Starting volume |
| NaCl | 20 uL | Starting volume |
| MgSO4 | 20 uL | Starting volume |
| Novel_Bio | 120 uL | Auto-calculated (180 - sum of supplements) |

### Objective Function

What are you optimizing? The default is **delta OD600** (endpoint - baseline), which measures actual growth over the monitoring window rather than absolute optical density.

Alternatives to consider:
- Absolute OD600 at endpoint
- Growth rate (slope of OD600 over time)
- Max OD600 achieved during monitoring

### Constraints

| Parameter | Value | Why |
|-----------|-------|-----|
| `WELL_VOLUME_UL` | 180 uL | Fixed by plate format |
| `MIN_SUPPLEMENT_UL` | 1 uL | Minimum pipettable volume (Flex) |
| `MAX_SUPPLEMENT_UL` | 90 uL | Prevents any single supplement from dominating |
| `MIN_NOVEL_BIO_UL` | 90 uL | Ensures cells have base nutrients |

### Algorithm Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `DELTA_UL` | 10 | Base perturbation size (uL) |
| `ALPHA` | 1.0 | Learning rate (decays by half when growth decreases) |
| `MAX_ITERATIONS` | 8 | Limited by plate columns (1 seed column + 8 experiment columns) |
| `MONITORING_INTERVAL_MINUTES` | 5 | Time between OD600 readings |
| `MONITORING_READINGS` | 18 | Number of readings (18 x 5 min = 90 min) |

All parameters are defined at the top of `experimental_execution/gradient_descent.py`.
