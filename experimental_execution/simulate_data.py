#!/usr/bin/env python3
"""Generate simulated experiment data for dashboard testing.

Prescribes a convergence trajectory toward the optimal composition
and generates consistent OD values, transfer arrays, and workflow
definitions for each iteration.

Usage:
    python simulate_data.py [plate_barcode]
"""

import json
import math
import random
import sys
from copy import deepcopy
from datetime import datetime, timedelta
from pathlib import Path

# Import constants and helpers from the daemon
from gradient_descent import (
    ALPHA,
    DELTA_UL,
    INITIAL_COMPOSITION,
    MAX_ITERATIONS,
    PERTURBATION_COLS,
    REAGENT_VOLUME_UL,
    ROWS,
    SUPPLEMENT_NAMES,
    WELL_VOLUME_UL,
    compute_novel_bio,
    generate_transfer_array,
    get_seed_params,
    make_perturbed,
    write_workflow_definition,
)

# =============================================================================
# SIMULATION PARAMETERS
# =============================================================================

# Optimal composition the optimizer should converge toward
OPTIMAL = {
    "Glucose_100mg_mL": 12,
    "MOPS_1M": 7,
    "DiH2O": 5,
}

# Growth function parameters
BASELINE_OD = 0.08        # Starting OD (cells + media background)
MAX_GROWTH = 0.55         # Maximum achievable delta OD at optimal
POS_CONTROL_GROWTH = 0.30 # Growth in plain Novel_Bio + cells
NOISE_STD = 0.008         # Measurement noise (stdev)

PLATE_BARCODE = "SIM_PLATE_001"

# Prescribed trajectory: compositions the optimizer will visit.
# Designed to show realistic convergence behavior:
#   1-2: Exploration phase (big steps, alpha=1.0)
#   3-4: Recovery after overshoot (alpha halves to 0.5)
#   5-6: Narrowing in (alpha halves to 0.25, step=2)
#   7-8: Fine-tuning (alpha at 0.125, step=1)
TRAJECTORY = [
    # iter 1: start
    {"Glucose_100mg_mL": 20, "MOPS_1M": 20, "DiH2O": 20},
    # iter 2: big step toward optimal (all gradients negative, step=10)
    {"Glucose_100mg_mL": 10, "MOPS_1M": 10, "DiH2O": 10},
    # iter 3: overshoot recovery (alpha halved to 0.5, step=5)
    {"Glucose_100mg_mL": 15, "MOPS_1M": 5,  "DiH2O": 5},
    # iter 4: alpha stays 0.5, step=5 — Glucose overshoots slightly
    {"Glucose_100mg_mL": 10, "MOPS_1M": 10, "DiH2O": 5},
    # iter 5: alpha halves to 0.25, step=2 — narrowing
    {"Glucose_100mg_mL": 12, "MOPS_1M": 8,  "DiH2O": 5},
    # iter 6: step=2, fine adjustment
    {"Glucose_100mg_mL": 12, "MOPS_1M": 6,  "DiH2O": 5},
    # iter 7: alpha halves to 0.125, step=1 — near optimal
    {"Glucose_100mg_mL": 12, "MOPS_1M": 7,  "DiH2O": 5},
    # iter 8: converged — same as optimal
    {"Glucose_100mg_mL": 12, "MOPS_1M": 7,  "DiH2O": 5},
]

# Alpha schedule (prescribed to match the trajectory logic)
ALPHA_SCHEDULE = [1.0, 1.0, 0.5, 0.5, 0.25, 0.25, 0.125, 0.125]


def growth_function(composition: dict) -> float:
    """Simulate growth (delta OD) as a function of media composition.

    Asymmetric Gaussian: steep penalty below optimal (deficiency),
    gentle above (excess tolerated). Biologically realistic.
    """
    widths = {
        "Glucose_100mg_mL": (8, 25),
        "MOPS_1M": (7, 20),
        "DiH2O": (5, 15),
    }

    exponent = 0.0
    for name in SUPPLEMENT_NAMES:
        x = composition.get(name, 0)
        opt = OPTIMAL[name]
        diff = x - opt
        w_below, w_above = widths[name]
        w = w_below if diff < 0 else w_above
        exponent += (diff ** 2) / (2 * w ** 2)

    growth = MAX_GROWTH * math.exp(-exponent)
    return max(growth, 0.05)


def simulate_od(composition: dict) -> float:
    """Simulate a delta OD measurement with noise."""
    growth = growth_function(composition)
    growth += random.gauss(0, NOISE_STD)
    return round(growth, 4)


def run_simulation():
    """Generate simulated experiment data following a prescribed trajectory."""
    barcode = sys.argv[1] if len(sys.argv) > 1 else PLATE_BARCODE
    data_dir = Path(__file__).parent.parent / "data" / barcode
    data_dir.mkdir(parents=True, exist_ok=True)

    print(f"Simulating {MAX_ITERATIONS} iterations for plate {barcode}")
    print(f"  Optimal: {OPTIMAL}")
    print(f"  Starting: {TRAJECTORY[0]}")
    print(f"  Output: {data_dir}")
    print()

    random.seed(42)

    state = {
        "current_iteration": 0,
        "current_composition": deepcopy(TRAJECTORY[0]),
        "alpha": ALPHA,
        "best_od": None,
        "prev_center_od": None,
        "no_improvement_count": 0,
        "converged": False,
        "history": [],
    }

    base_time = datetime(2026, 2, 10, 9, 0, 0)

    for iteration in range(1, MAX_ITERATIONS + 1):
        row_letter = ROWS[iteration - 1]
        center = TRAJECTORY[iteration - 1]
        alpha = ALPHA_SCHEDULE[iteration - 1]
        perturbation_delta = max(1, int(alpha * DELTA_UL))

        # Next composition (for logging)
        next_comp = TRAJECTORY[iteration] if iteration < MAX_ITERATIONS else TRAJECTORY[-1]

        iter_time = base_time + timedelta(hours=2 * (iteration - 1))
        iter_dir = data_dir / f"iteration_{iteration}"
        iter_dir.mkdir(parents=True, exist_ok=True)

        print(f"{'='*60}")
        print(f"  ITERATION {iteration} / {MAX_ITERATIONS}  (row {row_letter}, alpha={alpha})")
        print(f"{'='*60}")

        # Generate real transfer array and workflow definition
        transfer_array = generate_transfer_array(center, row_letter, delta=perturbation_delta)
        (iter_dir / "transfer_array.json").write_text(json.dumps(transfer_array, indent=2))

        seed_params = get_seed_params(iteration)
        write_workflow_definition(transfer_array, row_letter, iteration, seed_params)

        # Simulate OD measurements
        neg_delta = round(random.gauss(0.005, 0.003), 4)
        pos_delta = round(POS_CONTROL_GROWTH + random.gauss(0, NOISE_STD), 4)
        center_delta = simulate_od(center)

        # Perturbation wells
        perturbed_ods = {}
        for _, _, supplement in PERTURBATION_COLS:
            perturbed_comp = make_perturbed(center, supplement, perturbation_delta)
            rep1 = simulate_od(perturbed_comp)
            rep2 = simulate_od(perturbed_comp)
            perturbed_ods[supplement] = [rep1, rep2]

        extra1 = simulate_od(center)
        extra2 = simulate_od(center)

        od_results = {
            "neg_control_od": neg_delta,
            "control_od": pos_delta,
            "center_od": center_delta,
            "perturbed_ods": perturbed_ods,
            "extra_ods": [extra1, extra2],
            "abs_control_od": round(BASELINE_OD + pos_delta, 4),
            "abs_center_od": round(BASELINE_OD + center_delta, 4),
        }

        # Compute gradient (for logging — doesn't drive the trajectory)
        gradient = {}
        for name in SUPPLEMENT_NAMES:
            mean_pert = sum(perturbed_ods[name]) / len(perturbed_ods[name])
            gradient[name] = round((mean_pert - center_delta) / perturbation_delta, 6)

        # Print summary
        novel_bio = compute_novel_bio(center)
        supps = ", ".join(f"{k}={v}" for k, v in center.items())
        print(f"  Composition: Novel_Bio={novel_bio}, {supps}")
        print(f"  Center OD: {center_delta:+.4f}")
        print(f"  Gradient: {', '.join(f'{k}={v:+.6f}' for k, v in gradient.items())}")
        next_supps = ", ".join(f"{k}={v}" for k, v in next_comp.items())
        print(f"  Next: Novel_Bio={compute_novel_bio(next_comp)}, {next_supps}")
        print()

        # Write iteration log
        log = {
            "iteration": iteration,
            "timestamp": iter_time.isoformat(),
            "composition": {"Novel_Bio": novel_bio, **center},
            "transfer_count": len(transfer_array),
            "alpha": alpha,
            "od_results": od_results,
            "gradient": gradient,
            "next_composition": {"Novel_Bio": compute_novel_bio(next_comp), **next_comp},
        }
        (iter_dir / "iteration_log.json").write_text(json.dumps(log, indent=2))

        # Write simulated workflow IDs
        (iter_dir / "workflow_ids.json").write_text(json.dumps({
            "workflow_definition_id": 100 + iteration,
            "workflow_instance_uuid": f"sim-{barcode}-iter{iteration}",
        }, indent=2))

        # Check improvement
        any_improvement = any(
            sum(perturbed_ods[name]) / 2 > center_delta
            for name in SUPPLEMENT_NAMES
        )
        if not any_improvement:
            state["no_improvement_count"] += 1
        else:
            state["no_improvement_count"] = 0

        # Update state
        state["current_iteration"] = iteration
        state["current_composition"] = deepcopy(next_comp)
        state["alpha"] = alpha
        state["prev_center_od"] = center_delta
        if state["best_od"] is None or center_delta > state["best_od"]:
            state["best_od"] = center_delta

        state["history"].append({
            "iteration": iteration,
            "composition": {"Novel_Bio": novel_bio, **center},
            "center_od": center_delta,
            "control_od": pos_delta,
            "neg_control_od": neg_delta,
            "extra_ods": [extra1, extra2],
            "abs_center_od": round(BASELINE_OD + center_delta, 4),
            "abs_control_od": round(BASELINE_OD + pos_delta, 4),
            "gradient": gradient,
            "alpha": alpha,
        })

    state["converged"] = True
    (data_dir / "state.json").write_text(json.dumps(state, indent=2, default=str))

    # Summary
    print(f"\n{'='*60}")
    print(f"  SIMULATION COMPLETE")
    print(f"{'='*60}")
    print(f"  Best center OD: {state['best_od']:.4f}")
    final = state["current_composition"]
    print(f"  Final: Novel_Bio={compute_novel_bio(final)}, "
          + ", ".join(f"{k}={v}" for k, v in final.items()))
    print(f"\n  Iteration history:")
    for h in state["history"]:
        c = h["composition"]
        supps = ", ".join(f"{name}={c.get(name, 0)}" for name in SUPPLEMENT_NAMES)
        print(f"    {h['iteration']}: NB={c['Novel_Bio']}, {supps} "
              f"-> OD={h['center_od']:.4f} (a={h['alpha']})")

    print(f"\n  Data written to: {data_dir}")


if __name__ == "__main__":
    run_simulation()
