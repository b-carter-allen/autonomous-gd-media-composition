#!/usr/bin/env python3
"""Generate simulated experiment data for dashboard testing.

Simulates 8 iterations of gradient descent converging to an optimal
media composition. Creates state.json and iteration_N/ directories
with iteration_log.json and transfer_array.json files.

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
    MIN_NOVEL_BIO_UL,
    MIN_SUPPLEMENT_UL,
    MAX_SUPPLEMENT_UL,
    PERTURBATION_COLS,
    REAGENT_VOLUME_UL,
    REAGENT_WELLS,
    ROWS,
    SEED_TRANSFER_VOLUME,
    SUPPLEMENT_NAMES,
    WELL_VOLUME_UL,
    apply_constraints,
    compute_novel_bio,
    compute_tip_consumption,
    generate_transfer_array,
    get_seed_params,
    gradient_step,
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
    "NaCl_2M": 5,
}

# Growth function parameters
BASELINE_OD = 0.08       # Starting OD (cells + media background)
MAX_GROWTH = 0.55         # Maximum achievable delta OD at optimal
NEG_CONTROL_GROWTH = 0.005  # Neg control drift (no cells)
POS_CONTROL_GROWTH = 0.30   # Growth in plain Novel_Bio + cells
NOISE_STD = 0.008         # Measurement noise (stdev) — low for clean convergence

PLATE_BARCODE = "SIM_PLATE_001"


def growth_function(composition: dict) -> float:
    """Simulate growth (delta OD) as a function of media composition.

    Uses a Gaussian-like response surface centered on OPTIMAL.
    Widths are broad enough that the sign-based gradient consistently
    points toward the optimum, even from 10+ uL away.

    Note: the sign-based stepper (step = alpha * delta * sign(gradient))
    can only reach compositions on a grid (multiples of 5 uL at alpha=0.5).
    The algorithm converges to the nearest reachable point to OPTIMAL.
    """
    # Width parameters — broad enough for clear gradient direction
    widths = {
        "Glucose_100mg_mL": 20,  # Broad
        "MOPS_1M": 18,           # Broad
        "NaCl_2M": 15,           # Slightly tighter (salt stress)
    }

    # Compute distance from optimal (weighted by sensitivity)
    exponent = 0.0
    for name in SUPPLEMENT_NAMES:
        diff = composition.get(name, 0) - OPTIMAL[name]
        exponent += (diff ** 2) / (2 * widths[name] ** 2)

    # Growth peaks at OPTIMAL, falls off as a Gaussian
    growth = MAX_GROWTH * math.exp(-exponent)

    # Floor: even bad compositions show some growth with cells
    growth = max(growth, 0.05)

    return growth


def simulate_od(composition: dict, noise: bool = True) -> float:
    """Simulate a delta OD measurement for a given composition."""
    growth = growth_function(composition)
    if noise:
        growth += random.gauss(0, NOISE_STD)
    return round(growth, 4)


def run_simulation():
    """Run a full simulated experiment and write data files."""
    barcode = sys.argv[1] if len(sys.argv) > 1 else PLATE_BARCODE
    data_dir = Path(__file__).parent.parent / "data" / barcode
    data_dir.mkdir(parents=True, exist_ok=True)

    print(f"Simulating {MAX_ITERATIONS} iterations for plate {barcode}")
    print(f"  Optimal: {OPTIMAL}")
    print(f"  Starting: {INITIAL_COMPOSITION}")
    print(f"  Output: {data_dir}")
    print()

    random.seed(42)  # Reproducible results

    state = {
        "current_iteration": 0,
        "current_composition": deepcopy(INITIAL_COMPOSITION),
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
        center = state["current_composition"]
        alpha = state["alpha"]
        perturbation_delta = max(1, int(alpha * DELTA_UL))

        iter_time = base_time + timedelta(hours=2 * (iteration - 1))
        iter_dir = data_dir / f"iteration_{iteration}"
        iter_dir.mkdir(parents=True, exist_ok=True)

        print(f"{'='*60}")
        print(f"  ITERATION {iteration} / {MAX_ITERATIONS}  (row {row_letter})")
        print(f"{'='*60}")

        # Generate transfer array
        transfer_array = generate_transfer_array(center, row_letter, delta=perturbation_delta)
        (iter_dir / "transfer_array.json").write_text(json.dumps(transfer_array, indent=2))

        # Write workflow definition
        seed_params = get_seed_params(iteration)
        workflow_path = write_workflow_definition(
            transfer_array, row_letter, iteration, seed_params
        )

        # Simulate OD measurements
        novel_bio_center = compute_novel_bio(center)

        # Neg control (col 2): no cells
        neg_delta = round(random.gauss(NEG_CONTROL_GROWTH, 0.003), 4)

        # Pos control (col 3): 180 uL Novel_Bio + cells
        pos_delta = simulate_od({"Glucose_100mg_mL": 0, "MOPS_1M": 0, "NaCl_2M": 0})
        pos_delta = round(POS_CONTROL_GROWTH + random.gauss(0, NOISE_STD), 4)

        # Center (col 4)
        center_delta = simulate_od(center)

        # Perturbation wells
        perturbed_ods = {}
        for col1, col2, supplement in PERTURBATION_COLS:
            perturbed_comp = make_perturbed(center, supplement, perturbation_delta)
            rep1 = simulate_od(perturbed_comp)
            rep2 = simulate_od(perturbed_comp)
            perturbed_ods[supplement] = [rep1, rep2]

        # Extra wells (same as center)
        extra1 = simulate_od(center)
        extra2 = simulate_od(center)

        # Build od_results
        od_results = {
            "neg_control_od": neg_delta,
            "control_od": pos_delta,
            "center_od": center_delta,
            "perturbed_ods": perturbed_ods,
            "extra_ods": [extra1, extra2],
            "abs_control_od": round(BASELINE_OD + pos_delta, 4),
            "abs_center_od": round(BASELINE_OD + center_delta, 4),
        }

        # Compute gradient
        gradient = {}
        for name in SUPPLEMENT_NAMES:
            mean_perturbed = sum(perturbed_ods[name]) / len(perturbed_ods[name])
            gradient[name] = (mean_perturbed - center_delta) / perturbation_delta

        # Adaptive alpha: halve if center OD decreased
        if state["prev_center_od"] is not None:
            if center_delta < state["prev_center_od"]:
                alpha = max(0.5, alpha / 2)
                print(f"  Alpha halved to {alpha} (center OD decreased)")

        # Step composition
        new_composition = gradient_step(center, gradient, alpha)

        # Print summary
        novel_bio = compute_novel_bio(center)
        supps = ", ".join(f"{k}={v}" for k, v in center.items())
        print(f"  Composition: Novel_Bio={novel_bio}, {supps}")
        print(f"  Center OD: {center_delta:+.4f}")
        print(f"  Gradient: {', '.join(f'{k}={v:+.4f}' for k, v in gradient.items())}")
        new_supps = ", ".join(f"{k}={v}" for k, v in new_composition.items())
        print(f"  Next: Novel_Bio={compute_novel_bio(new_composition)}, {new_supps}")
        print()

        # Write iteration log
        log = {
            "iteration": iteration,
            "timestamp": iter_time.isoformat(),
            "composition": {
                "Novel_Bio": novel_bio,
                **center,
            },
            "transfer_count": len(transfer_array),
            "alpha": alpha,
            "od_results": od_results,
            "gradient": gradient,
            "next_composition": {
                "Novel_Bio": compute_novel_bio(new_composition),
                **new_composition,
            },
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
        state["current_composition"] = new_composition
        state["alpha"] = alpha
        state["prev_center_od"] = center_delta
        if state["best_od"] is None or center_delta > state["best_od"]:
            state["best_od"] = center_delta

        state["history"].append({
            "iteration": iteration,
            "composition": {
                "Novel_Bio": novel_bio,
                **center,
            },
            "center_od": center_delta,
            "control_od": pos_delta,
            "neg_control_od": neg_delta,
            "extra_ods": [extra1, extra2],
            "abs_center_od": round(BASELINE_OD + center_delta, 4),
            "abs_control_od": round(BASELINE_OD + pos_delta, 4),
            "gradient": gradient,
            "alpha": alpha,
        })

    # Mark as converged (max iterations reached)
    state["converged"] = True

    # Save final state
    (data_dir / "state.json").write_text(json.dumps(state, indent=2, default=str))

    # Print summary table
    print(f"\n{'='*60}")
    print(f"  SIMULATION COMPLETE")
    print(f"{'='*60}")
    print(f"  Best center OD: {state['best_od']:.4f}")
    final = state["current_composition"]
    print(f"  Final composition: Novel_Bio={compute_novel_bio(final)}, "
          + ", ".join(f"{k}={v}" for k, v in final.items()))
    print(f"\n  Iteration history:")
    for h in state["history"]:
        c = h["composition"]
        supps = ", ".join(f"{name}={c.get(name, 0)}" for name in SUPPLEMENT_NAMES)
        print(f"    {h['iteration']}: NB={c['Novel_Bio']}, {supps} -> OD={h['center_od']:.4f} (a={h['alpha']})")

    print(f"\n  Data written to: {data_dir}")


if __name__ == "__main__":
    run_simulation()
