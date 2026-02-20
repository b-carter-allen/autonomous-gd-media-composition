#!/usr/bin/env python3
"""Gradient descent media optimization for V. natriegens growth.

Supplements Novel_Bio with 3 configurable reagents and optimizes
concentrations via gradient descent to maximize OD600 growth after 1.5 hours.

Each iteration uses 1 row (12 wells), with column 1 as the seed well:
  Col 1:    Seed well (NM+Cells stock, pre-loaded)
  Col 2:    Negative control (200 uL Novel_Bio, no cells)
  Col 3:    Positive control (180 uL Novel_Bio + 20 uL cells)
  Col 4:    Center point (current best composition + cells)
  Col 5-6:  +delta Supplement 1 (2 replicates + cells)
  Col 7-8:  +delta Supplement 2 (2 replicates + cells)
  Col 9-10: +delta Supplement 3 (2 replicates + cells)
  Col 11-12: Extra wells (center duplicates + cells)

Autonomous daemon loop:
  1. Generate transfer array from current composition
  2. Write workflow definition with iteration-specific constants
  3. Upload + register workflow via MCP
  4. Instantiate workflow via MCP (auto-approved)
  5. Poll for completion via MCP
  6. Fetch OD600 results from datasets REST API (baseline + endpoint)
  7. Compute growth delta (endpoint - baseline) and gradient
  8. Repeat for up to 8 iterations (one per row, A-H)

Objective function: delta OD600 (growth = endpoint - baseline), not absolute OD.
This normalizes for varying initial cell densities across wells.

Usage:
  python gradient_descent.py run <plate_barcode>     # Start from iteration 1
  python gradient_descent.py resume <plate_barcode>  # Resume from last state
  python gradient_descent.py dry-run <plate_barcode> # Preview next iteration
  python gradient_descent.py status                  # Show current state
  python gradient_descent.py reset                   # Reset state
"""

import json
import os
import re
import sys
import time
from copy import deepcopy
from datetime import datetime
from pathlib import Path

import requests
from dotenv import load_dotenv

load_dotenv()

# =============================================================================
# CONFIGURATION
# =============================================================================

# Workcell connection
# Set WORKCELL_BASE_URL to override the full base URL (e.g. for Tailscale/HTTPS).
# Otherwise falls back to http://{WORKCELL_HOST}:{WORKCELL_PORT}.
_host = os.getenv("WORKCELL_HOST", "192.168.68.55")
_port = int(os.getenv("WORKCELL_PORT", "8080"))
WORKCELL_API_BASE = os.getenv("WORKCELL_BASE_URL", f"http://{_host}:{_port}")

# API request headers (ClientIdentifierMiddleware requires desktop-frontend)
_auth_header = os.getenv("WORKCELL_AUTH_HEADER", "")
API_HEADERS = {
    "Content-Type": "application/json",
    "X-Monomer-Client": "desktop-frontend",
    **( {"Authorization": _auth_header} if _auth_header else {} ),
}

# Backend MCP connection (for fetching absorbance observations)
BACKEND_MCP_URL = os.getenv("BACKEND_MCP_URL", "https://backend-staging.monomerbio.com/mcp")
_backend_auth = os.getenv("BACKEND_AUTH_HEADER", _auth_header)
BACKEND_MCP_HEADERS = {
    "Content-Type": "application/json",
    "Accept": "application/json, text/event-stream",
    **( {"Authorization": _backend_auth} if _backend_auth else {} ),
}

# Reagent plate well map (24-well deep well, modeled as 96-well in the system)
# UPDATE THESE with your chosen supplements before running.
REAGENT_WELLS = {
    "Glucose_100mg_mL": "A1",
    "MOPS_1M": "B1",
    "DiH2O": "C1",
    "Novel_Bio": "D1",  # Starting well; switches to D2, D3... as volume is consumed
}

# Novel Bio wells used sequentially as each fills. Fill each to NOVEL_BIO_WELL_CAPACITY_UL.
# The 24-well deep well plate holds 10.4 mL/well; 8000 uL is a safe fill target (~77%).
NOVEL_BIO_WELLS = ["D1", "D2", "D3", "D4"]
NOVEL_BIO_WELL_CAPACITY_UL = 8000  # uL of Novel Bio to load per well

# NM+Cells well on the reagent plate (24-well deep well, well A2)
NM_CELLS_REAGENT_WELL = "A2"
NM_CELLS_REAGENT_NAME = "NM+Cells"
NM_CELLS_FILL_UL = 2000  # Fill target: covers 7-8 rounds × 220 uL + margin

# Supplement fill target per well (conservative upper bound for one experiment)
SUPPLEMENT_FILL_UL = 3000

# Plate type tag used when registering AGD reagent plates
AGD_PLATE_TYPE = "AGD Stock Plate"

# Poll interval when waiting for the reagent plate to be checked in
REAGENT_PLATE_POLL_INTERVAL = 30   # seconds
REAGENT_PLATE_POLL_TIMEOUT = 1800  # 30 minutes


def get_novel_bio_well(cumulative_ul: float) -> str:
    """Return which Novel Bio well to use based on cumulative volume consumed so far."""
    well_index = int(cumulative_ul // NOVEL_BIO_WELL_CAPACITY_UL)
    if well_index >= len(NOVEL_BIO_WELLS):
        raise RuntimeError(
            f"Novel Bio exhausted: {cumulative_ul:.0f} uL used across "
            f"{len(NOVEL_BIO_WELLS)} wells x {NOVEL_BIO_WELL_CAPACITY_UL} uL = "
            f"{len(NOVEL_BIO_WELLS) * NOVEL_BIO_WELL_CAPACITY_UL} uL total capacity."
        )
    return NOVEL_BIO_WELLS[well_index]

# Reagent names in order (matches cols 5/6, 7/8, 9/10)
SUPPLEMENT_NAMES = ["Glucose_100mg_mL", "MOPS_1M", "DiH2O"]

# Plate layout: column -> purpose (row-wise iteration)
# Cols 2 and 11 are intentionally empty (buffer wells) to isolate:
#   - seed well (col 1) from experiments (cols 3-10)
#   - experiments (cols 3-10) from negative control (col 12)
COL_LABELS = {
    1: "seed",
    2: "empty_buffer",
    3: "pos_control",
    4: "center",
    5: "glucose_rep1",
    6: "glucose_rep2",
    7: "mops_rep1",
    8: "mops_rep2",
    9: "diH2O_rep1",
    10: "diH2O_rep2",
    11: "empty_buffer",
    12: "neg_control",
}

# Perturbation column pairs: (col1, col2, supplement_name)
PERTURBATION_COLS = [
    (5, 6, "Glucose_100mg_mL"),
    (7, 8, "MOPS_1M"),
    (9, 10, "DiH2O"),
]

# Volumes
REAGENT_VOLUME_UL = 180   # Reagent mix volume per well (before cells)
SEED_TRANSFER_VOLUME = 20  # Cells added from seed well
WELL_VOLUME_UL = 200       # Total volume (reagent + cells)

# Constraints
MIN_SUPPLEMENT_UL = 5   # Minimum if included (can be 0)
MAX_SUPPLEMENT_UL = 90  # Maximum per supplement
MIN_NOVEL_BIO_UL = 90   # Minimum Novel_Bio volume

# Algorithm parameters
DELTA_UL = 10            # Perturbation step size
ALPHA = 1.0              # Learning rate (multiplied by delta for step size)
MAX_ITERATIONS = 8       # Max iterations (one per row, A-H)
CONVERGENCE_ROUNDS = 2   # Stop if no improvement for this many consecutive rounds

# Starting point (iteration 1 center)
# UPDATE THESE to match your SUPPLEMENT_NAMES.
INITIAL_COMPOSITION = {
    "Glucose_100mg_mL": 20,
    "MOPS_1M": 20,
    "DiH2O": 20,
}

# Monitoring
MONITORING_INTERVAL_MINUTES = 5
MONITORING_READINGS = 18  # 1.5 hours

# Paths
# DATA_DIR is set dynamically per experiment in main() based on plate barcode.
# Default here for status/reset commands that run before barcode is known.
DATA_DIR = Path(__file__).parent.parent / "data" / "gradient_descent"
WORKFLOW_TEMPLATE_PATH = Path(__file__).parent / "workflow_template.py"

# Iteration rows (one row per iteration)
ROWS = ["A", "B", "C", "D", "E", "F", "G", "H"]

# Daemon config
DAEMON_POLL_INTERVAL = 30  # seconds between workflow status checks
WORKFLOW_TIMEOUT_MINUTES = 180  # max wait per workflow (3 hours)

# =============================================================================
# STATE MANAGEMENT
# =============================================================================


def load_state() -> dict:
    """Load experiment state from disk."""
    state_path = DATA_DIR / "state.json"
    if state_path.exists():
        return json.loads(state_path.read_text())
    return {
        "current_iteration": 0,
        "current_composition": deepcopy(INITIAL_COMPOSITION),
        "alpha": ALPHA,
        "best_od": None,
        "prev_center_od": None,
        "no_improvement_count": 0,
        "converged": False,
        "history": [],
        "novel_bio_used_ul": 0.0,
        "reagent_plate_barcode": None,
        "reagent_plate_volumes": {},  # {well: volume_uL} — current tracked volumes
    }


def save_state(state: dict):
    """Save experiment state to disk."""
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    state_path = DATA_DIR / "state.json"
    state_path.write_text(json.dumps(state, indent=2, default=str))


# =============================================================================
# COMPOSITION HELPERS
# =============================================================================


def compute_novel_bio(supplements: dict) -> int:
    """Compute Novel_Bio volume to fill remaining reagent volume (before cells)."""
    return REAGENT_VOLUME_UL - sum(supplements.values())


def apply_constraints(supplements: dict) -> dict:
    """Apply volume constraints to a supplement composition.

    Rules:
    - Each supplement in {0} union [MIN_SUPPLEMENT_UL, MAX_SUPPLEMENT_UL]
    - All values are integers
    - Novel_Bio must be >= MIN_NOVEL_BIO_UL
    - Sum of all components = WELL_VOLUME_UL
    """
    result = {}
    for name in SUPPLEMENT_NAMES:
        vol = int(round(supplements.get(name, 0)))
        if vol < MIN_SUPPLEMENT_UL:
            vol = 0  # Below minimum -> drop to 0
        vol = min(vol, MAX_SUPPLEMENT_UL)
        result[name] = vol

    # Ensure Novel_Bio >= MIN_NOVEL_BIO_UL
    novel_bio = compute_novel_bio(result)
    while novel_bio < MIN_NOVEL_BIO_UL:
        # Reduce the largest supplement by DELTA_UL
        largest = max(SUPPLEMENT_NAMES, key=lambda n: result[n])
        if result[largest] <= 0:
            break
        result[largest] = max(0, result[largest] - DELTA_UL)
        novel_bio = compute_novel_bio(result)

    return result


def make_perturbed(center: dict, supplement_name: str, delta: int) -> dict:
    """Create a perturbed composition: center + delta on one axis."""
    perturbed = deepcopy(center)
    perturbed[supplement_name] = center[supplement_name] + delta
    return apply_constraints(perturbed)


# =============================================================================
# TRANSFER ARRAY GENERATION
# =============================================================================


def generate_transfer_array(
    center: dict,
    row_letter: str,
    delta: int = DELTA_UL,
    novel_bio_well: str = None,
) -> list:
    """Generate a transfer array for one iteration (9 filled wells in 1 row).

    Returns: [[source_well, dest_well, volume_uL], ...]

    Column layout (all in the same row):
      Col 1:  Seed well (pre-loaded, not part of transfer array)
      Col 2:  Empty buffer (no transfer — isolates seed well from experiments)
      Col 3:  Positive control (180 uL Novel_Bio, cells added by seeding)
      Col 4:  Center point (current best recipe)
      Col 5-6:  +delta Supplement 1 (2 reps)
      Col 7-8:  +delta Supplement 2 (2 reps)
      Col 9-10: +delta Supplement 3 (2 reps)
      Col 11: Empty buffer (no transfer — isolates experiments from neg control)
      Col 12: Negative control (200 uL Novel_Bio, no cells)
    """
    if novel_bio_well is None:
        novel_bio_well = REAGENT_WELLS["Novel_Bio"]

    transfers = []
    row = row_letter

    # Col 3: Positive control — 180 uL Novel_Bio (cells added by seeding step)
    transfers.append([novel_bio_well, f"{row}3", REAGENT_VOLUME_UL])

    # Col 4: Center point
    novel_bio_center = compute_novel_bio(center)
    if novel_bio_center > 0:
        transfers.append([novel_bio_well, f"{row}4", novel_bio_center])
    for name in SUPPLEMENT_NAMES:
        if center[name] > 0:
            transfers.append([REAGENT_WELLS[name], f"{row}4", center[name]])

    # Cols 5-10: Perturbations (2 reps each for 3 supplements)
    for col1, col2, supplement in PERTURBATION_COLS:
        perturbed = make_perturbed(center, supplement, delta)
        novel_bio_pert = compute_novel_bio(perturbed)

        for col in (col1, col2):
            if novel_bio_pert > 0:
                transfers.append([novel_bio_well, f"{row}{col}", novel_bio_pert])
            for name in SUPPLEMENT_NAMES:
                if perturbed[name] > 0:
                    transfers.append([REAGENT_WELLS[name], f"{row}{col}", perturbed[name]])

    # Col 12: Negative control — 200 uL Novel_Bio (no cells seeded)
    transfers.append([novel_bio_well, f"{row}12", WELL_VOLUME_UL])

    # Sort: Novel_Bio well first (reuse tip), then supplements in reverse order
    source_order = [
        novel_bio_well,
        REAGENT_WELLS["DiH2O"],
        REAGENT_WELLS["MOPS_1M"],
        REAGENT_WELLS["Glucose_100mg_mL"],
    ]
    order_map = {well: i for i, well in enumerate(source_order)}
    transfers.sort(key=lambda t: order_map.get(t[0], 99))

    return transfers


def novel_bio_volume_for_array(transfer_array: list, novel_bio_well: str) -> float:
    """Sum total Novel Bio volume drawn from novel_bio_well in a transfer array."""
    return sum(float(vol) for src, _, vol in transfer_array if src == novel_bio_well)


def compute_tip_consumption(transfer_array: list, novel_bio_well: str = None) -> dict:
    """Compute tip consumption for hybrid tip mode.

    Novel_Bio well: 1 tip reused for all transfers (grouped by pipette type).
    All other source wells: 1 new tip per transfer.
    """
    nb_well = novel_bio_well or REAGENT_WELLS["Novel_Bio"]
    reuse_set = {nb_well}

    # Reuse wells: 1 tip per unique (source_well, pipette_type)
    reuse_by_pipette = {"p50": set(), "p200": set(), "p1000": set()}
    # Single-use wells: count each transfer
    single_counts = {"p50": 0, "p200": 0, "p1000": 0}

    for source_well, _, volume in transfer_array:
        if volume <= 50:
            pip = "p50"
        elif volume <= 200:
            pip = "p200"
        else:
            pip = "p1000"

        if source_well in reuse_set:
            reuse_by_pipette[pip].add(source_well)
        else:
            single_counts[pip] += 1

    return {
        "p50": len(reuse_by_pipette["p50"]) + single_counts["p50"],
        "p200": len(reuse_by_pipette["p200"]) + single_counts["p200"],
        "p1000": len(reuse_by_pipette["p1000"]) + single_counts["p1000"],
    }


# =============================================================================
# GRADIENT COMPUTATION
# =============================================================================


def compute_gradient(center_od: float, perturbed_ods: dict, delta: int = DELTA_UL) -> dict:
    """Compute gradient via forward finite differences.

    :param center_od: OD600 of center point (row B)
    :param perturbed_ods: {supplement_name: [rep1_od, rep2_od]}
    :param delta: perturbation step size
    :returns: {supplement_name: gradient_value}
    """
    gradient = {}
    for name in SUPPLEMENT_NAMES:
        mean_perturbed = sum(perturbed_ods[name]) / len(perturbed_ods[name])
        gradient[name] = (mean_perturbed - center_od) / delta

    return gradient


def gradient_step(
    center: dict,
    gradient: dict,
    alpha: float,
    delta: int = DELTA_UL,
) -> dict:
    """Apply sign-based gradient ascent step.

    x += alpha * delta * sign(gradient)

    Sign-based stepping is natural at the 10 uL resolution — each
    supplement moves by exactly one delta in the direction of improvement.
    """
    new_composition = deepcopy(center)
    for name in SUPPLEMENT_NAMES:
        if gradient[name] > 0:
            step = int(alpha * delta)
        elif gradient[name] < 0:
            step = -int(alpha * delta)
        else:
            step = 0
        new_composition[name] = center[name] + step

    return apply_constraints(new_composition)


# =============================================================================
# SEED PARAMETERS (combined routine)
# =============================================================================


def get_seed_params(iteration: int) -> dict:
    """Get seed well parameters for an iteration of the combined routine.

    Row-wise layout: each iteration uses one row (A-H).
    Seed well is column 1 of the iteration's row.
    Cells are seeded into cols 3-10 only:
      - Col 2 skipped: empty buffer (isolates seed well from experiments)
      - Cols 11-12 skipped: col 11 is empty buffer, col 12 is neg control (no cells)

    Iteration 1: seed from A1, experiments in row A, warm up B1
    Iteration 8: seed from H1, experiments in row H, no warmup
    """
    row = ROWS[iteration - 1]
    seed_well = f"{row}1"
    is_last = iteration >= MAX_ITERATIONS
    next_seed_well = f"{ROWS[iteration]}1" if not is_last else "B1"  # dummy for last
    nm_cells_volume = 0 if is_last else 220

    # Wells to seed with cells: cols 3-10 (skip col 2 = buffer, skip cols 11-12 = buffer/neg ctrl)
    seed_dest_wells = [f"{row}{col}" for col in range(3, 11)]

    return {
        "seed_well": seed_well,
        "next_seed_well": next_seed_well,
        "nm_cells_volume": nm_cells_volume,
        "seed_dest_wells": seed_dest_wells,
    }


# =============================================================================
# REAGENT PLATE SETUP
# =============================================================================


def _reagent_plate_fill_plan(state: dict) -> dict:
    """Compute fill volumes (uL) for each well of the AGD reagent plate.

    Returns {well: volume_uL} covering all wells that need to be pre-loaded.
    Novel Bio: fills whichever wells are needed based on cumulative usage so far.
    Supplements: fixed SUPPLEMENT_FILL_UL per well (conservative upper bound).
    NM+Cells: fixed NM_CELLS_FILL_UL.
    """
    nb_used = state.get("novel_bio_used_ul", 0.0)
    remaining_iters = MAX_ITERATIONS - state["current_iteration"]

    # Estimate Novel Bio volume needed for remaining iterations
    sample_ta = generate_transfer_array(
        state["current_composition"], ROWS[0], novel_bio_well=NOVEL_BIO_WELLS[0]
    )
    nb_per_iter = novel_bio_volume_for_array(sample_ta, NOVEL_BIO_WELLS[0])
    nb_needed = nb_per_iter * remaining_iters

    # Determine which Novel Bio wells to fill
    fill = {}
    first_well_idx = int(nb_used // NOVEL_BIO_WELL_CAPACITY_UL)
    remaining_in_current = NOVEL_BIO_WELL_CAPACITY_UL - (nb_used % NOVEL_BIO_WELL_CAPACITY_UL)
    accounted = 0.0
    for i in range(first_well_idx, len(NOVEL_BIO_WELLS)):
        w = NOVEL_BIO_WELLS[i]
        if i == first_well_idx:
            # Partial fill: only remaining capacity in the current well
            fill[w] = int(remaining_in_current)
        else:
            fill[w] = NOVEL_BIO_WELL_CAPACITY_UL
        accounted += fill[w]
        if accounted >= nb_needed:
            break

    # Supplements
    for name in SUPPLEMENT_NAMES:
        fill[REAGENT_WELLS[name]] = SUPPLEMENT_FILL_UL

    # NM+Cells
    fill[NM_CELLS_REAGENT_WELL] = NM_CELLS_FILL_UL

    return fill


def print_reagent_plate_setup(state: dict):
    """Print instructions for filling the AGD reagent plate for the next experiment."""
    fill = _reagent_plate_fill_plan(state)
    nb_used = state.get("novel_bio_used_ul", 0.0)
    remaining_iters = MAX_ITERATIONS - state["current_iteration"]

    print(f"\n{'='*60}")
    print(f"  AGD REAGENT PLATE SETUP")
    print(f"  (Plate type: \"{AGD_PLATE_TYPE}\")")
    print(f"{'='*60}")
    print(f"  {remaining_iters} iterations remaining. Fill a 24-well deep well plate:\n")

    # Supplements
    for name in SUPPLEMENT_NAMES:
        well = REAGENT_WELLS[name]
        print(f"  {well:3s}  {name:<22s}  →  {fill[well]:>5} uL")

    # NM+Cells
    print(f"  {NM_CELLS_REAGENT_WELL:3s}  {'NM+Cells':<22s}  →  {fill[NM_CELLS_REAGENT_WELL]:>5} uL")

    # Novel Bio wells
    for w in NOVEL_BIO_WELLS:
        if w in fill:
            already_used = nb_used % NOVEL_BIO_WELL_CAPACITY_UL if w == get_novel_bio_well(nb_used) and nb_used > 0 else 0
            note = f"  (~{already_used:.0f} uL already consumed)" if already_used > 0 else ""
            print(f"  {w:3s}  {'Novel_Bio':<22s}  →  {fill[w]:>5} uL{note}")

    print(f"\n  Load the plate, check it into the workcell, then run the experiment.")
    print(f"{'='*60}\n")


def find_and_setup_reagent_plate(state: dict) -> str | None:
    """Find a freshly loaded AGD Stock Plate and register its reagents via MCP.

    Polls list_reagent_plates until a plate with:
      - initial_media_type == AGD_PLATE_TYPE
      - reagents_by_well == {} (not yet configured)
    ...appears. Then calls set_reagents_by_well to populate it.

    Returns the plate barcode on success, None on timeout.
    """
    fill = _reagent_plate_fill_plan(state)

    # Build reagents_by_well payload for set_reagents_by_well
    reagents_payload = {}

    # Supplements
    for name in SUPPLEMENT_NAMES:
        well = REAGENT_WELLS[name]
        vol = fill.get(well, 0)
        if vol > 0:
            reagents_payload[well] = [{"name": name, "volume": f"{vol} uL"}]

    # NM+Cells
    nm_vol = fill.get(NM_CELLS_REAGENT_WELL, 0)
    if nm_vol > 0:
        reagents_payload[NM_CELLS_REAGENT_WELL] = [
            {"name": NM_CELLS_REAGENT_NAME, "volume": f"{nm_vol} uL"}
        ]

    # Novel Bio wells
    for w in NOVEL_BIO_WELLS:
        vol = fill.get(w, 0)
        if vol > 0:
            reagents_payload[w] = [{"name": "Novel_Bio", "volume": f"{vol} uL"}]

    print(f"\n  Waiting for an empty {AGD_PLATE_TYPE} to be checked in...")
    print(f"  (Load + check in the plate as shown above, then the system will configure it)")

    start = time.time()
    while time.time() - start < REAGENT_PLATE_POLL_TIMEOUT:
        try:
            plates = _mcp_client.call_tool(
                "list_reagent_plates",
                {"is_checked_in": True, "include_reagents": True},
            )
        except Exception as e:
            print(f"    WARNING: Could not list reagent plates: {e}. Retrying...")
            time.sleep(REAGENT_PLATE_POLL_INTERVAL)
            continue

        for plate in (plates if isinstance(plates, list) else []):
            if plate.get("initial_media_type") != AGD_PLATE_TYPE:
                continue
            if plate.get("reagents_by_well"):  # already configured — skip
                continue

            barcode = plate["barcode"]
            print(f"  Found empty {AGD_PLATE_TYPE}: {barcode}")

            try:
                _mcp_client.call_tool(
                    "set_reagents_by_well",
                    {"plate_barcode": barcode, "reagents_by_well": reagents_payload},
                )
            except Exception as e:
                print(f"  ERROR: Could not set reagents on {barcode}: {e}")
                return None

            print(f"  Reagents registered on {barcode}:")
            for well in sorted(reagents_payload):
                r = reagents_payload[well][0]
                print(f"    {well}: {r['name']} ({r['volume']})")
            return barcode

        elapsed = int(time.time() - start)
        print(f"    [{elapsed}s] No empty {AGD_PLATE_TYPE} found. Retrying in {REAGENT_PLATE_POLL_INTERVAL}s...")
        time.sleep(REAGENT_PLATE_POLL_INTERVAL)

    print(f"  TIMEOUT: No {AGD_PLATE_TYPE} found within {REAGENT_PLATE_POLL_TIMEOUT // 60} minutes.")
    return None


def update_reagent_plate_volumes(state: dict, transfer_array: list, novel_bio_well: str, nm_cells_volume: float):
    """Deduct this iteration's consumed volumes from the tracked reagent plate volumes.

    Reads state["reagent_plate_volumes"], subtracts what was consumed, updates
    state in place, and calls set_reagents_by_well so the MCP stays in sync.
    """
    barcode = state.get("reagent_plate_barcode")
    if not barcode:
        return

    tracked = dict(state.get("reagent_plate_volumes", {}))

    # Compute volumes consumed per source well this iteration
    consumed: dict[str, float] = {}
    for src, _, vol in transfer_array:
        consumed[src] = consumed.get(src, 0.0) + float(vol)
    if nm_cells_volume > 0:
        consumed[NM_CELLS_REAGENT_WELL] = consumed.get(NM_CELLS_REAGENT_WELL, 0.0) + nm_cells_volume

    # Deduct and warn on unexpected negatives
    for well, used in consumed.items():
        prev = tracked.get(well, 0.0)
        tracked[well] = max(0.0, prev - used)
        status = f"{prev:.0f} → {tracked[well]:.0f} uL"
        if prev > 0 and prev - used < 0:
            print(f"  WARNING: {well} over-consumed ({prev:.0f} uL available, {used:.0f} uL drawn)")
        else:
            print(f"  Reagent {well}: consumed {used:.0f} uL  ({status})")

    # Build full reagents_by_well payload for set_reagents_by_well
    # Map well -> reagent name
    well_to_name = {v: k for k, v in REAGENT_WELLS.items()}
    well_to_name[NM_CELLS_REAGENT_WELL] = NM_CELLS_REAGENT_NAME
    for w in NOVEL_BIO_WELLS:
        well_to_name[w] = "Novel_Bio"

    reagents_payload = {}
    for well, vol in tracked.items():
        name = well_to_name.get(well, well)
        reagents_payload[well] = [{"name": name, "volume": f"{int(vol)} uL"}]

    try:
        _mcp_client.call_tool(
            "set_reagents_by_well",
            {"plate_barcode": barcode, "reagents_by_well": reagents_payload},
        )
        print(f"  Reagent plate {barcode} volumes updated in MCP.")
    except Exception as e:
        print(f"  WARNING: Could not update reagent volumes in MCP: {e}")

    state["reagent_plate_volumes"] = tracked


# =============================================================================
# WORKCELL COMMUNICATION (MCP HTTP + REST)
# =============================================================================


def check_workcell_connectivity():
    """Verify we can reach the workcell MCP server."""
    try:
        _mcp_client.connect()
        result = _mcp_client.call_tool("list_workflow_definitions", {})
        count = len(result) if isinstance(result, list) else "?"
        print(f"  Workcell MCP reachable at {WORKCELL_API_BASE} ({count} definitions)")
        return True
    except Exception as e:
        print(f"  ERROR: Cannot reach workcell MCP at {WORKCELL_API_BASE}: {e}")
        return False


class McpClient:
    """Lightweight MCP client that calls tools via HTTP Streamable Transport.

    The workcell's FastMCP server is mounted at /mcp and accepts JSON-RPC 2.0
    tool calls over HTTP POST. Each session requires an initialize handshake.
    """

    def __init__(self, base_url: str = WORKCELL_API_BASE):
        self.mcp_url = f"{base_url}/mcp"
        self.session_id = None
        self._next_id = 1

    def _get_id(self) -> int:
        id_ = self._next_id
        self._next_id += 1
        return id_

    def connect(self):
        """Initialize the MCP session."""
        # Step 1: Initialize
        resp = requests.post(
            self.mcp_url,
            headers={
                "Content-Type": "application/json",
                "Accept": "application/json, text/event-stream",
            },
            json={
                "jsonrpc": "2.0",
                "id": self._get_id(),
                "method": "initialize",
                "params": {
                    "protocolVersion": "2024-11-05",
                    "capabilities": {},
                    "clientInfo": {"name": "gradient_descent_daemon", "version": "1.0"},
                },
            },
            timeout=15,
        )
        resp.raise_for_status()
        self.session_id = resp.headers.get("mcp-session-id")
        if not self.session_id:
            raise RuntimeError("MCP server did not return a session ID")

        # Step 2: Send initialized notification
        requests.post(
            self.mcp_url,
            headers={
                "Content-Type": "application/json",
                "Mcp-Session-Id": self.session_id,
            },
            json={"jsonrpc": "2.0", "method": "notifications/initialized"},
            timeout=10,
        )

    def call_tool(self, tool_name: str, arguments: dict, timeout: int = 30):
        """Call an MCP tool and return the parsed result.

        Handles both structuredContent and text content responses.
        """
        if not self.session_id:
            self.connect()

        resp = requests.post(
            self.mcp_url,
            headers={
                "Content-Type": "application/json",
                "Accept": "application/json, text/event-stream",
                "Mcp-Session-Id": self.session_id,
            },
            json={
                "jsonrpc": "2.0",
                "id": self._get_id(),
                "method": "tools/call",
                "params": {"name": tool_name, "arguments": arguments},
            },
            timeout=timeout,
        )
        resp.raise_for_status()

        # Parse SSE response (event: message\ndata: {...})
        body = resp.text
        for line in body.split("\n"):
            if line.startswith("data: "):
                payload = json.loads(line[6:])
                result = payload.get("result", {})
                if result.get("isError"):
                    error_text = result.get("content", [{}])[0].get("text", "Unknown error")
                    raise RuntimeError(f"MCP tool error: {error_text}")

                # Prefer structuredContent, fall back to content[0].text
                sc = result.get("structuredContent", {}).get("result")
                if sc is not None:
                    return sc
                content = result.get("content", [])
                if content and content[0].get("text"):
                    try:
                        return json.loads(content[0]["text"])
                    except (json.JSONDecodeError, KeyError):
                        return content[0]["text"]
                return result

        raise RuntimeError(f"Could not parse MCP response: {body[:500]}")


# Shared MCP client instance (lazy-connected)
_mcp_client = McpClient()


# =============================================================================
# WORKFLOW MANAGEMENT
# =============================================================================


def write_workflow_definition(
    transfer_array: list,
    row_letter: str,
    iteration: int,
    seed_params: dict | None = None,
    novel_bio_well: str = None,
) -> Path:
    """Write an iteration-specific workflow definition file.

    Copies the template and replaces iteration-specific constants using regex
    to handle any default values in the template.
    """
    template = WORKFLOW_TEMPLATE_PATH.read_text()

    nb_well = novel_bio_well or REAGENT_WELLS["Novel_Bio"]
    tip_counts = compute_tip_consumption(transfer_array, novel_bio_well=nb_well)
    # Reagent well counter: count supplement source wells only (exclude Novel Bio,
    # which is tracked separately via novel_bio_used_ul in state).
    supplement_wells_used = len(set(t[0] for t in transfer_array) - {nb_well})

    def replace_const(name: str, value):
        """Replace a module-level constant in the template by name."""
        nonlocal template
        template = re.sub(
            rf"^({re.escape(name)}\s*=\s*).*$",
            rf"\g<1>{value}",
            template,
            flags=re.MULTILINE,
        )

    # Core parameters
    replace_const("TRANSFER_ARRAY", json.dumps(json.dumps(transfer_array)))
    replace_const("DEST_ROW", f'"{row_letter}"')

    # Seed parameters (for combined routine template)
    if seed_params:
        replace_const("SEED_WELL", f'"{seed_params["seed_well"]}"')
        replace_const("NEXT_SEED_WELL", f'"{seed_params["next_seed_well"]}"')
        replace_const("NM_CELLS_VOLUME", str(seed_params["nm_cells_volume"]))
        replace_const("SEED_DEST_WELLS", json.dumps(seed_params["seed_dest_wells"]))

    # Tip consumption — add extras for combined routine (seed mixing, seeding, NM warmup)
    # Seeding uses a unique P50 tip per destination well (cols 3-10 = 8 wells)
    n_seed_wells = len(seed_params["seed_dest_wells"]) if seed_params else 0
    p50_extra = n_seed_wells   # unique tip per seed dest well
    p200_extra = 1 if seed_params else 0  # +1 P200 for seed well mixing
    p1000_count = 1 if seed_params and seed_params["nm_cells_volume"] > 0 else 0
    reagent_extra = 1 if seed_params and seed_params["nm_cells_volume"] > 0 else 0

    replace_const("P50_TIPS_TO_CONSUME", str(tip_counts["p50"] + p50_extra))
    replace_const("P200_TIPS_TO_CONSUME", str(tip_counts["p200"] + p200_extra))
    replace_const("P1000_TIPS_TO_CONSUME", str(tip_counts.get("p1000", 0) + p1000_count))
    replace_const("REAGENT_WELLS_TO_CONSUME", str(supplement_wells_used + reagent_extra))

    # Write iteration-specific file
    iter_dir = DATA_DIR / f"iteration_{iteration}"
    iter_dir.mkdir(parents=True, exist_ok=True)
    output_path = iter_dir / "workflow_definition.py"
    output_path.write_text(template)

    return output_path


def register_workflow(workflow_path: Path, iteration: int) -> int:
    """Register a workflow definition via MCP (upload file + create DB record).

    Returns the workflow definition database ID.
    """
    file_name = f"gradient_descent_iteration_r{iteration}.py"
    workflow_name = f"Gradient Descent Iteration {iteration}"

    # Step 1: Upload the workflow definition file via MCP
    code_content = workflow_path.read_text()
    _mcp_client.call_tool("create_workflow_definition_file", {
        "file_name": file_name,
        "code_content": code_content,
    })
    print(f"    MCP: uploaded {file_name}")

    # Step 2: Register the workflow definition (creates DB record)
    _mcp_client.call_tool("register_workflow_definition", {
        "name": workflow_name,
        "file_name": file_name,
    })
    print(f"    MCP: registered '{workflow_name}'")

    # Step 3: Look up the definition ID
    definitions = _mcp_client.call_tool("list_workflow_definitions", {})
    for d in definitions:
        if d["name"] == workflow_name:
            return d["id"]

    raise RuntimeError(f"Definition '{workflow_name}' not found after registration")


def instantiate_workflow(definition_id: int, plate_barcode: str, iteration: int) -> str:
    """Create a workflow instance via MCP HTTP. Returns the instance UUID.

    Calls the instantiate_workflow MCP tool over HTTP Streamable Transport.
    With auto_approve_pending_instances=True, the workflow starts immediately.
    """
    reason = f"Gradient descent iteration {iteration} for plate {plate_barcode}"

    result = _mcp_client.call_tool("instantiate_workflow", {
        "definition_id": definition_id,
        "inputs": {"plate_barcode": plate_barcode},
        "reason": reason,
    })

    uuid = result.get("uuid")
    if not uuid:
        raise RuntimeError(f"MCP instantiate_workflow did not return UUID: {result}")

    status = result.get("status", "unknown")
    print(f"    MCP instantiate: uuid={uuid}, status={status}")
    return uuid


def poll_workflow_completion(
    uuid: str, timeout_minutes: int = WORKFLOW_TIMEOUT_MINUTES
) -> dict:
    """Poll via MCP until a workflow instance completes. Returns the instance data."""
    start_time = time.time()
    deadline = start_time + timeout_minutes * 60

    while time.time() < deadline:
        instance = _mcp_client.call_tool(
            "get_workflow_instance_details", {"instance_uuid": uuid}
        )
        status = instance.get("status", "unknown")

        if status in ("completed", "failed", "cancelled", "canceled"):
            print()  # Clear the \r line
            return instance

        elapsed = int(time.time() - start_time)
        print(
            f"    Status: {status} ({elapsed // 60}m{elapsed % 60:02d}s elapsed)",
            end="\r",
        )
        time.sleep(DAEMON_POLL_INTERVAL)

    raise TimeoutError(
        f"Workflow {uuid} did not complete within {timeout_minutes} minutes"
    )


def _get_plate_uuid(plate_barcode: str) -> str:
    """Look up a culture plate's UUID from its barcode via MCP."""
    plates = _mcp_client.call_tool("list_culture_plates", {})
    for p in plates:
        if p.get("barcode") == plate_barcode:
            return p.get("uuid", "")
    raise RuntimeError(f"Plate '{plate_barcode}' not found on workcell")


def _mcp_call(tool_name: str, arguments: dict) -> dict:
    """Call a tool on the backend MCP server. Returns the structured result."""
    # Initialize session
    init_resp = requests.post(
        BACKEND_MCP_URL,
        headers=BACKEND_MCP_HEADERS,
        json={
            "jsonrpc": "2.0", "id": 1, "method": "initialize",
            "params": {
                "protocolVersion": "2024-11-05",
                "capabilities": {},
                "clientInfo": {"name": "gd-daemon", "version": "1.0"},
            },
        },
        timeout=15,
    )
    init_resp.raise_for_status()
    session_id = init_resp.headers.get("mcp-session-id")
    if not session_id:
        raise RuntimeError("Backend MCP did not return a session ID")

    session_headers = {**BACKEND_MCP_HEADERS, "mcp-session-id": session_id}

    # Send initialized notification
    requests.post(
        BACKEND_MCP_URL,
        headers=session_headers,
        json={"jsonrpc": "2.0", "method": "notifications/initialized", "params": {}},
        timeout=5,
    )

    # Call the tool
    tool_resp = requests.post(
        BACKEND_MCP_URL,
        headers=session_headers,
        json={
            "jsonrpc": "2.0", "id": 2, "method": "tools/call",
            "params": {"name": tool_name, "arguments": arguments},
        },
        timeout=30,
    )
    tool_resp.raise_for_status()

    # Parse SSE response
    for line in tool_resp.text.splitlines():
        if line.startswith("data:"):
            data = json.loads(line[5:].strip())
            if "error" in data:
                raise RuntimeError(f"MCP tool error: {data['error']}")
            structured = data.get("result", {}).get("structuredContent", {}).get("result")
            if structured is not None:
                return structured
            # Fall back to parsing text content
            content = data.get("result", {}).get("content", [])
            for item in content:
                if item.get("type") == "text":
                    return json.loads(item["text"])
    raise RuntimeError(f"No result from MCP tool {tool_name}")


def fetch_absorbance_results(plate_barcode: str, row_letter: str) -> dict:
    """Fetch OD600 readings for a plate via the backend MCP.

    Uses get_plate_observations to get all datasets for the plate, filters
    to datasets containing target row wells, then returns earliest (baseline)
    and latest (endpoint) readings.

    Returns: {
        "baseline": {well: od600_value},   # First reading post-seeding
        "endpoint": {well: od600_value},   # Latest reading (post-growth)
    }
    """
    result = _mcp_call("get_plate_observations", {"plate_name": plate_barcode, "limit": 500})

    datasets = result.get("datasets", [])
    if not datasets:
        raise RuntimeError(f"No datasets found for plate {plate_barcode}")

    # Target wells: cols 3-10 (experiments) + col 12 (neg control). Cols 2 and 11 are empty.
    target_wells = [f"{row_letter}{col}" for col in list(range(3, 11)) + [12]]
    row_readings = {}
    for ds in datasets:
        obs = ds.get("observations_by_well", {})
        if any(w in obs for w in target_wells):
            ts = ds["timestamp"]
            row_readings[ts] = {
                well: float(data["absorbance"])
                for well, data in obs.items()
                if data.get("absorbance") is not None
            }

    if not row_readings:
        raise RuntimeError(
            f"No OD600 readings found for row {row_letter} wells on plate {plate_barcode}"
        )

    sorted_timestamps = sorted(row_readings.keys())
    earliest_well_data = row_readings[sorted_timestamps[0]]
    latest_well_data = row_readings[sorted_timestamps[-1]]

    print(f"    Row {row_letter} readings: {sorted_timestamps[0]} -> {sorted_timestamps[-1]} "
          f"({len(sorted_timestamps)} readings)")

    baseline = {}
    endpoint = {}
    for col in list(range(3, 11)) + [12]:
        well = f"{row_letter}{col}"
        baseline[well] = earliest_well_data.get(well, 0.0)
        endpoint[well] = latest_well_data.get(well, 0.0)

    return {"baseline": baseline, "endpoint": endpoint}


# =============================================================================
# RESULT PARSING
# =============================================================================


def parse_od_results(absorbance_results: dict, row_letter: str) -> dict:
    """Parse baseline + endpoint absorbance into growth deltas for gradient computation.

    Uses delta OD (endpoint - baseline) as the objective function, which normalizes
    for varying initial cell densities across wells.

    Row-wise layout:
      Col 2:  Empty buffer (ignored)
      Col 3:  Positive control (cells in plain media)
      Col 4:  Center point
      Cols 5-10: Perturbation wells (3 supplements x 2 reps)
      Col 11: Empty buffer (ignored)
      Col 12: Negative control (no cells)

    Returns:
        {
            "neg_control_od": float,   # delta OD for negative control (no cells, col 12)
            "control_od": float,       # delta OD for positive control (col 3)
            "center_od": float,        # delta OD for center well (col 4)
            "perturbed_ods": {supplement: [rep1_delta, rep2_delta]},
            "abs_control_od": float,   # absolute endpoint (for logging)
            "abs_center_od": float,    # absolute endpoint (for logging)
        }
    """
    row = row_letter
    baseline = absorbance_results.get("baseline", {})
    endpoint = absorbance_results.get("endpoint", {})

    def delta(well):
        return endpoint.get(well, 0.0) - baseline.get(well, 0.0)

    neg_control_delta = delta(f"{row}12")
    pos_control_delta = delta(f"{row}3")
    center_delta = delta(f"{row}4")

    perturbed_deltas = {}
    for col1, col2, supplement in PERTURBATION_COLS:
        perturbed_deltas[supplement] = [delta(f"{row}{col1}"), delta(f"{row}{col2}")]

    return {
        "neg_control_od": neg_control_delta,
        "control_od": pos_control_delta,
        "center_od": center_delta,
        "perturbed_ods": perturbed_deltas,
        "abs_control_od": endpoint.get(f"{row}3", 0.0),
        "abs_center_od": endpoint.get(f"{row}4", 0.0),
    }


# =============================================================================
# LOGGING
# =============================================================================


def log_iteration(
    iteration: int,
    center: dict,
    transfer_array: list,
    od_results: dict = None,
    gradient: dict = None,
    new_composition: dict = None,
    alpha: float = None,
):
    """Log iteration details to disk and stdout."""
    iter_dir = DATA_DIR / f"iteration_{iteration}"
    iter_dir.mkdir(parents=True, exist_ok=True)

    novel_bio = compute_novel_bio(center)
    log = {
        "iteration": iteration,
        "timestamp": datetime.now().isoformat(),
        "composition": {
            "Novel_Bio": novel_bio,
            **center,
        },
        "transfer_count": len(transfer_array),
        "alpha": alpha,
    }

    if od_results:
        log["od_results"] = od_results
    if gradient:
        log["gradient"] = gradient
    if new_composition:
        log["next_composition"] = {
            "Novel_Bio": compute_novel_bio(new_composition),
            **new_composition,
        }

    (iter_dir / "iteration_log.json").write_text(json.dumps(log, indent=2))

    # Also save transfer array
    (iter_dir / "transfer_array.json").write_text(json.dumps(transfer_array, indent=2))

    # Print summary
    print(f"\n{'='*60}")
    print(f"  ITERATION {iteration} / {MAX_ITERATIONS}")
    print(f"{'='*60}")
    print(f"  Composition: Novel_Bio={novel_bio}, "
          + ", ".join(f"{k}={v}" for k, v in center.items()))
    print(f"  Transfers: {len(transfer_array)}")

    if od_results:
        print(f"\n  OD600 Growth (delta OD = endpoint - baseline):")
        if "neg_control_od" in od_results:
            print(f"    Neg Ctrl: {od_results['neg_control_od']:+.4f}  (no cells)")
        print(f"    Pos Ctrl: {od_results['control_od']:+.4f}"
              + (f"  (abs: {od_results['abs_control_od']:.4f})" if "abs_control_od" in od_results else ""))
        print(f"    Center:   {od_results['center_od']:+.4f}"
              + (f"  (abs: {od_results['abs_center_od']:.4f})" if "abs_center_od" in od_results else ""))
        for name in SUPPLEMENT_NAMES:
            reps = od_results["perturbed_ods"][name]
            print(f"    +d {name}: {reps[0]:+.4f}, {reps[1]:+.4f} (mean={sum(reps)/2:+.4f})")
        if "extra_ods" in od_results:
            extras = od_results["extra_ods"]
            print(f"    Extra:  {extras[0]:+.4f}, {extras[1]:+.4f} (mean={sum(extras)/2:+.4f})")

    if gradient:
        print(f"\n  Gradient:")
        for name in SUPPLEMENT_NAMES:
            direction = "+" if gradient[name] > 0 else ("-" if gradient[name] < 0 else "=")
            print(f"    {name}: {gradient[name]:+.6f} ({direction})")

    if new_composition:
        new_nb = compute_novel_bio(new_composition)
        print(f"\n  Next composition: Novel_Bio={new_nb}, "
              + ", ".join(f"{k}={v}" for k, v in new_composition.items()))


# =============================================================================
# MAIN ORCHESTRATION LOOP
# =============================================================================


def run_iteration(
    state: dict,
    plate_barcode: str,
    dry_run: bool = False,
) -> dict:
    """Run a single iteration of gradient descent.

    1. Generate transfer array from current composition
    2. Write workflow definition with iteration-specific constants
    3. Upload + register workflow via MCP
    4. Instantiate workflow via MCP (auto-approved)
    5. Poll for completion via MCP
    6. Fetch and parse OD600 results (REST API)
    7. Compute gradient and update composition
    8. Check convergence
    """
    iteration = state["current_iteration"] + 1
    row_letter = ROWS[iteration - 1]  # Iteration 1 = row A, 2 = row B, ...
    center = state["current_composition"]
    alpha = state["alpha"]

    # Seed parameters for combined routine
    seed_params = get_seed_params(iteration)

    print(f"\nPreparing iteration {iteration} (row {row_letter})...")
    print(f"  Seed well: {seed_params['seed_well']}")
    if seed_params["nm_cells_volume"] > 0:
        print(f"  NM+Cells warmup: -> {seed_params['next_seed_well']}")
    else:
        print(f"  Last round — no NM+Cells warmup")

    # --- Resume logic: check for existing workflow from a previous run ---
    iter_dir = DATA_DIR / f"iteration_{iteration}"
    wf_ids_path = iter_dir / "workflow_ids.json"
    resume_phase = None  # None = fresh start, "poll" = resume polling, "results" = skip to results

    if wf_ids_path.exists():
        wf_ids = json.loads(wf_ids_path.read_text())
        existing_uuid = wf_ids.get("workflow_instance_uuid")
        existing_def_id = wf_ids.get("workflow_definition_id")

        if existing_uuid:
            print(f"\n  Found existing workflow for iteration {iteration}: {existing_uuid}")
            try:
                instance = _mcp_client.call_tool(
                    "get_workflow_instance_details", {"instance_uuid": existing_uuid}
                )
                wf_status = instance.get("status", "unknown")
                print(f"  Existing workflow status: {wf_status}")

                if wf_status in ("completed",):
                    # Already done — skip to unlink + results
                    resume_phase = "results"
                    print(f"  Resuming from results phase (workflow already completed)")
                elif wf_status in ("in_progress", "pending", "approved"):
                    # Still running — skip to polling
                    resume_phase = "poll"
                    print(f"  Resuming polling (workflow still in progress)")
                else:
                    # Failed/cancelled — re-run from scratch
                    print(f"  Workflow {wf_status} — will re-register and re-run")
            except Exception as e:
                print(f"  WARNING: Could not check workflow status: {e}")
                print(f"  Will re-register and re-run from scratch")

    # --- Fresh start: generate transfer array + workflow ---
    if resume_phase is None:
        # Step 1: Determine Novel Bio well based on cumulative volume consumed so far
        novel_bio_used_ul = state.get("novel_bio_used_ul", 0.0)
        novel_bio_well = get_novel_bio_well(novel_bio_used_ul)
        prev_well = get_novel_bio_well(max(0, novel_bio_used_ul - 1)) if novel_bio_used_ul > 0 else novel_bio_well
        if novel_bio_well != prev_well:
            print(f"\n  Novel Bio well switch: {prev_well} -> {novel_bio_well} "
                  f"({novel_bio_used_ul:.0f} uL consumed)")
        elif novel_bio_used_ul > 0:
            remaining = NOVEL_BIO_WELL_CAPACITY_UL - (novel_bio_used_ul % NOVEL_BIO_WELL_CAPACITY_UL)
            print(f"  Novel Bio: using {novel_bio_well} "
                  f"({novel_bio_used_ul:.0f} uL used, ~{remaining:.0f} uL remaining in well)")

        # Step 2: Generate transfer array
        # Perturbation delta scales with alpha so we test the actual step size
        perturbation_delta = max(1, int(alpha * DELTA_UL))
        transfer_array = generate_transfer_array(
            center, row_letter, delta=perturbation_delta, novel_bio_well=novel_bio_well
        )
        nb_this_iter = novel_bio_volume_for_array(transfer_array, novel_bio_well)

        # Step 3: Write workflow definition
        workflow_path = write_workflow_definition(
            transfer_array, row_letter, iteration, seed_params,
            novel_bio_well=novel_bio_well,
        )
        print(f"  Workflow definition: {workflow_path}")

        # Log the setup
        log_iteration(iteration, center, transfer_array, alpha=alpha)

        if dry_run:
            print(f"\n  [DRY RUN] Would register and run workflow for row {row_letter}")
            print(f"  Transfer array: {DATA_DIR / f'iteration_{iteration}' / 'transfer_array.json'}")
            return state

        # Step 3: Upload + register workflow via MCP
        print(f"\n  Registering workflow...")
        workflow_def_id = register_workflow(workflow_path, iteration)
        print(f"  Definition ID: {workflow_def_id}")

        # Step 4: Instantiate workflow via MCP (auto-approved if config flag is set)
        print(f"  Instantiating workflow via MCP...")
        workflow_uuid = instantiate_workflow(workflow_def_id, plate_barcode, iteration)
        print(f"  Instance UUID: {workflow_uuid}")

        # Save workflow IDs
        iter_dir.mkdir(parents=True, exist_ok=True)
        (iter_dir / "workflow_ids.json").write_text(json.dumps({
            "workflow_definition_id": workflow_def_id,
            "workflow_instance_uuid": workflow_uuid,
        }, indent=2))

        existing_uuid = workflow_uuid
    else:
        # Resuming — load transfer array from disk for logging
        ta_path = iter_dir / "transfer_array.json"
        if ta_path.exists():
            transfer_array = json.loads(ta_path.read_text())
        else:
            transfer_array = []
        # Recompute Novel Bio volume for state tracking
        novel_bio_used_ul = state.get("novel_bio_used_ul", 0.0)
        novel_bio_well = get_novel_bio_well(novel_bio_used_ul)
        nb_this_iter = novel_bio_volume_for_array(transfer_array, novel_bio_well)

    # --- Poll for completion (fresh or resumed) ---
    if resume_phase in (None, "poll"):
        print(f"\n  Waiting for workflow completion (timeout: {WORKFLOW_TIMEOUT_MINUTES} min)...")
        result = poll_workflow_completion(existing_uuid)
        status = result.get("status", "unknown")
        if status != "completed":
            print(f"  ERROR: Workflow ended with status: {status}")
            print(f"  Check workcell dashboard for details.")
            return state

        print(f"  Workflow completed!")

    # Step 6b: Unlink plate from completed workflow (required before next iteration)
    print(f"  Unlinking plate from workflow...")
    unlink_result = _mcp_client.call_tool(
        "unlink_culture_plate_from_workflow",
        {"plate_barcode": plate_barcode, "cancel_workflow": False},
    )
    if isinstance(unlink_result, dict) and unlink_result.get("success"):
        print(f"    Plate unlinked from {unlink_result.get('workflow_name', 'workflow')}")
    else:
        print(f"    WARNING: Unlink returned: {unlink_result}")

    # Step 7: Fetch results
    print(f"  Fetching OD600 results...")
    absorbance_results = fetch_absorbance_results(plate_barcode, row_letter)
    od_results = parse_od_results(absorbance_results, row_letter)

    # Step 8: Compute gradient
    perturbation_delta = max(1, int(alpha * DELTA_UL))
    gradient = compute_gradient(od_results["center_od"], od_results["perturbed_ods"], delta=perturbation_delta)

    # Adaptive alpha: halve if center OD decreased from previous round
    if state["prev_center_od"] is not None:
        if od_results["center_od"] < state["prev_center_od"]:
            alpha = max(0.1, alpha / 2)
            print(f"\n  Center OD decreased ({state['prev_center_od']:.4f} -> "
                  f"{od_results['center_od']:.4f}), halving alpha to {alpha}")

    # Step 9: Update composition
    new_composition = gradient_step(center, gradient, alpha)

    # Log full iteration with results
    log_iteration(iteration, center, transfer_array, od_results, gradient, new_composition, alpha)

    # Update reagent plate volume tracking in MCP
    seed_nm_vol = seed_params["nm_cells_volume"] if seed_params else 0
    print(f"\n  Updating reagent plate volumes...")
    update_reagent_plate_volumes(state, transfer_array, novel_bio_well, seed_nm_vol)

    # Check convergence: no perturbation improved on center
    any_improvement = any(
        sum(od_results["perturbed_ods"][name]) / 2 > od_results["center_od"]
        for name in SUPPLEMENT_NAMES
    )

    if not any_improvement:
        state["no_improvement_count"] += 1
        print(f"\n  No improvement this round ({state['no_improvement_count']}/{CONVERGENCE_ROUNDS})")
    else:
        state["no_improvement_count"] = 0

    # Update state
    state["current_iteration"] = iteration
    state["current_composition"] = new_composition
    state["alpha"] = alpha
    state["prev_center_od"] = od_results["center_od"]
    state["novel_bio_used_ul"] = state.get("novel_bio_used_ul", 0.0) + nb_this_iter

    if state["best_od"] is None or od_results["center_od"] > state["best_od"]:
        state["best_od"] = od_results["center_od"]

    state["history"].append({
        "iteration": iteration,
        "composition": {
            "Novel_Bio": compute_novel_bio(center),
            **center,
        },
        "center_od": od_results["center_od"],
        "control_od": od_results["control_od"],
        "neg_control_od": od_results.get("neg_control_od"),
        "extra_ods": od_results.get("extra_ods"),
        "abs_center_od": od_results.get("abs_center_od"),
        "abs_control_od": od_results.get("abs_control_od"),
        "gradient": gradient,
        "alpha": alpha,
    })

    if iteration >= MAX_ITERATIONS:
        state["converged"] = True
        print(f"\n  MAX ITERATIONS REACHED ({MAX_ITERATIONS}).")

    save_state(state)
    return state


def run_experiment(
    plate_barcode: str,
    dry_run: bool = False,
):
    """Run the full gradient descent optimization loop."""
    state = load_state()

    if state["converged"]:
        print("Experiment already converged. Use 'status' to see results.")
        return

    if not dry_run:
        print("Checking workcell connectivity...")
        if not check_workcell_connectivity():
            print("Aborting. Check network connection and WORKCELL_HOST env var.")
            return

    print(f"\nStarting gradient descent media optimization")
    print(f"  Plate: {plate_barcode}")
    print(f"  Starting from iteration: {state['current_iteration'] + 1}")
    print(f"  Current composition: {state['current_composition']}")
    print(f"  Max iterations: {MAX_ITERATIONS}")

    # Pre-run Novel Bio volume estimate
    remaining_iters = MAX_ITERATIONS - state["current_iteration"]
    nb_used_so_far = state.get("novel_bio_used_ul", 0.0)
    sample_row = ROWS[state["current_iteration"]]
    sample_transfer = generate_transfer_array(state["current_composition"], sample_row, novel_bio_well=NOVEL_BIO_WELLS[0])
    nb_per_iter = novel_bio_volume_for_array(sample_transfer, NOVEL_BIO_WELLS[0])
    nb_total_estimate = nb_per_iter * remaining_iters
    first_well_idx = int(nb_used_so_far // NOVEL_BIO_WELL_CAPACITY_UL)
    last_well_idx = int((nb_used_so_far + nb_total_estimate) // NOVEL_BIO_WELL_CAPACITY_UL)
    if last_well_idx < len(NOVEL_BIO_WELLS):
        wells_needed = NOVEL_BIO_WELLS[first_well_idx:last_well_idx + 1]
        wells_str = ", ".join(wells_needed)
    else:
        wells_needed = NOVEL_BIO_WELLS[first_well_idx:]
        wells_str = ", ".join(wells_needed) + " (WARNING: may not be enough!)"
    print(f"\n  Novel Bio estimate: ~{nb_per_iter:.0f} uL/iter × {remaining_iters} remaining = ~{nb_total_estimate:.0f} uL total")
    print(f"  Already consumed: {nb_used_so_far:.0f} uL  |  Wells needed: {wells_str}")

    # Reagent plate setup: print instructions, then find and configure the plate via MCP
    if not dry_run:
        if not state.get("reagent_plate_barcode"):
            print_reagent_plate_setup(state)
            reagent_barcode = find_and_setup_reagent_plate(state)
            if reagent_barcode:
                state["reagent_plate_barcode"] = reagent_barcode
                # Initialize volume tracking from the fill plan
                fill = _reagent_plate_fill_plan(state)
                state["reagent_plate_volumes"] = {w: float(v) for w, v in fill.items()}
                save_state(state)
            else:
                print("  WARNING: Proceeding without reagent plate tracking.")
        else:
            print(f"\n  Resuming with reagent plate: {state['reagent_plate_barcode']}")

    while not state["converged"]:
        try:
            state = run_iteration(state, plate_barcode, dry_run)
        except Exception as e:
            print(f"\n  ERROR during iteration: {e}")
            print(f"  State saved. Use 'resume' to retry.")
            save_state(state)
            raise

        if dry_run:
            print("\n[DRY RUN] Stopping after 1 iteration preview.")
            break

    if state["converged"] and not dry_run:
        print_final_summary(state)


def print_final_summary(state: dict):
    """Print a summary of the optimization results."""
    print(f"\n{'='*60}")
    print(f"  GRADIENT DESCENT OPTIMIZATION COMPLETE")
    print(f"{'='*60}")
    print(f"  Total iterations: {state['current_iteration']}")
    print(f"  Best center OD600: {state['best_od']:.4f}")

    final = state["current_composition"]
    final_nb = compute_novel_bio(final)
    print(f"\n  Final composition:")
    print(f"    Novel_Bio: {final_nb} uL")
    for name in SUPPLEMENT_NAMES:
        print(f"    {name}: {final[name]} uL")

    print(f"\n  History:")
    # Build dynamic header from supplement names
    supp_headers = " | ".join(f"{name:>10}" for name in SUPPLEMENT_NAMES)
    supp_dashes = " | ".join(f"{'----------':>10}" for _ in SUPPLEMENT_NAMES)
    print(f"  {'Iter':>4} | {'Novel_Bio':>9} | {supp_headers} | {'Center OD':>9} | {'Control':>7}")
    print(f"  {'----':>4} | {'---------':>9} | {supp_dashes} | {'---------':>9} | {'-------':>7}")
    for h in state["history"]:
        c = h["composition"]
        supp_vals = " | ".join(f"{c.get(name, 0):10d}" for name in SUPPLEMENT_NAMES)
        print(f"  {h['iteration']:4d} | {c['Novel_Bio']:9d} | {supp_vals} | {h['center_od']:9.4f} | {h['control_od']:7.4f}")


def print_status():
    """Print current experiment status."""
    state = load_state()

    if state["current_iteration"] == 0:
        print("No iterations completed yet.")
        print(f"Starting composition: {state['current_composition']}")
        return

    print(f"Gradient Descent Status")
    print(f"  Iterations completed: {state['current_iteration']}")
    print(f"  Converged: {state['converged']}")
    print(f"  Current alpha: {state['alpha']}")
    print(f"  Best center OD600: {state['best_od']}")
    print(f"  No-improvement streak: {state['no_improvement_count']}")
    print(f"\n  Next composition: {state['current_composition']}")
    print(f"  Novel_Bio: {compute_novel_bio(state['current_composition'])} uL")

    if state["history"]:
        print(f"\n  History:")
        for h in state["history"]:
            c = h["composition"]
            supps = ", ".join(f"{name}={c.get(name, 0)}" for name in SUPPLEMENT_NAMES)
            print(f"    Iter {h['iteration']}: NB={c['Novel_Bio']}, "
                  f"{supps} -> OD={h['center_od']:.4f}")


# =============================================================================
# CLI
# =============================================================================


def main():
    if len(sys.argv) < 2:
        print("Usage:")
        print("  python gradient_descent.py run <plate_barcode>     # Fresh start")
        print("  python gradient_descent.py resume <plate_barcode>  # Continue from last state")
        print("  python gradient_descent.py dry-run <plate_barcode> # Preview next iteration")
        print("  python gradient_descent.py status                  # Show current state")
        print("  python gradient_descent.py reset                   # Reset state")
        print()
        print("Environment variables:")
        print(f"  WORKCELL_HOST={WORKCELL_HOST}")
        print(f"  WORKCELL_PORT={WORKCELL_PORT}")
        return

    command = sys.argv[1]
    global DATA_DIR

    # For status/reset, optionally accept a plate barcode to target a specific experiment
    if command in ("status", "reset") and len(sys.argv) >= 3:
        DATA_DIR = Path(__file__).parent / "data" / sys.argv[2]

    if command == "status":
        print_status()

    elif command == "reset":
        import shutil

        state_path = DATA_DIR / "state.json"
        if state_path.exists():
            state_path.unlink()

        # Also clean iteration directories so resume logic doesn't find stale workflows
        removed = 0
        for d in sorted(DATA_DIR.glob("iteration_*")):
            if d.is_dir():
                shutil.rmtree(d)
                removed += 1

        print(f"State reset ({removed} iteration dirs removed). Ready for a new experiment.")

    elif command in ("run", "resume", "dry-run"):
        if len(sys.argv) < 3:
            print(f"Usage: python gradient_descent.py {command} <plate_barcode>")
            return

        plate_barcode = sys.argv[2]
        dry_run = command == "dry-run"

        # Use plate barcode as experiment data directory name
        DATA_DIR = Path(__file__).parent / "data" / plate_barcode
        DATA_DIR.mkdir(parents=True, exist_ok=True)

        if command == "run":
            # Fresh start — check for existing state
            state_path = DATA_DIR / "state.json"
            if state_path.exists():
                print("WARNING: Existing state found. Use 'resume' to continue or 'reset' first.")
                return

        run_experiment(plate_barcode, dry_run)

    else:
        print(f"Unknown command: {command}")


if __name__ == "__main__":
    main()
