#!/usr/bin/env python3
"""Gradient descent media optimization for V. natriegens growth.

Supplements Novel_Bio with 3 reagents (Glucose, NaCl, MgSO4) and optimizes
concentrations via gradient descent to maximize OD600 growth after 1.5 hours.

Each iteration uses 1 row (12 wells), with column 1 as the seed well:
  Col 1:    Seed well (NM+Cells stock, pre-loaded)
  Col 2:    Negative control (200 uL Novel_Bio, no cells)
  Col 3:    Positive control (180 uL Novel_Bio + 20 uL cells)
  Col 4:    Center point (current best composition + cells)
  Col 5-6:  +delta Glucose (2 replicates + cells)
  Col 7-8:  +delta NaCl (2 replicates + cells)
  Col 9-10: +delta MgSO4 (2 replicates + cells)
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
WORKCELL_HOST = os.getenv("WORKCELL_HOST", "192.168.68.55")
WORKCELL_PORT = int(os.getenv("WORKCELL_PORT", "8080"))
WORKCELL_API_BASE = f"http://{WORKCELL_HOST}:{WORKCELL_PORT}"

# API request headers (ClientIdentifierMiddleware requires desktop-frontend)
API_HEADERS = {
    "Content-Type": "application/json",
    "X-Monomer-Client": "desktop-frontend",
}

# Reagent plate well map (24-well deep well, modeled as 96-well in the system)
REAGENT_WELLS = {
    "Glucose": "A1",
    "NaCl": "B1",
    "MgSO4": "C1",
    "Novel_Bio": "D1",
}

# Reagent names in order (matches cols 5/6, 7/8, 9/10)
SUPPLEMENT_NAMES = ["Glucose", "NaCl", "MgSO4"]

# Plate layout: column -> purpose (row-wise iteration)
COL_LABELS = {
    1: "seed",
    2: "neg_control",
    3: "pos_control",
    4: "center",
    5: "glucose_rep1",
    6: "glucose_rep2",
    7: "nacl_rep1",
    8: "nacl_rep2",
    9: "mgso4_rep1",
    10: "mgso4_rep2",
    11: "extra1",
    12: "extra2",
}

# Perturbation column pairs: (col1, col2, supplement_name)
PERTURBATION_COLS = [
    (5, 6, "Glucose"),
    (7, 8, "NaCl"),
    (9, 10, "MgSO4"),
]

# Volumes
REAGENT_VOLUME_UL = 180   # Reagent mix volume per well (before cells)
SEED_TRANSFER_VOLUME = 20  # Cells added from seed well
WELL_VOLUME_UL = 200       # Total volume (reagent + cells)

# Constraints
MIN_SUPPLEMENT_UL = 1   # Minimum if included (can be 0)
MAX_SUPPLEMENT_UL = 90  # Maximum per supplement
MIN_NOVEL_BIO_UL = 90   # Minimum Novel_Bio volume

# Algorithm parameters
DELTA_UL = 10            # Perturbation step size
ALPHA = 1.0              # Learning rate (multiplied by delta for step size)
MAX_ITERATIONS = 8       # Max iterations (one per row, A-H)
CONVERGENCE_ROUNDS = 2   # Stop if no improvement for this many consecutive rounds

# Starting point (iteration 1 center)
INITIAL_COMPOSITION = {
    "Glucose": 20,
    "NaCl": 20,
    "MgSO4": 20,
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
) -> list:
    """Generate a transfer array for one iteration (11 wells in 1 row).

    Returns: [[source_well, dest_well, volume_uL], ...]

    Column layout (all in the same row):
      Col 1:  Seed well (pre-loaded, not part of transfer array)
      Col 2:  Negative control (200 uL Novel_Bio, no cells)
      Col 3:  Positive control (180 uL Novel_Bio, cells added by seeding)
      Col 4:  Center point (current best recipe)
      Col 5-6:  +delta Glucose (2 reps)
      Col 7-8:  +delta NaCl (2 reps)
      Col 9-10: +delta MgSO4 (2 reps)
      Col 11-12: Extra wells (same as center for now)
    """
    transfers = []
    row = row_letter

    # Col 2: Negative control — 200 uL Novel_Bio (no cells seeded)
    transfers.append([REAGENT_WELLS["Novel_Bio"], f"{row}2", WELL_VOLUME_UL])

    # Col 3: Positive control — 180 uL Novel_Bio (cells added by seeding step)
    transfers.append([REAGENT_WELLS["Novel_Bio"], f"{row}3", REAGENT_VOLUME_UL])

    # Col 4: Center point
    novel_bio_center = compute_novel_bio(center)
    if novel_bio_center > 0:
        transfers.append([REAGENT_WELLS["Novel_Bio"], f"{row}4", novel_bio_center])
    for name in SUPPLEMENT_NAMES:
        if center[name] > 0:
            transfers.append([REAGENT_WELLS[name], f"{row}4", center[name]])

    # Cols 5-10: Perturbations (2 reps each for 3 supplements)
    for col1, col2, supplement in PERTURBATION_COLS:
        perturbed = make_perturbed(center, supplement, delta)
        novel_bio_pert = compute_novel_bio(perturbed)

        for col in (col1, col2):
            if novel_bio_pert > 0:
                transfers.append([REAGENT_WELLS["Novel_Bio"], f"{row}{col}", novel_bio_pert])
            for name in SUPPLEMENT_NAMES:
                if perturbed[name] > 0:
                    transfers.append([REAGENT_WELLS[name], f"{row}{col}", perturbed[name]])

    # Cols 11-12: Extra wells (duplicate center for additional data)
    for col in (11, 12):
        if novel_bio_center > 0:
            transfers.append([REAGENT_WELLS["Novel_Bio"], f"{row}{col}", novel_bio_center])
        for name in SUPPLEMENT_NAMES:
            if center[name] > 0:
                transfers.append([REAGENT_WELLS[name], f"{row}{col}", center[name]])

    # Sort: Novel_Bio (D1) first, then MgSO4 (C1), NaCl (B1), Glucose (A1)
    source_order = [
        REAGENT_WELLS["Novel_Bio"],   # D1 — reuse tip, bottom dispense
        REAGENT_WELLS["MgSO4"],       # C1
        REAGENT_WELLS["NaCl"],        # B1
        REAGENT_WELLS["Glucose"],     # A1
    ]
    order_map = {well: i for i, well in enumerate(source_order)}
    transfers.sort(key=lambda t: order_map.get(t[0], 99))

    return transfers


# Novel_Bio well is the only one that reuses tips
REUSE_TIP_SOURCE_WELLS = [REAGENT_WELLS["Novel_Bio"]]


def compute_tip_consumption(transfer_array: list) -> dict:
    """Compute tip consumption for hybrid tip mode.

    Novel_Bio (D1): 1 tip reused for all transfers (grouped by pipette type).
    All other source wells: 1 new tip per transfer.
    """
    reuse_set = set(REUSE_TIP_SOURCE_WELLS)

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
    Cells are seeded into cols 3-12 (skip col 2 = negative control).

    Iteration 1: seed from A1, experiments in row A (A2-A12), warm up B1
    Iteration 8: seed from H1, experiments in row H (H2-H12), no warmup
    """
    row = ROWS[iteration - 1]
    seed_well = f"{row}1"
    is_last = iteration >= MAX_ITERATIONS
    next_seed_well = f"{ROWS[iteration]}1" if not is_last else "B1"  # dummy for last
    nm_cells_volume = 0 if is_last else 220

    # Wells to seed with cells: cols 3-12 (skip col 2 = neg control)
    seed_dest_wells = [f"{row}{col}" for col in range(3, 13)]

    return {
        "seed_well": seed_well,
        "next_seed_well": next_seed_well,
        "nm_cells_volume": nm_cells_volume,
        "seed_dest_wells": seed_dest_wells,
    }


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
) -> Path:
    """Write an iteration-specific workflow definition file.

    Copies the template and replaces iteration-specific constants using regex
    to handle any default values in the template.
    """
    template = WORKFLOW_TEMPLATE_PATH.read_text()

    tip_counts = compute_tip_consumption(transfer_array)
    # Reagent well counter: the 24-well reagent plate tracks usage via this
    # count. Each iteration "consumes" N wells — the plate allows 8 uses
    # (24 wells / 3 supplements = 8). Only count supplement wells (exclude
    # Novel_Bio D1 which is not a consumed reagent slot).
    supplement_wells_used = len(
        set(t[0] for t in transfer_array) - {REAGENT_WELLS["Novel_Bio"]}
    )

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
    p50_extra = 1 if seed_params else 0   # +1 P50 for seeding (reused for 10 wells)
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


def fetch_absorbance_results(plate_barcode: str, row_letter: str) -> dict:
    """Fetch OD600 readings for a plate: both baseline (earliest) and endpoint (latest).

    Queries the datasets REST API (camelCase response), filters by plate UUID
    and OD600 wavelength, and extracts well values for the target row.

    Returns: {
        "baseline": {well: od600_value},   # First reading (pre-growth)
        "endpoint": {well: od600_value},   # Latest reading (post-growth)
    }
    """
    plate_uuid = _get_plate_uuid(plate_barcode)

    # Fetch all datasets, ordered newest first
    resp = requests.get(
        f"{WORKCELL_API_BASE}/api/datasets/",
        headers=API_HEADERS,
        params={"verbose": "1", "ordering": "-createdAt"},
        timeout=30,
    )
    resp.raise_for_status()
    data = resp.json()
    datasets = data.get("results", data) if isinstance(data, dict) else data

    # Filter for OD600 datasets matching our plate UUID
    absorbance_datasets = []
    for ds in datasets:
        meta = ds.get("metadata", {})
        rm = meta.get("resultMetadata", {})
        pm = meta.get("plateMetadata", {})
        if rm.get("measurementWavelength") == 600 and pm.get("uuid") == plate_uuid:
            absorbance_datasets.append(ds)

    if not absorbance_datasets:
        raise RuntimeError(f"No OD600 datasets found for plate {plate_barcode}")

    # Collect timestamps that have data for our target row wells.
    # Pre-absorbance readings cover prior rows only (current row is empty),
    # so we filter to timestamps where at least one target well has a nonzero value.
    target_wells = [f"{row_letter}{col}" for col in range(2, 13)]
    row_readings = {}
    for ds in absorbance_datasets:
        sd = ds.get("structuredData", {})
        results = sd.get("resultsByWell", {})
        for ts, wells in results.items():
            if any(wells.get(w) for w in target_wells):
                row_readings[ts] = wells

    if not row_readings:
        raise RuntimeError(
            f"No OD600 readings found for row {row_letter} wells on plate {plate_barcode}"
        )

    sorted_timestamps = sorted(row_readings.keys())
    earliest_well_data = row_readings[sorted_timestamps[0]]
    latest_well_data = row_readings[sorted_timestamps[-1]]

    print(f"    Row {row_letter} readings: {sorted_timestamps[0]} -> {sorted_timestamps[-1]} "
          f"({len(sorted_timestamps)} readings)")

    # Extract the 11 wells in our target row (cols 2-12) for both timepoints
    baseline = {}
    endpoint = {}
    for col in range(2, 13):
        well = f"{row_letter}{col}"
        baseline[well] = float(earliest_well_data.get(well, 0.0))
        endpoint[well] = float(latest_well_data.get(well, 0.0))

    return {"baseline": baseline, "endpoint": endpoint}


# =============================================================================
# RESULT PARSING
# =============================================================================


def parse_od_results(absorbance_results: dict, row_letter: str) -> dict:
    """Parse baseline + endpoint absorbance into growth deltas for gradient computation.

    Uses delta OD (endpoint - baseline) as the objective function, which normalizes
    for varying initial cell densities across wells.

    Row-wise layout:
      Col 2: Negative control (no cells)
      Col 3: Positive control (cells in plain media)
      Col 4: Center point
      Cols 5-10: Perturbation wells (3 supplements x 2 reps)
      Cols 11-12: Extra wells

    Returns:
        {
            "neg_control_od": float,   # delta OD for negative control (no cells)
            "control_od": float,       # delta OD for positive control
            "center_od": float,        # delta OD for center well
            "perturbed_ods": {supplement: [rep1_delta, rep2_delta]},
            "extra_ods": [extra1_delta, extra2_delta],
            "abs_control_od": float,   # absolute endpoint (for logging)
            "abs_center_od": float,    # absolute endpoint (for logging)
        }
    """
    row = row_letter
    baseline = absorbance_results.get("baseline", {})
    endpoint = absorbance_results.get("endpoint", {})

    def delta(well):
        return endpoint.get(well, 0.0) - baseline.get(well, 0.0)

    neg_control_delta = delta(f"{row}2")
    pos_control_delta = delta(f"{row}3")
    center_delta = delta(f"{row}4")

    perturbed_deltas = {}
    for col1, col2, supplement in PERTURBATION_COLS:
        perturbed_deltas[supplement] = [delta(f"{row}{col1}"), delta(f"{row}{col2}")]

    extra_deltas = [delta(f"{row}11"), delta(f"{row}12")]

    return {
        "neg_control_od": neg_control_delta,
        "control_od": pos_control_delta,
        "center_od": center_delta,
        "perturbed_ods": perturbed_deltas,
        "extra_ods": extra_deltas,
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
        # Step 1: Generate transfer array
        # Perturbation delta scales with alpha so we test the actual step size
        perturbation_delta = max(1, int(alpha * DELTA_UL))
        transfer_array = generate_transfer_array(center, row_letter, delta=perturbation_delta)

        # Step 2: Write workflow definition
        workflow_path = write_workflow_definition(
            transfer_array, row_letter, iteration, seed_params
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
            alpha = max(0.5, alpha / 2)
            print(f"\n  Center OD decreased ({state['prev_center_od']:.4f} -> "
                  f"{od_results['center_od']:.4f}), halving alpha to {alpha}")

    # Step 9: Update composition
    new_composition = gradient_step(center, gradient, alpha)

    # Log full iteration with results
    log_iteration(iteration, center, transfer_array, od_results, gradient, new_composition, alpha)

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
    print(f"  {'Iter':>4} | {'Novel_Bio':>9} | {'Glucose':>7} | {'NaCl':>5} | {'MgSO4':>6} | {'Center OD':>9} | {'Control':>7}")
    print(f"  {'----':>4} | {'---------':>9} | {'-------':>7} | {'-----':>5} | {'------':>6} | {'---------':>9} | {'-------':>7}")
    for h in state["history"]:
        c = h["composition"]
        print(f"  {h['iteration']:4d} | {c['Novel_Bio']:9d} | {c['Glucose']:7d} | "
              f"{c['NaCl']:5d} | {c['MgSO4']:6d} | {h['center_od']:9.4f} | {h['control_od']:7.4f}")


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
            print(f"    Iter {h['iteration']}: NB={c['Novel_Bio']}, "
                  f"Glc={c['Glucose']}, NaCl={c['NaCl']}, Mg={c['MgSO4']} "
                  f"-> OD={h['center_od']:.4f}")


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
