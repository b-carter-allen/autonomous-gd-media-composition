from src.platform.core_domain.units import Time
from src.workflows.workflow_definition_dsl.workflow_definition_descriptor import (
    MoreThanConstraint,
    RoutineReference,
    WorkflowDefinitionDescriptor,
)

# =============================================================================
# ITERATION-SPECIFIC CONSTANTS
# These values are baked in per iteration by the gradient descent orchestrator.
# The orchestrator generates a new copy of this file each iteration
# with updated values below, then registers and instantiates it.
# =============================================================================

# Transfer array: JSON string of [[source_well, dest_well, volume_uL], ...]
# Source wells are on the reagent plate, dest wells on the experiment plate.
# Generated from 8 wells (1 column) x 4 reagents per iteration.
TRANSFER_ARRAY = "[]"

# Target experimental column on the 96-well plate (2-9)
DEST_COLUMN_INDEX = 2

# Seed well on experiment plate for this round (A1=round 1, B1=round 2, ...)
SEED_WELL = "A1"

# Next seed well to warm up with NM+Cells (B1=round 1, C1=round 2, ...)
# Ignored when NM_CELLS_VOLUME = 0 (round 8)
NEXT_SEED_WELL = "B1"

# NM+Cells source well on reagent plate (24-well deep well)
NM_CELLS_SOURCE_WELL = "A2"

# Volume of NM+Cells to transfer to next seed well (uL). 0 to skip (round 8).
NM_CELLS_VOLUME = 220

# Seed transfer volume from seed well to each experimental well (uL)
SEED_TRANSFER_VOLUME = 20

# Seed well mix parameters
MIX_VOLUME = 100  # uL
MIX_REPS = 5

# Reagent plate type tag (matches reagent plate tags in the system)
REAGENT_TYPE = "GD Compound Stock Plate"

# Tip consumption for combined routine
# Reagent transfers: varies by source grouping (reuse_tips_for_same_source)
# + 1 P200 for seed well mixing
# + 1 P50 for seeding 8 wells (reused)
# + 1 P1000 for NM+Cells warmup (if nm_cells_volume > 0)
P50_TIPS_TO_CONSUME = 5  # 4 reagent + 1 seeding
P200_TIPS_TO_CONSUME = 1  # 1 seed well mixing
P1000_TIPS_TO_CONSUME = 1  # 1 NM+Cells warmup (0 for round 8)

# Number of reagent wells consumed (4 reagents + 1 NM+Cells = 5)
REAGENT_WELLS_TO_CONSUME = 5

# Reuse tips for consecutive transfers from same source well
REUSE_TIPS_FOR_SAME_SOURCE = True

# OD600 monitoring: 5 min intervals for 1.5 hours (18 readings)
MONITORING_INTERVAL_MINUTES = 5
MONITORING_READINGS = 18


def build_definition(plate_barcode: str = "default_plate") -> WorkflowDefinitionDescriptor:
    """Gradient Descent Media Optimization — Single Iteration Workflow.

    One complete iteration of the gradient descent experiment:
      Phase 1: Pre-iteration absorbance of all filled wells
      Phase 2: Combined liquid handling (reagent transfers + on-plate seeding + NM warmup)
      Phase 3: Monitor OD600 absorbance at 5 min intervals for 1.5 hours (18 readings)

    Only plate_barcode is passed at instantiation time. All other parameters
    are set as module-level constants by the orchestration script before
    registering this workflow definition.
    """
    workflow = WorkflowDefinitionDescriptor(
        description=(
            f"Gradient Descent Media Optimization: Iteration column {DEST_COLUMN_INDEX}. "
            "Combined reagent transfer + on-plate seeding + NM warmup, "
            "then monitors growth via OD600 absorbance for 1.5 hours."
        ),
    )

    rows = ["A", "B", "C", "D", "E", "F", "G", "H"]

    # -------------------------------------------------------------------------
    # Compute well lists for absorbance measurements
    # -------------------------------------------------------------------------
    # Round number: column 2 = round 1, column 3 = round 2, ...
    round_number = DEST_COLUMN_INDEX - 1

    # Pre-iteration wells: seed wells + prior experimental columns
    pre_wells = []
    for i in range(round_number):
        pre_wells.append(f"{rows[i]}1")
    for col in range(2, DEST_COLUMN_INDEX):
        for row in rows:
            pre_wells.append(f"{row}{col}")

    # Post-iteration wells: all seed wells (incl. next) + all exp columns (incl. current)
    post_wells = []
    seed_count = round_number + (1 if NM_CELLS_VOLUME > 0 else 0)
    for i in range(seed_count):
        post_wells.append(f"{rows[i]}1")
    for col in range(2, DEST_COLUMN_INDEX + 1):
        for row in rows:
            post_wells.append(f"{row}{col}")

    # === PHASE 1: PRE-ITERATION ABSORBANCE ===
    if pre_wells:
        pre_absorbance = RoutineReference(
            routine_name="Measure Absorbance",
            routine_parameters={
                "culture_plate_barcode": plate_barcode,
                "method_name": "96wp_od600",
                "wells_to_process": pre_wells,
            },
        )
        workflow.add_routine("pre_absorbance", pre_absorbance)

    # === PHASE 2: GD ITERATION COMBINED ===
    # Single OT Flex routine: reagent transfers + seed from on-plate well + NM warmup
    gd_iteration = RoutineReference(
        routine_name="GD Iteration Combined",
        routine_parameters={
            "experiment_plate_barcode": plate_barcode,
            "reagent_type": REAGENT_TYPE,
            "transfer_array": TRANSFER_ARRAY,
            "seed_well": SEED_WELL,
            "dest_column_index": DEST_COLUMN_INDEX,
            "seed_transfer_volume": SEED_TRANSFER_VOLUME,
            "nm_cells_source_well": NM_CELLS_SOURCE_WELL,
            "nm_cells_volume": NM_CELLS_VOLUME,
            "next_seed_well": NEXT_SEED_WELL,
            "mix_volume": MIX_VOLUME,
            "mix_reps": MIX_REPS,
            "p50_tips_to_consume": P50_TIPS_TO_CONSUME,
            "p200_tips_to_consume": P200_TIPS_TO_CONSUME,
            "p1000_tips_to_consume": P1000_TIPS_TO_CONSUME,
            "reuse_tips_for_same_source": REUSE_TIPS_FOR_SAME_SOURCE,
            "reagent_wells_to_consume": REAGENT_WELLS_TO_CONSUME,
        },
    )
    workflow.add_routine("gd_iteration", gd_iteration)

    # Iteration starts after pre-absorbance (if it exists)
    if pre_wells:
        workflow.add_time_constraint(
            MoreThanConstraint(
                from_start="pre_absorbance",
                to_start="gd_iteration",
                value=Time("0 minutes"),
            )
        )

    # === PHASE 3: OD600 MONITORING ===
    # Read absorbance at 5-minute intervals for 1.5 hours (18 readings)
    # Reads ALL filled wells (seed wells + experimental columns)
    monitoring_absorbance = RoutineReference(
        routine_name="Measure Absorbance",
        routine_parameters={
            "culture_plate_barcode": plate_barcode,
            "method_name": "96wp_od600",
            "wells_to_process": post_wells,
        },
    )

    monitoring_routines = []
    for i in range(MONITORING_READINGS):
        routine_key = f"od600_reading_{i + 1}"
        workflow.add_routine(routine_key, monitoring_absorbance)
        monitoring_routines.append(routine_key)

    # Space out readings by the monitoring interval
    workflow.space_out_routines(monitoring_routines, Time(f"{MONITORING_INTERVAL_MINUTES} minutes"))

    # First monitoring reading starts 30s after iteration completes
    workflow.add_time_constraint(
        MoreThanConstraint(
            from_start="gd_iteration",
            to_start=monitoring_routines[0],
            value=Time("30 seconds"),
        )
    )

    return workflow
