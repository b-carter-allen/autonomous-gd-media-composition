# Autonomous Gradient Descent of Existing Media Composition

An autonomous optimizer that tunes growth media composition for *Vibrio natriegens* using gradient descent on a [Monomer Bio](https://monomerbio.com) robotic workcell. The daemon designs media recipes, the robot executes them, OD600 growth is measured, and the algorithm iterates — no human in the loop.

## Current Experiment

**Hypothesis**: Novel_Bio base media may be oversalted for optimal *V. natriegens* growth. We supplement with Glucose (carbon source), MOPS (pH buffer), and DiH2O (to dilute excess salt) and let the optimizer find the best concentrations.

**Supplements**:
| Reagent | Well | Concentration |
|---------|------|---------------|
| Glucose | A1 | 100 mg/mL |
| MOPS | B1 | 1 M |
| DiH2O | C1 | Pure |
| Novel_Bio | D1 | Base media |

**Starting composition**: 20 uL each supplement, 120 uL Novel_Bio, 20 uL cells = 200 uL total.

## How It Works

Each iteration uses one row (12 wells) of a 96-well plate. Column 1 holds the seed well (cells stock), and columns 2-12 are experimental wells:

| Column | Role | Composition |
|--------|------|-------------|
| 1 | Seed well | NM+Cells stock (pre-loaded) |
| 2 | Negative control | 200 uL Novel_Bio, no cells |
| 3 | Positive control | 180 uL Novel_Bio + 20 uL cells |
| 4 | Center | Current best recipe + 20 uL cells |
| 5-6 | +Glucose | Perturbation wells (2 reps) + cells |
| 7-8 | +MOPS | Perturbation wells (2 reps) + cells |
| 9-10 | +DiH2O | Perturbation wells (2 reps) + cells |
| 11-12 | Extra | Additional sample wells + cells |

The daemon loop:

1. Generates a transfer array from the current composition + perturbations
2. Writes and uploads a workflow definition to the workcell via MCP
3. Instantiates the workflow (compound plate generation, seeding, OD600 monitoring)
4. Polls for completion (~90 min per iteration)
5. Fetches OD600 results from the datasets REST API
6. Computes delta OD (endpoint - baseline) and the gradient
7. Steps the composition in the gradient direction
8. Repeats for up to 8 iterations (one per row, A-H)

## Workflow Diagram

```mermaid
flowchart TD
    START([Start Daemon]) --> INIT[Load state\ncurrent_composition\niteration count]
    INIT --> CHECK{iteration ≤ 8\nand not converged?}
    CHECK -- No --> DONE([Experiment Complete\nOptimal composition found])

    CHECK -- Yes --> ROW["Select row\n(A=iter 1 · · · H=iter 8)"]
    ROW --> COMPOSE[Generate transfer array\nfrom center + perturbations]
    COMPOSE --> UPLOAD[Register workflow definition\nvia MCP]
    UPLOAD --> INST[Instantiate workflow\non workcell]

    INST --> WF

    subgraph WF["  Workcell Workflow — ~95 min per iteration  "]
        direction TB
        P1["📊 Phase 1 — Pre-iteration Absorbance\nBaseline OD600 of seed wells + prior rows"]
        P2["🤖 Phase 2 — GD Iteration Combined  (OT Flex)\n① Remove lids\n② Transfer reagents: Novel_Bio → cols 2–12\n   Supplements (Glc/MOPS/DiH₂O) → cols 4–12\n③ Mix seed well, seed cols 3–12 (20 µL each)\n④ Warm up next seed well with NM+Cells\n⑤ Replace lids"]
        P3["📈 Phase 3 — OD600 Growth Monitoring\n18 readings × 5 min intervals = 90 min"]
        P1 --> P2 --> P3
    end

    WF --> FETCH[Fetch OD600 results\nfrom REST API]
    FETCH --> DELTA[Compute Δ OD\nendpoint − baseline per well]
    DELTA --> GRAD["Compute gradient\n∂(growth)/∂(supplement) via finite differences"]
    GRAD --> CONV{Center better\nthan last round?}
    CONV -- No --> DECAY[Decay learning rate α\nα = max(α/2, 0.1)]
    CONV -- Yes --> STEP
    DECAY --> STEP[Step composition\nstep = α × Δ × sign(gradient)]
    STEP --> SAVE[Save state to disk]
    SAVE --> CHECK
```

### Plate Layout per Iteration

```
     Col 1      Col 2       Col 3      Col 4     Col 5-6    Col 7-8    Col 9-10   Col 11-12
Row A ║ Seed A1 ║ Neg Ctrl  ║ Pos Ctrl ║ Center  ║ +Glucose ║ +MOPS    ║ +DiH2O   ║ Extra
Row B ║ Seed B1 ║ Neg Ctrl  ║ Pos Ctrl ║ Center  ║ +Glucose ║ +MOPS    ║ +DiH2O   ║ Extra
 ...  ║  ...    ║   ...     ║   ...    ║  ...    ║   ...    ║  ...     ║   ...    ║  ...
Row H ║ Seed H1 ║ Neg Ctrl  ║ Pos Ctrl ║ Center  ║ +Glucose ║ +MOPS    ║ +DiH2O   ║ Extra

Seed col 1: pre-loaded NM+Cells (next row warmed up each iteration)
Neg ctrl:   200 µL Novel_Bio only, no cells — growth floor reference
Pos ctrl:   180 µL Novel_Bio + 20 µL cells — unoptimized growth reference
Center:     Current best composition + cells
Perturbations: center ± 10 µL on one supplement (2 replicates each)
```

## Repo Structure

```
autonomous-gd-media-composition/
├── experimental_design/          # Define WHAT to test
│   ├── README.md                 # Template for hypothesis, reagents, parameters
│   └── EXPERIMENT_INFO.md        # Current experiment details
│
├── experimental_execution/       # HOW to run the experiment
│   ├── gradient_descent.py       # Autonomous daemon
│   ├── dashboard.py              # Streamlit visualization
│   ├── workflow_template.py      # Workcell workflow definition
│   └── simulate_data.py          # Generate test data for dashboard
│
└── data/                         # Created at runtime (gitignored)
```

**`experimental_design/`** — Start here. Define your hypothesis, choose your reagents and starting composition, and set algorithm parameters. See `EXPERIMENT_INFO.md` for the current experiment.

**`experimental_execution/`** — The autonomous system. Once your design is set, the daemon handles everything: pipetting, incubation, measurement, and optimization.

## Setup

### Prerequisites

- Python 3.11+
- A Monomer Bio workcell on the local network
- A 96-well culture plate loaded on the workcell
- A reagent plate with your supplements and base media

### Install

```bash
git clone https://github.com/b-carter-allen/autonomous-gd-media-composition.git
cd autonomous-gd-media-composition

python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

cp .env.example .env
# Edit .env with your workcell IP address
```

### Physical Setup

Before running, prepare:

1. **Culture plate** — 96-well plate with a barcode, loaded on the workcell
2. **Reagent plate** — 24-well deep well plate (type: "AGD Stock Plate") with:
   - **A1**: Glucose 100 mg/mL (5 mL)
   - **B1**: MOPS 1 M (5 mL)
   - **C1**: DiH2O (5 mL)
   - **D1**: Novel_Bio base media (9 mL)
   - **D2**: Novel_Bio base media (9 mL, overflow)
   - **A2**: NM+Cells (5 mL, refrigerated)
3. **Tip racks** — ~200 P200 tips, ~60 P50 tips across all iterations

## Running

```bash
cd experimental_execution

# Check workcell connectivity
python gradient_descent.py status

# Preview the first iteration without running anything
python gradient_descent.py dry-run <plate_barcode>

# Start a fresh experiment
python gradient_descent.py run <plate_barcode>

# Resume after a restart (picks up mid-iteration)
python gradient_descent.py resume <plate_barcode>

# Reset state for a new experiment
python gradient_descent.py reset
```

The daemon runs in the foreground and logs progress to stdout. For unattended operation:

```bash
PYTHONUNBUFFERED=1 nohup python gradient_descent.py run <plate_barcode> > gd.log 2>&1 &
```

### Simulate Data

Generate test data for the dashboard without a workcell:

```bash
cd experimental_execution
python simulate_data.py [plate_barcode]
```

## Dashboard

A Streamlit dashboard provides live visualization of experiment progress:

```bash
cd experimental_execution
streamlit run dashboard.py
```

Features:
- 96-well plate heatmap with delta OD color scale
- OD600 growth progress across iterations
- Composition trajectory over time
- Gradient direction indicators per supplement
- Live workflow progress when connected to the workcell
- Auto-refresh mode

## Algorithm

**Sign-based gradient ascent** with adaptive learning rate:

- `step = int(alpha * delta) * sign(gradient)`
- Perturbation delta scales with alpha: `perturbation = max(1, int(alpha * DELTA_UL))`
- Alpha starts at 1.0 and halves when center growth decreases (floor: 0.1)
- Step sizes narrow as alpha decays: 10 → 5 → 2 → 1 uL
- Objective function: delta OD600 (endpoint - baseline), normalizing for initial cell density

All tunable parameters are at the top of `gradient_descent.py`:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `DELTA_UL` | 10 | Base perturbation size (uL) |
| `ALPHA` | 1.0 | Initial learning rate |
| `MAX_ITERATIONS` | 8 | Max iterations per plate |
| `REAGENT_VOLUME_UL` | 180 | Reagent mix volume per well (uL, before cells) |
| `SEED_TRANSFER_VOLUME` | 20 | Cells added from seed well (uL) |
| `WELL_VOLUME_UL` | 200 | Total volume per well (reagent + cells) |
| `MIN_SUPPLEMENT_UL` | 5 | Minimum supplement volume if included |
| `MAX_SUPPLEMENT_UL` | 90 | Maximum supplement volume |
| `MIN_NOVEL_BIO_UL` | 90 | Minimum base media volume |
| `MONITORING_INTERVAL_MINUTES` | 5 | OD600 reading interval |
| `MONITORING_READINGS` | 18 | Number of readings (90 min total) |

## Workcell Communication

The daemon uses two protocols to communicate with the Monomer Bio workcell:

- **MCP (Streamable HTTP)** at `/mcp` — Workflow registration, instantiation, status polling, plate operations. JSON-RPC 2.0 over HTTP POST with SSE responses.
- **REST API** at `/api/` — Dataset retrieval (OD600 absorbance results). Requires `X-Monomer-Client: desktop-frontend` header.

No API keys needed — the daemon communicates directly with the workcell over the local network.
