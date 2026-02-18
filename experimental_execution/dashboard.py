"""Gradient Descent Media Optimization — Live Dashboard.

Usage:
    streamlit run dashboard.py
    streamlit run dashboard.py --server.address 0.0.0.0 --server.port 8501
"""

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

# ---------------------------------------------------------------------------
# Constants (mirrored from gradient_descent.py)
# ---------------------------------------------------------------------------
ROWS = ["A", "B", "C", "D", "E", "F", "G", "H"]
COLS = list(range(1, 13))
SUPPLEMENT_NAMES = ["Glucose_100mg_mL", "MOPS_1M", "NaCl_2M"]
REAGENT_VOLUME_UL = 180
WELL_VOLUME_UL = 200
MAX_ITERATIONS = 8

# Row-wise layout: column -> role (each iteration is one row)
COL_ROLES = {
    1: "Seed",
    2: "Neg Control",
    3: "Pos Control",
    4: "Center",
    5: f"{SUPPLEMENT_NAMES[0]} +d",
    6: f"{SUPPLEMENT_NAMES[0]} +d",
    7: f"{SUPPLEMENT_NAMES[1]} +d",
    8: f"{SUPPLEMENT_NAMES[1]} +d",
    9: f"{SUPPLEMENT_NAMES[2]} +d",
    10: f"{SUPPLEMENT_NAMES[2]} +d",
    11: "Extra",
    12: "Extra",
}

ROLE_COLORS = {
    "Neg Control": "#64748b",
    "Pos Control": "#94a3b8",
    "Center": "#3b82f6",
    f"{SUPPLEMENT_NAMES[0]} +d": "#f59e0b",
    f"{SUPPLEMENT_NAMES[1]} +d": "#10b981",
    f"{SUPPLEMENT_NAMES[2]} +d": "#8b5cf6",
    "Extra": "#6b7280",
    "Seed": "#d1d5db",
}

# Reagent source wells on the compound plate
REAGENT_SOURCE = {
    "D1": "Novel_Bio",
    "A1": SUPPLEMENT_NAMES[0],
    "B1": SUPPLEMENT_NAMES[1],
    "C1": SUPPLEMENT_NAMES[2],
}

BASE_DATA_DIR = Path(__file__).parent.parent / "data"

# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------


def load_state() -> dict:
    state_path = DATA_DIR / "state.json"
    if state_path.exists():
        return json.loads(state_path.read_text())
    return {
        "current_iteration": 0,
        "current_composition": {name: 20 for name in SUPPLEMENT_NAMES},
        "alpha": 1.0,
        "best_od": None,
        "prev_center_od": None,
        "no_improvement_count": 0,
        "converged": False,
        "history": [],
    }


def load_iteration_log(iteration: int) -> dict | None:
    path = DATA_DIR / f"iteration_{iteration}" / "iteration_log.json"
    if path.exists():
        return json.loads(path.read_text())
    return None


def load_transfer_array(iteration: int) -> list | None:
    path = DATA_DIR / f"iteration_{iteration}" / "transfer_array.json"
    if path.exists():
        return json.loads(path.read_text())
    return None


def load_workflow_ids(iteration: int) -> dict | None:
    path = DATA_DIR / f"iteration_{iteration}" / "workflow_ids.json"
    if path.exists():
        return json.loads(path.read_text())
    return None


def get_live_workflow_status(workflow_ids: dict) -> dict | None:
    """Try to get live workflow status from workcell MCP."""
    try:
        from gradient_descent import McpClient

        client = McpClient()
        client.connect()
        return client.call_tool(
            "get_workflow_instance_details",
            {"instance_uuid": workflow_ids["workflow_instance_uuid"]},
        )
    except Exception:
        return None


def compute_well_volumes(transfer_array: list) -> dict:
    """Aggregate transfer array into per-well volume breakdowns.

    Returns: {dest_well: {reagent_name: volume, ...}, ...}
    """
    well_volumes = {}
    for src, dest, vol in transfer_array:
        reagent = REAGENT_SOURCE.get(src, src)
        if dest not in well_volumes:
            well_volumes[dest] = {}
        well_volumes[dest][reagent] = well_volumes[dest].get(reagent, 0) + vol
    return well_volumes


def uses_delta_od(history: list) -> bool:
    """Detect whether experiment uses delta OD (has abs_center_od in history)."""
    return any("abs_center_od" in h for h in history)


# ---------------------------------------------------------------------------
# Plate heatmap
# ---------------------------------------------------------------------------


def build_plate_data(state: dict, is_delta: bool) -> tuple[np.ndarray, list[list[str]]]:
    """Build 8x12 arrays for OD values and hover text (row-wise layout)."""
    od_grid = np.full((8, 12), np.nan)
    hover_grid = [["" for _ in range(12)] for _ in range(8)]

    for h in state.get("history", []):
        iteration = h["iteration"]
        row_i = iteration - 1  # iteration 1 = row A (index 0)
        if row_i >= 8:
            continue
        row = ROWS[row_i]

        od_results = h.get("od_results", {})
        if not od_results:
            log = load_iteration_log(iteration)
            if log:
                od_results = log.get("od_results", {})

        # Load transfer array for volume info
        ta = load_transfer_array(iteration)
        well_vols = compute_well_volumes(ta) if ta else {}

        if od_results:
            neg_control_od = od_results.get("neg_control_od", 0)
            control_od = od_results.get("control_od", 0)
            center_od = od_results.get("center_od", 0)
            perturbed = od_results.get("perturbed_ods", {})
            extra_ods = od_results.get("extra_ods", [0, 0])

            cols_od = {
                2: neg_control_od,
                3: control_od,
                4: center_od,
                5: perturbed.get(SUPPLEMENT_NAMES[0], [0, 0])[0],
                6: perturbed.get(SUPPLEMENT_NAMES[0], [0, 0])[1],
                7: perturbed.get(SUPPLEMENT_NAMES[1], [0, 0])[0],
                8: perturbed.get(SUPPLEMENT_NAMES[1], [0, 0])[1],
                9: perturbed.get(SUPPLEMENT_NAMES[2], [0, 0])[0],
                10: perturbed.get(SUPPLEMENT_NAMES[2], [0, 0])[1],
                11: extra_ods[0] if extra_ods else 0,
                12: extra_ods[1] if len(extra_ods) > 1 else 0,
            }

            for col_num, od_val in cols_od.items():
                col_i = col_num - 1  # 0-indexed
                od_grid[row_i][col_i] = od_val
                role = COL_ROLES.get(col_num, "")
                well_name = f"{row}{col_num}"

                # Build volume breakdown for hover
                vols = well_vols.get(well_name, {})
                vol_lines = "".join(
                    f"{name}: {v} uL<br>"
                    for name, v in sorted(vols.items())
                ) if vols else ""

                if is_delta:
                    od_label = f"Growth: {od_val:+.4f}"
                else:
                    od_label = f"OD600: {od_val:.4f}"

                hover_grid[row_i][col_i] = (
                    f"<b>{well_name}</b><br>"
                    f"{od_label}<br>"
                    f"Role: {role}<br>"
                    f"Iteration {iteration}<br>"
                    f"---<br>"
                    f"{vol_lines}"
                    f"Total: {sum(vols.values()):.0f} uL"
                    if vols else
                    f"<b>{well_name}</b><br>"
                    f"{od_label}<br>"
                    f"Role: {role}<br>"
                    f"Iteration {iteration}"
                )

    # Mark seed wells in column 1
    for row_i, row in enumerate(ROWS):
        col_i = 0
        iteration_for_seed = row_i + 1
        if iteration_for_seed <= state.get("current_iteration", 0):
            hover_grid[row_i][col_i] = (
                f"<b>{row}1</b><br>Seed well (iteration {iteration_for_seed})"
            )
            od_grid[row_i][col_i] = -0.01  # marker value
        else:
            hover_grid[row_i][col_i] = f"<b>{row}1</b><br>Empty seed well"

    # Mark unfilled rows
    current_iter = state.get("current_iteration", 0)
    for row_i in range(current_iter, 8):
        for col_i in range(1, 12):
            if not hover_grid[row_i][col_i]:
                hover_grid[row_i][col_i] = f"<b>{ROWS[row_i]}{col_i + 1}</b><br>Empty"

    return od_grid, hover_grid


def plate_heatmap(state: dict, is_delta: bool) -> go.Figure:
    od_grid, hover_grid = build_plate_data(state, is_delta)

    if is_delta:
        # Diverging color scale: red (negative growth) -> white (zero) -> green (positive)
        colorscale = [
            [0.0, "#ef4444"],
            [0.3, "#fca5a5"],
            [0.45, "#ffffff"],
            [0.55, "#bbf7d0"],
            [0.7, "#86efac"],
            [1.0, "#15803d"],
        ]
        # Symmetric range around 0
        max_abs = max(abs(np.nanmin(od_grid)), abs(np.nanmax(od_grid)), 0.01)
        zmin, zmax = -max_abs, max_abs
        colorbar_title = "Growth (dOD)"
    else:
        colorscale = [
            [0.0, "#dbeafe"],
            [0.05, "#ffffff"],
            [0.5, "#86efac"],
            [1.0, "#15803d"],
        ]
        zmin, zmax = -0.02, None
        colorbar_title = "OD600"

    fig = go.Figure(
        data=go.Heatmap(
            z=od_grid,
            x=[str(c) for c in COLS],
            y=ROWS,
            text=hover_grid,
            hoverinfo="text",
            colorscale=colorscale,
            zmin=zmin,
            zmax=zmax,
            colorbar=dict(title=colorbar_title, thickness=15),
            xgap=3,
            ygap=3,
        )
    )

    current_iter = state.get("current_iteration", 0)
    annotations = []
    annotations.append(dict(x="1", y=-0.7, text="Seed", showarrow=False, font=dict(size=10, color="#64748b")))
    for i in range(1, min(current_iter + 1, 9)):
        annotations.append(
            dict(x=str(i + 1), y=-0.7, text=f"Iter {i}", showarrow=False, font=dict(size=10, color="#64748b"))
        )

    fig.update_layout(
        height=320,
        margin=dict(l=30, r=30, t=10, b=40),
        xaxis=dict(title=None, side="top", dtick=1),
        yaxis=dict(title=None, autorange="reversed"),
        annotations=annotations,
        plot_bgcolor="#f8fafc",
    )

    return fig


# ---------------------------------------------------------------------------
# Charts
# ---------------------------------------------------------------------------


def od_progress_chart(history: list, is_delta: bool) -> go.Figure:
    if not history:
        fig = go.Figure()
        fig.update_layout(height=300, margin=dict(l=40, r=20, t=30, b=40))
        fig.add_annotation(text="No data yet", xref="paper", yref="paper", x=0.5, y=0.5, showarrow=False)
        return fig

    iters = [h["iteration"] for h in history]
    center = [h.get("center_od", 0) for h in history]
    control = [h.get("control_od", 0) for h in history]

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=iters, y=control, mode="lines+markers", name="Control", line=dict(color="#94a3b8", dash="dash"), marker=dict(size=6)))
    fig.add_trace(go.Scatter(x=iters, y=center, mode="lines+markers", name="Center", line=dict(color="#3b82f6", width=3), marker=dict(size=8)))

    # Perturbation means
    colors = {SUPPLEMENT_NAMES[0]: "#f59e0b", SUPPLEMENT_NAMES[1]: "#10b981", SUPPLEMENT_NAMES[2]: "#8b5cf6"}
    for name in SUPPLEMENT_NAMES:
        means = []
        for h in history:
            od_results = h.get("od_results", {})
            if not od_results:
                log = load_iteration_log(h["iteration"])
                od_results = log.get("od_results", {}) if log else {}
            reps = od_results.get("perturbed_ods", {}).get(name, [0, 0])
            means.append(sum(reps) / max(len(reps), 1))
        fig.add_trace(go.Scatter(x=iters, y=means, mode="lines+markers", name=f"+{name}", line=dict(color=colors[name], width=1.5), marker=dict(size=5)))

    y_label = "Growth (dOD)" if is_delta else "OD600"
    fig.update_layout(
        height=300,
        margin=dict(l=40, r=20, t=30, b=40),
        xaxis=dict(title="Iteration", dtick=1),
        yaxis=dict(title=y_label),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        plot_bgcolor="#f8fafc",
    )

    # Add zero line for delta OD
    if is_delta:
        fig.add_hline(y=0, line_dash="dot", line_color="#cbd5e1", line_width=1)

    return fig


def composition_trajectory_chart(history: list) -> go.Figure:
    if not history:
        fig = go.Figure()
        fig.update_layout(height=300, margin=dict(l=40, r=20, t=30, b=40))
        fig.add_annotation(text="No data yet", xref="paper", yref="paper", x=0.5, y=0.5, showarrow=False)
        return fig

    iters = [h["iteration"] for h in history]
    colors = {"Novel_Bio": "#64748b", SUPPLEMENT_NAMES[0]: "#f59e0b", SUPPLEMENT_NAMES[1]: "#10b981", SUPPLEMENT_NAMES[2]: "#8b5cf6"}

    fig = go.Figure()
    for component in ["Novel_Bio"] + SUPPLEMENT_NAMES:
        vals = [h["composition"].get(component, 0) for h in history]
        fig.add_trace(go.Bar(x=iters, y=vals, name=component, marker_color=colors[component]))

    fig.update_layout(
        barmode="stack",
        height=300,
        margin=dict(l=40, r=20, t=30, b=40),
        xaxis=dict(title="Iteration", dtick=1),
        yaxis=dict(title="Volume (uL)", range=[0, 200]),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        plot_bgcolor="#f8fafc",
    )
    return fig


def objective_chart(history: list, is_delta: bool) -> go.Figure:
    """Hero chart: objective function value vs iteration."""
    if not history:
        fig = go.Figure()
        fig.update_layout(height=300, margin=dict(l=40, r=20, t=30, b=40))
        fig.add_annotation(text="No data yet", xref="paper", yref="paper", x=0.5, y=0.5, showarrow=False)
        return fig

    iters = [h["iteration"] for h in history]
    scores = [h.get("center_od", 0) for h in history]

    # Running best
    best_so_far = []
    current_best = 0
    for s in scores:
        current_best = max(current_best, s)
        best_so_far.append(current_best)

    if is_delta:
        trace_name = "Growth (Center dOD)"
        y_title = "Growth (dOD)"
    else:
        trace_name = "Score (Center OD)"
        y_title = "Objective Score (OD600)"

    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=iters, y=scores, mode="lines+markers",
        name=trace_name,
        line=dict(color="#3b82f6", width=3),
        marker=dict(size=10, symbol="circle"),
    ))
    fig.add_trace(go.Scatter(
        x=iters, y=best_so_far, mode="lines",
        name="Best so far",
        line=dict(color="#15803d", width=2, dash="dot"),
    ))

    fig.update_layout(
        height=300,
        margin=dict(l=40, r=20, t=30, b=40),
        xaxis=dict(title="Iteration", dtick=1, range=[0.5, MAX_ITERATIONS + 0.5]),
        yaxis=dict(title=y_title),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        plot_bgcolor="#f8fafc",
    )

    if is_delta:
        fig.add_hline(y=0, line_dash="dot", line_color="#cbd5e1", line_width=1)

    return fig


def composition_bar(composition: dict) -> go.Figure:
    novel_bio = WELL_VOLUME_UL - sum(composition.get(n, 0) for n in SUPPLEMENT_NAMES)
    components = ["Novel_Bio"] + SUPPLEMENT_NAMES
    values = [novel_bio] + [composition.get(n, 0) for n in SUPPLEMENT_NAMES]
    colors = ["#64748b", "#f59e0b", "#10b981", "#8b5cf6"]

    fig = go.Figure(
        go.Bar(
            y=components,
            x=values,
            orientation="h",
            marker_color=colors,
            text=[f"{v} uL" for v in values],
            textposition="inside",
            textfont=dict(color="white", size=14),
        )
    )
    fig.update_layout(
        height=260,
        margin=dict(l=10, r=20, t=10, b=10),
        xaxis=dict(title="Volume (uL)", range=[0, 200]),
        yaxis=dict(autorange="reversed"),
        plot_bgcolor="#f8fafc",
    )
    return fig


# ---------------------------------------------------------------------------
# Main app
# ---------------------------------------------------------------------------

st.set_page_config(page_title="GD Media Optimization", layout="wide", page_icon="🧫")

# Custom CSS
st.markdown(
    """
    <style>
    .status-pill {
        display: inline-block;
        padding: 4px 16px;
        border-radius: 20px;
        font-weight: 600;
        font-size: 14px;
    }
    .status-running { background: #dbeafe; color: #1d4ed8; }
    .status-converged { background: #dcfce7; color: #15803d; }
    .status-idle { background: #f1f5f9; color: #64748b; }
    .metric-card {
        background: #f8fafc;
        border-radius: 12px;
        padding: 16px;
        text-align: center;
    }
    .metric-value { font-size: 28px; font-weight: 700; color: #1e293b; }
    .metric-label { font-size: 12px; color: #64748b; text-transform: uppercase; letter-spacing: 0.05em; }
    .explanation {
        font-size: 13px;
        color: #64748b;
        margin-top: -8px;
        margin-bottom: 12px;
        line-height: 1.4;
    }
    </style>
    """,
    unsafe_allow_html=True,
)

# ---------------------------------------------------------------------------
# Experiment selector — discover all experiment directories under data/
# ---------------------------------------------------------------------------
def discover_experiments() -> list[str]:
    """Find experiment data dirs (any dir with a state.json)."""
    experiments = []
    for d in sorted(BASE_DATA_DIR.iterdir()):
        if d.is_dir() and (d / "state.json").exists():
            experiments.append(d.name)
    return experiments

experiments = discover_experiments()
if not experiments:
    st.error("No experiments found. Run `python gradient_descent.py run <plate>` first.")
    st.stop()

# Sidebar
with st.sidebar:
    selected_experiment = st.selectbox("Experiment", experiments, index=len(experiments) - 1)
    auto_refresh = st.toggle("Auto-refresh (5s)", value=True)
    st.divider()
    st.caption("Gradient Descent Media Optimization")
    st.caption(f"Data: `{BASE_DATA_DIR / selected_experiment}`")

DATA_DIR = BASE_DATA_DIR / selected_experiment

# Load data
state = load_state()
current_iter = state.get("current_iteration", 0)
history = state.get("history", [])
is_delta = uses_delta_od(history)

# Enrich history with od_results from iteration logs if missing
for h in history:
    if "od_results" not in h or not h["od_results"]:
        log = load_iteration_log(h["iteration"])
        if log and "od_results" in log:
            h["od_results"] = log["od_results"]

# Sidebar info (after data load so we can show objective type)
with st.sidebar:
    st.divider()
    obj_label = "Delta OD (endpoint - baseline)" if is_delta else "Absolute OD600"
    st.markdown(f"**Objective:** {obj_label}")
    st.divider()
    st.markdown("**How it works**")
    if is_delta:
        st.markdown(
            "Each iteration tests the current media recipe (center) plus "
            "small perturbations of each supplement. Growth is measured as "
            "**delta OD** (endpoint - baseline) to normalize for initial cell density. "
            "The gradient tells us which supplements improve growth. "
            "The optimizer steps toward higher growth, converging on the best recipe.",
            help="Gradient descent on a 96-well plate"
        )
    else:
        st.markdown(
            "Each iteration tests the current media recipe (center) plus "
            "small perturbations of each supplement. The gradient of OD600 "
            "tells us which direction to adjust. The optimizer steps toward "
            "higher growth, converging on the best recipe.",
            help="Gradient descent on a 96-well plate"
        )
    st.markdown("**Well layout per row**")
    st.markdown(
        "- **Col 1** -- Seed well (cells stock)\n"
        "- **Col 2** -- Neg control (no cells)\n"
        "- **Col 3** -- Pos control (cells only)\n"
        "- **Col 4** -- Center (current best)\n"
        f"- **Cols 5-6** -- +{SUPPLEMENT_NAMES[0]} perturbation\n"
        f"- **Cols 7-8** -- +{SUPPLEMENT_NAMES[1]} perturbation\n"
        f"- **Cols 9-10** -- +{SUPPLEMENT_NAMES[2]} perturbation\n"
        "- **Cols 11-12** -- Extra wells"
    )

# Detect status
workflow_ids = load_workflow_ids(current_iter + 1) if current_iter < MAX_ITERATIONS else None
live_status = None
if workflow_ids:
    live_status = get_live_workflow_status(workflow_ids)

is_running = live_status and live_status.get("status") in ("in_progress", "pending_approval", "approved")
is_converged = state.get("converged", False)

# Header
col_title, col_status = st.columns([3, 1])
with col_title:
    st.title("GD Media Optimization")
with col_status:
    if is_running:
        running_iter = current_iter + 1
        st.markdown(f'<div class="status-pill status-running">Running Iteration {running_iter}</div>', unsafe_allow_html=True)
    elif is_converged:
        st.markdown('<div class="status-pill status-converged">Converged</div>', unsafe_allow_html=True)
    else:
        status_text = f"Idle -- {current_iter}/{MAX_ITERATIONS} done" if current_iter > 0 else "Not started"
        st.markdown(f'<div class="status-pill status-idle">{status_text}</div>', unsafe_allow_html=True)

# Metrics row
m1, m2, m3, m4 = st.columns(4)
with m1:
    st.metric("Iteration", f"{current_iter} / {MAX_ITERATIONS}")
with m2:
    best_od = state.get("best_od")
    if is_delta:
        st.metric("Best Growth (dOD)", f"{best_od:+.4f}" if best_od is not None else "--")
    else:
        st.metric("Best Center OD", f"{best_od:.4f}" if best_od else "--")
with m3:
    st.metric("Learning Rate (a)", f"{state.get('alpha', 1.0):.2f}")
with m4:
    streak = state.get("no_improvement_count", 0)
    st.metric("No-improvement streak", f"{streak} / 2")

st.divider()

# Objective function chart — the hero chart
if is_delta:
    st.subheader("Growth Optimization")
    st.markdown(
        '<div class="explanation">'
        "The objective function: center-point growth (delta OD = endpoint - baseline) at each iteration. "
        "Delta OD normalizes for initial cell density, so we're measuring actual growth from the recipe. "
        "If the optimization is working, this line should trend upward. "
        "The dotted green line tracks the best growth seen so far."
        '</div>',
        unsafe_allow_html=True,
    )
else:
    st.subheader("Optimization Score")
    st.markdown(
        '<div class="explanation">'
        "The objective function: center-point OD600 at each trial. "
        "If the optimization is working, this line should trend upward and flatten as it converges. "
        "The dotted green line tracks the best score seen so far."
        '</div>',
        unsafe_allow_html=True,
    )
st.plotly_chart(objective_chart(history, is_delta), use_container_width=True, key="objective")

st.divider()

# Row 1: Plate + Composition
col_plate, col_comp = st.columns([3, 2])
with col_plate:
    st.subheader("Plate View")
    if is_delta:
        st.markdown(
            '<div class="explanation">'
            "96-well plate colored by growth (delta OD). Green = positive growth, red = negative growth, "
            "white = no change. Column 1 holds seed wells. Each row (A-H) is one iteration. "
            "Hover over any well to see growth value and reagent volumes."
            '</div>',
            unsafe_allow_html=True,
        )
    else:
        st.markdown(
            '<div class="explanation">'
            "96-well plate colored by OD600 (darker green = higher growth). "
            "Column 1 holds seed wells. Each row (A-H) is one iteration. "
            "Hover over any well to see exact reagent volumes."
            '</div>',
            unsafe_allow_html=True,
        )
    st.plotly_chart(plate_heatmap(state, is_delta), use_container_width=True, key="plate")

with col_comp:
    st.subheader("Current Composition")
    st.markdown(
        '<div class="explanation">'
        "The media recipe the optimizer will test next. Each well is 180 uL total -- "
        f"Novel_Bio is the base media, and the three supplements ({', '.join(SUPPLEMENT_NAMES)}) "
        "are the variables being optimized. As supplements increase, Novel_Bio decreases."
        '</div>',
        unsafe_allow_html=True,
    )
    comp = state.get("current_composition", {name: 20 for name in SUPPLEMENT_NAMES})
    st.plotly_chart(composition_bar(comp), use_container_width=True, key="comp_bar")

# Row 2: Progress charts
col_od, col_traj = st.columns(2)
with col_od:
    if is_delta:
        st.subheader("Growth by Well Role")
        st.markdown(
            '<div class="explanation">'
            "Growth (delta OD) across iterations for each well role. The <b>center</b> (blue) is the "
            "current recipe — we want this to go up. The <b>control</b> (gray dashed) is Novel_Bio only. "
            "Colored lines show mean growth when each supplement is added. If a colored line is above "
            "center, that supplement improves growth. The dotted line at zero separates growth from decline."
            '</div>',
            unsafe_allow_html=True,
        )
    else:
        st.subheader("OD600 Progress")
        st.markdown(
            '<div class="explanation">'
            "Tracks bacterial growth (OD600) across iterations. The <b>center</b> line (blue) "
            "is the current best recipe — we want this to go up. The <b>control</b> (gray dashed) "
            "is Novel_Bio only, serving as a baseline. Colored lines show the mean growth when each "
            "supplement was increased. If a colored line is above center, that supplement helps growth."
            '</div>',
            unsafe_allow_html=True,
        )
    st.plotly_chart(od_progress_chart(history, is_delta), use_container_width=True, key="od_progress")
with col_traj:
    st.subheader("Composition Trajectory")
    st.markdown(
        '<div class="explanation">'
        "How the media recipe changes over time. Each stacked bar shows the volume breakdown "
        "for that iteration's center point. The optimizer increases supplements with positive "
        "gradients and decreases those with negative gradients, always keeping the total at 180 uL."
        '</div>',
        unsafe_allow_html=True,
    )
    st.plotly_chart(composition_trajectory_chart(history), use_container_width=True, key="comp_traj")

# Row 3: Current iteration detail
if is_running and live_status:
    st.divider()
    st.subheader(f"Iteration {current_iter + 1} -- Live")
    routines = live_status.get("workflow_routines", [])
    completed = sum(1 for r in routines if r.get("status") == "completed")
    total = len(routines)
    st.progress(completed / max(total, 1), text=f"{completed}/{total} routines completed")

    routine_df = pd.DataFrame(
        [{"Routine": r.get("routine_name", "?"), "Status": r.get("status", "?")} for r in routines]
    )
    st.dataframe(routine_df, use_container_width=True, hide_index=True, height=200)

elif history:
    st.divider()
    latest = history[-1]
    latest_log = load_iteration_log(latest["iteration"]) or latest
    gradient = latest.get("gradient", latest_log.get("gradient"))

    if gradient:
        st.subheader(f"Iteration {latest['iteration']} -- Gradient Direction")
        st.markdown(
            '<div class="explanation">'
            "Direction the optimizer will move next. An up arrow means increasing that supplement "
            "improved growth; a down arrow means it hurt growth. "
            "The perturbation mean shows the average growth for wells with extra supplement vs the center."
            '</div>',
            unsafe_allow_html=True,
        )
        g1, g2, g3 = st.columns(3)

        # Get perturbation means for display
        od_results = latest.get("od_results", {})
        if not od_results:
            log = load_iteration_log(latest["iteration"])
            od_results = log.get("od_results", {}) if log else {}
        perturbed_ods = od_results.get("perturbed_ods", {})
        center_val = latest.get("center_od", 0)

        for col, name in zip([g1, g2, g3], SUPPLEMENT_NAMES):
            val = gradient.get(name, 0)
            arrow = "^" if val > 0 else ("v" if val < 0 else "-")
            color = "#15803d" if val > 0 else ("#ef4444" if val < 0 else "#94a3b8")

            # Perturbation mean
            reps = perturbed_ods.get(name, [])
            pert_mean = sum(reps) / max(len(reps), 1) if reps else 0
            diff = pert_mean - center_val

            with col:
                if is_delta:
                    detail = f"Mean: {pert_mean:+.4f} ({diff:+.4f} vs center)"
                else:
                    detail = f"Mean: {pert_mean:.4f} ({diff:+.4f} vs center)"
                st.markdown(
                    f"<div class='metric-card'>"
                    f"<div class='metric-value' style='color:{color}'>{arrow}</div>"
                    f"<div class='metric-label'>{name} (grad: {val:+.4f})</div>"
                    f"<div style='font-size:11px;color:#94a3b8;margin-top:4px'>{detail}</div>"
                    f"</div>",
                    unsafe_allow_html=True,
                )

# History table
if history:
    st.divider()
    st.subheader("Iteration History")
    rows = []
    for h in history:
        c = h.get("composition", {})
        row = {
            "Iter": h["iteration"],
            "Novel_Bio": c.get("Novel_Bio", "--"),
            **{name: c.get(name, "--") for name in SUPPLEMENT_NAMES},
        }
        if is_delta:
            row["Growth (dOD)"] = f"{h['center_od']:+.4f}" if "center_od" in h else "--"
            row["Control Growth"] = f"{h['control_od']:+.4f}" if "control_od" in h else "--"
            # Show absolute OD if available
            if "abs_center_od" in h:
                row["Abs. Center OD"] = f"{h['abs_center_od']:.4f}"
                row["Abs. Control OD"] = f"{h['abs_control_od']:.4f}"
        else:
            row["Center OD"] = f"{h['center_od']:.4f}" if "center_od" in h else "--"
            row["Control OD"] = f"{h['control_od']:.4f}" if "control_od" in h else "--"
        row["a"] = h.get("alpha", "--")
        rows.append(row)
    st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)

# Auto-refresh
if auto_refresh:
    time.sleep(5)
    st.rerun()
else:
    st.divider()
    if st.button("Refresh data"):
        st.rerun()
