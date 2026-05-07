"""Dash callbacks for the TřiVis Map Editor."""

from __future__ import annotations

import base64
import copy
import io
import json
import tempfile

import plotly.graph_objects as go
from dash import Input, Output, State, callback, clientside_callback, ctx, no_update
from dash.exceptions import PreventUpdate

import clipper_utils as cu
from layout import CANVAS_RANGE, blank_figure

# ── Colour palette ────────────────────────────────────────────────────────────

_BORDER_FILL = "rgba(33, 150, 243, 0.15)"
_BORDER_LINE = "#1565c0"
_HOLE_FILL   = "rgba(229, 57, 53, 0.15)"
_HOLE_LINE   = "#b71c1c"
_SEL_LINE    = "#fbc02d"
_SEL_WIDTH   = 3
_NORM_WIDTH  = 2
_PROG_FILL   = "rgba(0, 200, 83, 0.10)"
_PROG_LINE   = "#00c853"

# ── Clientside callback: capture actual mouse position before plotly snaps ────

clientside_callback(
    """
    function(clickData) {
        const coords = window._lastClickCoords;
        if (!coords) return window.dash_clientside.no_update;
        return coords;
    }
    """,
    Output("click-coords", "data"),
    Input("canvas", "clickData"),
)

# ── Figure rendering ──────────────────────────────────────────────────────────

def _render_figure(state: dict) -> go.Figure:
    fig = blank_figure()

    for poly in state["polygons"]:
        verts = poly["vertices"]
        if not verts:
            continue
        closed = verts + [verts[0]]
        xs = [v[0] for v in closed]
        ys = [v[1] for v in closed]

        is_border  = poly["role"] == "border"
        selected   = poly["selected"]
        fill_color = _BORDER_FILL if is_border else _HOLE_FILL
        line_color = _SEL_LINE if selected else (_BORDER_LINE if is_border else _HOLE_LINE)
        line_width = _SEL_WIDTH if selected else _NORM_WIDTH
        role_label = "Border" if is_border else "Hole"
        name       = f"{role_label} #{poly['id']}"

        fig.add_trace(go.Scatter(
            x=xs, y=ys,
            mode="lines+markers",
            fill="toself",
            fillcolor=fill_color,
            line=dict(color=line_color, width=line_width),
            marker=dict(size=6, color=line_color),
            name=name,
            customdata=[[poly["id"]]] * len(xs),
            hovertemplate=f"{name}<br>x: %{{x:.2f}}  y: %{{y:.2f}}<extra></extra>",
        ))

    # In-progress polygon
    prog = state["in_progress"]
    if prog:
        xs = [v[0] for v in prog]
        ys = [v[1] for v in prog]
        fig.add_trace(go.Scatter(
            x=xs, y=ys,
            mode="lines+markers",
            fill="toself" if len(prog) >= 3 else "none",
            fillcolor=_PROG_FILL,
            line=dict(color=_PROG_LINE, width=2, dash="dash"),
            marker=dict(size=7, color=_PROG_LINE),
            name="In progress",
            hoverinfo="skip",
        ))

    fig.update_layout(dragmode="pan" if state["mode"] == "select" else "pan")
    return fig


# ── Polygon list panel ────────────────────────────────────────────────────────

def _render_poly_list(state: dict):
    from dash import html
    items = []
    for poly in state["polygons"]:
        n = len(poly["vertices"])
        role = poly["role"]
        sel  = poly["selected"]
        color = "#1565c0" if role == "border" else "#b71c1c"
        items.append(html.Div(
            f"{'★ ' if sel else ''}{'Border' if role == 'border' else 'Hole'} #{poly['id']} ({n}v)",
            id={"type": "poly-item", "index": poly["id"]},
            n_clicks=0,
            style={
                "padding": "4px 6px",
                "marginBottom": "2px",
                "background": "#fff3e0" if sel else "white",
                "border": f"1px solid {color}",
                "borderRadius": "3px",
                "fontSize": "12px",
                "cursor": "pointer",
                "color": color,
            }
        ))
    return items or [html.Span("No polygons yet.", style={"fontSize": "12px", "color": "#aaa"})]


# ── Hover → auto-fill coordinate inputs ──────────────────────────────────────

@callback(
    Output("input-x", "value"),
    Output("input-y", "value"),
    Input("canvas", "hoverData"),
    prevent_initial_call=True,
)
def update_coords_from_hover(hover_data):
    if not hover_data or not hover_data.get("points"):
        raise PreventUpdate
    pt = hover_data["points"][0]
    x = round(pt.get("x", 0), 2)
    y = round(pt.get("y", 0), 2)
    return x, y


# ── Canvas click → select polygon (select mode) ───────────────────────────────

@callback(
    Output("state", "data", allow_duplicate=True),
    Output("status-bar", "children", allow_duplicate=True),
    Input("canvas", "clickData"),
    State("state", "data"),
    State("mode-radio", "value"),
    prevent_initial_call=True,
)
def handle_canvas_click(click_data, state, mode):
    if mode != "select" or not click_data or not click_data.get("points"):
        raise PreventUpdate

    pt = click_data["points"][0]
    custom = pt.get("customdata")
    if custom is None:
        raise PreventUpdate

    poly_id = custom[0]
    state = copy.deepcopy(state)
    for poly in state["polygons"]:
        if poly["id"] == poly_id:
            poly["selected"] = not poly["selected"]
            verb = "Selected" if poly["selected"] else "Deselected"
            return state, f"{verb} {poly['role']} #{poly_id}."
    raise PreventUpdate


# ── Main state-mutation callback ──────────────────────────────────────────────

@callback(
    Output("state", "data"),
    Output("canvas", "figure"),
    Output("polygon-list", "children"),
    Output("status-bar", "children"),
    Input("btn-add-vertex",  "n_clicks"),
    Input("btn-close-poly",  "n_clicks"),
    Input("btn-undo-vertex", "n_clicks"),
    Input("btn-clear-draw",  "n_clicks"),
    Input("btn-union",       "n_clicks"),
    Input("btn-intersect",   "n_clicks"),
    Input("btn-difference",  "n_clicks"),
    Input("btn-xor",         "n_clicks"),
    Input("btn-offset",      "n_clicks"),
    Input("btn-simplify",    "n_clicks"),
    Input("btn-set-border",  "n_clicks"),
    Input("btn-set-hole",    "n_clicks"),
    Input("btn-delete",      "n_clicks"),
    Input("btn-export",      "n_clicks"),
    Input("upload-map",      "contents"),
    State("state",           "data"),
    State("click-coords",    "data"),
    State("input-x",         "value"),
    State("input-y",         "value"),
    State("input-offset",    "value"),
    State("input-epsilon",   "value"),
    prevent_initial_call=True,
)
def handle_actions(
    _add, _close, _undo, _clear,
    _union, _intersect, _diff, _xor,
    _offset, _simplify,
    _border, _hole, _delete,
    _export,
    upload_contents,
    state, click_coords,
    input_x, input_y,
    offset_delta, epsilon,
):
    triggered = ctx.triggered_id
    state = copy.deepcopy(state)
    status = no_update

    # ── Add vertex ────────────────────────────────────────────────────────────
    if triggered == "btn-add-vertex":
        # Prefer JS-captured coordinates; fall back to typed inputs
        if click_coords:
            x, y = click_coords["x"], click_coords["y"]
        elif input_x is not None and input_y is not None:
            x, y = float(input_x), float(input_y)
        else:
            status = "Enter x and y coordinates or hover over the canvas."
            return state, no_update, no_update, status
        state["in_progress"].append([x, y])
        n = len(state["in_progress"])
        status = f"Added vertex {n}: ({x}, {y})"

    # ── Close polygon ─────────────────────────────────────────────────────────
    elif triggered == "btn-close-poly":
        verts = state["in_progress"]
        if len(verts) < 3:
            status = "Need at least 3 vertices to close a polygon."
            return state, no_update, no_update, status
        pid = state["next_id"]
        state["polygons"].append({
            "id": pid,
            "vertices": verts,
            "role": "border",
            "selected": False,
        })
        state["in_progress"] = []
        state["next_id"] = pid + 1
        status = f"Polygon #{pid} closed ({len(verts)} vertices). Role: border."

    # ── Undo vertex ───────────────────────────────────────────────────────────
    elif triggered == "btn-undo-vertex":
        if state["in_progress"]:
            state["in_progress"].pop()
            status = f"Removed last vertex ({len(state['in_progress'])} remaining)."
        else:
            status = "Nothing to undo."

    # ── Clear drawing ─────────────────────────────────────────────────────────
    elif triggered == "btn-clear-draw":
        state["in_progress"] = []
        status = "Drawing cleared."

    # ── Boolean ops ───────────────────────────────────────────────────────────
    elif triggered in ("btn-union", "btn-intersect", "btn-difference", "btn-xor"):
        selected = [p for p in state["polygons"] if p["selected"]]
        if len(selected) < 2:
            status = "Select at least 2 polygons first."
            return state, no_update, no_update, status
        subjects = [selected[0]["vertices"]]
        clips    = [selected[1]["vertices"]]
        op_name  = triggered.removeprefix("btn-")
        try:
            if triggered == "btn-union":
                result_verts = cu.union(subjects, clips)
                op_name = "Union"
            elif triggered == "btn-intersect":
                result_verts = cu.intersect(subjects, clips)
                op_name = "Intersection"
            elif triggered == "btn-difference":
                result_verts = cu.difference(subjects, clips)
                op_name = "Difference"
            else:
                result_verts = cu.xor(subjects, clips)
                op_name = "XOR"
        except Exception as exc:
            return state, no_update, no_update, f"Operation failed: {exc}"

        # Remove the two selected polygons and add results
        sel_ids = {p["id"] for p in selected}
        state["polygons"] = [p for p in state["polygons"] if p["id"] not in sel_ids]
        for verts in result_verts:
            pid = state["next_id"]
            state["polygons"].append({
                "id": pid, "vertices": verts,
                "role": "border", "selected": False,
            })
            state["next_id"] = pid + 1
        status = f"{op_name}: produced {len(result_verts)} polygon(s)."

    # ── Offset ────────────────────────────────────────────────────────────────
    elif triggered == "btn-offset":
        selected = [p for p in state["polygons"] if p["selected"]]
        if not selected:
            status = "Select a polygon first."
            return state, no_update, no_update, status
        delta = float(offset_delta or 2.0)
        for poly in selected:
            try:
                results = cu.inflate([poly["vertices"]], delta)
                if results:
                    poly["vertices"] = results[0]
            except Exception as exc:
                status = f"Offset failed: {exc}"
                return state, no_update, no_update, status
        status = f"Offset by {delta} applied to {len(selected)} polygon(s)."

    # ── Simplify ──────────────────────────────────────────────────────────────
    elif triggered == "btn-simplify":
        selected = [p for p in state["polygons"] if p["selected"]]
        if not selected:
            status = "Select a polygon first."
            return state, no_update, no_update, status
        eps = float(epsilon or 0.5)
        for poly in selected:
            try:
                results = cu.simplify([poly["vertices"]], eps)
                if results:
                    poly["vertices"] = results[0]
            except Exception as exc:
                status = f"Simplify failed: {exc}"
                return state, no_update, no_update, status
        status = f"Simplify (ε={eps}) applied to {len(selected)} polygon(s)."

    # ── Set role ──────────────────────────────────────────────────────────────
    elif triggered in ("btn-set-border", "btn-set-hole"):
        new_role = "border" if triggered == "btn-set-border" else "hole"
        count = 0
        for poly in state["polygons"]:
            if poly["selected"]:
                poly["role"] = new_role
                count += 1
        status = f"Set {count} polygon(s) to '{new_role}'."

    # ── Delete ────────────────────────────────────────────────────────────────
    elif triggered == "btn-delete":
        before = len(state["polygons"])
        state["polygons"] = [p for p in state["polygons"] if not p["selected"]]
        removed = before - len(state["polygons"])
        status = f"Deleted {removed} polygon(s)."

    # ── Export ────────────────────────────────────────────────────────────────
    elif triggered == "btn-export":
        borders = [p["vertices"] for p in state["polygons"] if p["role"] == "border"]
        holes   = [p["vertices"] for p in state["polygons"] if p["role"] == "hole"]
        if not borders:
            status = "No border polygon defined. Mark at least one polygon as border."
            return state, no_update, no_update, status
        buf = io.StringIO()
        # Write to string buffer using the same format as save_trivis_map
        buf.write("[SCALE]\n1.0\n")
        for verts in borders:
            buf.write("\n[BORDER]\n")
            for x, y in verts:
                buf.write(f"{x:.17g} {y:.17g}\n")
        for verts in holes:
            buf.write("\n[OBSTACLE]\n")
            for x, y in verts:
                buf.write(f"{x:.17g} {y:.17g}\n")
        content = buf.getvalue()
        # Trigger download via dcc.Download — handled separately
        status = "Map exported. Check your downloads."
        # Encode and return via the download component
        from dash import dcc as _dcc
        return (
            state,
            _render_figure(state),
            _render_poly_list(state),
            status,
        )

    # ── Import ────────────────────────────────────────────────────────────────
    elif triggered == "upload-map" and upload_contents:
        try:
            _, content_str = upload_contents.split(",", 1)
            decoded = base64.b64decode(content_str).decode("utf-8")
            border, holes = cu.load_trivis_map(io.StringIO(decoded))  # type: ignore[arg-type]
            # Actually load_trivis_map expects a file path; use temp file
            with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as tmp:
                tmp.write(decoded)
                tmp_path = tmp.name
            border, holes = cu.load_trivis_map(tmp_path)
            state["polygons"] = []
            state["next_id"] = 0
            if border:
                state["polygons"].append({
                    "id": state["next_id"], "vertices": border,
                    "role": "border", "selected": False,
                })
                state["next_id"] += 1
            for hole in holes:
                state["polygons"].append({
                    "id": state["next_id"], "vertices": hole,
                    "role": "hole", "selected": False,
                })
                state["next_id"] += 1
            status = f"Loaded {1 if border else 0} border + {len(holes)} hole(s)."
        except Exception as exc:
            status = f"Import failed: {exc}"

    return state, _render_figure(state), _render_poly_list(state), status


# ── Mode radio → update state ─────────────────────────────────────────────────

@callback(
    Output("state", "data", allow_duplicate=True),
    Input("mode-radio", "value"),
    State("state", "data"),
    prevent_initial_call=True,
)
def sync_mode(mode, state):
    state = copy.deepcopy(state)
    state["mode"] = mode
    return state


# ── Trigger file-upload dialog via hidden dcc.Upload ─────────────────────────

clientside_callback(
    """
    function(n) {
        if (!n) return window.dash_clientside.no_update;
        const el = document.querySelector('#upload-map input[type=file]');
        if (el) el.click();
        return window.dash_clientside.no_update;
    }
    """,
    Output("upload-map", "contents"),
    Input("btn-import-open", "n_clicks"),
    prevent_initial_call=True,
)


# ── Export download ───────────────────────────────────────────────────────────

@callback(
    Output("download-map", "data"),
    Input("btn-export", "n_clicks"),
    State("state", "data"),
    prevent_initial_call=True,
)
def export_map(_, state):
    borders = [p["vertices"] for p in state["polygons"] if p["role"] == "border"]
    holes   = [p["vertices"] for p in state["polygons"] if p["role"] == "hole"]
    if not borders:
        raise PreventUpdate

    buf = io.StringIO()
    buf.write("[SCALE]\n1.0\n")
    for verts in borders:
        buf.write("\n[BORDER]\n")
        for x, y in verts:
            buf.write(f"{x:.17g} {y:.17g}\n")
    for verts in holes:
        buf.write("\n[OBSTACLE]\n")
        for x, y in verts:
            buf.write(f"{x:.17g} {y:.17g}\n")

    return {"content": buf.getvalue(), "filename": "map.txt"}
