"""Dash layout for the TřiVis Polygonal Map Editor."""

from dash import dcc, html
import plotly.graph_objects as go

CANVAS_RANGE = [0.0, 100.0]

# ── Canvas helpers ────────────────────────────────────────────────────────────

def blank_figure() -> go.Figure:
    """Return a blank canvas figure with an invisible background for click events."""
    r0, r1 = CANVAS_RANGE
    fig = go.Figure()

    # Transparent filled rectangle — clicking anywhere in the canvas fires clickData
    fig.add_trace(go.Scatter(
        x=[r0, r1, r1, r0, r0],
        y=[r0, r0, r1, r1, r0],
        mode="lines",
        fill="toself",
        fillcolor="rgba(255,255,255,0)",
        line=dict(color="rgba(0,0,0,0)", width=0),
        hovertemplate="x: %{x:.2f}  y: %{y:.2f}<extra></extra>",
        showlegend=False,
        name="__bg__",
    ))

    fig.update_layout(
        xaxis=dict(range=CANVAS_RANGE, showgrid=True, gridcolor="#e8e8e8",
                   zeroline=False, title="x"),
        yaxis=dict(range=CANVAS_RANGE, scaleanchor="x", showgrid=True,
                   gridcolor="#e8e8e8", zeroline=False, title="y"),
        margin=dict(l=40, r=10, t=10, b=40),
        plot_bgcolor="white",
        paper_bgcolor="white",
        dragmode="pan",
        hovermode="closest",
        showlegend=False,
        uirevision="keep",  # preserve zoom/pan across figure updates
    )
    return fig


# ── Button style helpers ──────────────────────────────────────────────────────

def _btn(label: str, btn_id: str, color: str = "#4a90d9") -> html.Button:
    return html.Button(label, id=btn_id, n_clicks=0, style={
        "display": "block",
        "width": "100%",
        "margin": "3px 0",
        "padding": "6px 8px",
        "background": color,
        "color": "white",
        "border": "none",
        "borderRadius": "4px",
        "cursor": "pointer",
        "fontSize": "13px",
        "textAlign": "left",
    })


def _section(title: str, *children) -> html.Div:
    return html.Div([
        html.P(title, style={
            "fontWeight": "bold",
            "fontSize": "12px",
            "color": "#555",
            "margin": "10px 0 4px",
            "textTransform": "uppercase",
            "letterSpacing": "0.5px",
        }),
        *children,
    ])


# ── Full layout ───────────────────────────────────────────────────────────────

layout = html.Div([

    # ── Hidden stores ─────────────────────────────────────────────────────────
    dcc.Store(id="state", data={
        "polygons": [],   # [{id, vertices:[[x,y],...], role:"border"|"hole", selected:bool}]
        "in_progress": [],
        "mode": "draw",
        "next_id": 0,
    }),
    dcc.Store(id="click-coords", data=None),  # set by clientside callback
    dcc.Download(id="download-map"),

    # ── Page header ───────────────────────────────────────────────────────────
    html.Div([
        html.H2("TřiVis Map Editor", style={
            "margin": "0",
            "fontSize": "18px",
            "fontWeight": "600",
        }),
        html.Span("Clipper2-powered polygonal map editor",
                  style={"fontSize": "12px", "color": "#888"}),
    ], style={
        "padding": "10px 16px",
        "background": "#1e2a3a",
        "color": "white",
        "display": "flex",
        "alignItems": "baseline",
        "gap": "12px",
    }),

    # ── Main area ─────────────────────────────────────────────────────────────
    html.Div([

        # ── Left toolbar ──────────────────────────────────────────────────────
        html.Div([

            _section("Mode",
                dcc.RadioItems(
                    id="mode-radio",
                    options=[
                        {"label": " Draw polygon", "value": "draw"},
                        {"label": " Select / edit", "value": "select"},
                    ],
                    value="draw",
                    labelStyle={"display": "block", "fontSize": "13px",
                                "margin": "3px 0", "cursor": "pointer"},
                ),
            ),

            _section("Drawing",
                html.Div([
                    html.Label("x", style={"fontSize": "12px", "marginRight": "4px"}),
                    dcc.Input(id="input-x", type="number", placeholder="0.0",
                              debounce=False, style={"width": "60px", "fontSize": "12px"}),
                    html.Label("y", style={"fontSize": "12px", "margin": "0 4px"}),
                    dcc.Input(id="input-y", type="number", placeholder="0.0",
                              debounce=False, style={"width": "60px", "fontSize": "12px"}),
                ], style={"display": "flex", "alignItems": "center", "marginBottom": "4px"}),
                html.Small("(hover canvas to auto-fill)", style={"color": "#999", "fontSize": "11px"}),
                _btn("➕  Add vertex",     "btn-add-vertex",   "#2e7d32"),
                _btn("✓  Close polygon",  "btn-close-poly",   "#1565c0"),
                _btn("↩  Undo vertex",    "btn-undo-vertex",  "#e65100"),
                _btn("✗  Clear drawing",  "btn-clear-draw",   "#b71c1c"),
            ),

            _section("Boolean ops (select 2)",
                _btn("∪  Union",          "btn-union",        "#6a1b9a"),
                _btn("∩  Intersection",   "btn-intersect",    "#4527a0"),
                _btn("−  Difference A∖B", "btn-difference",   "#283593"),
                _btn("△  XOR",            "btn-xor",          "#01579b"),
            ),

            _section("Offset selected",
                html.Div([
                    dcc.Input(id="input-offset", type="number", value=2.0,
                              placeholder="delta", debounce=False,
                              style={"width": "70px", "fontSize": "12px"}),
                    html.Button("Apply", id="btn-offset", n_clicks=0, style={
                        "marginLeft": "6px", "padding": "4px 8px",
                        "fontSize": "12px", "cursor": "pointer",
                    }),
                ], style={"display": "flex", "alignItems": "center"}),
            ),

            _section("Simplify selected",
                html.Div([
                    html.Label("ε", style={"fontSize": "12px", "marginRight": "4px"}),
                    dcc.Input(id="input-epsilon", type="number", value=0.5,
                              placeholder="epsilon", debounce=False,
                              style={"width": "70px", "fontSize": "12px"}),
                    html.Button("Apply", id="btn-simplify", n_clicks=0, style={
                        "marginLeft": "6px", "padding": "4px 8px",
                        "fontSize": "12px", "cursor": "pointer",
                    }),
                ], style={"display": "flex", "alignItems": "center"}),
            ),

            _section("Selected polygon role",
                _btn("🔵  Set as border",  "btn-set-border",   "#00695c"),
                _btn("🔴  Set as hole",    "btn-set-hole",     "#c62828"),
                _btn("🗑  Delete selected","btn-delete",       "#37474f"),
            ),

            _section("File",
                _btn("📂  Import .txt map", "btn-import-open", "#37474f"),
                dcc.Upload(id="upload-map", children=html.Div(), style={"display": "none"}),
                _btn("💾  Export .txt map", "btn-export",       "#37474f"),
            ),

        ], style={
            "width": "220px",
            "minWidth": "220px",
            "padding": "8px 12px",
            "overflowY": "auto",
            "background": "#f9f9f9",
            "borderRight": "1px solid #ddd",
        }),

        # ── Canvas + polygon list ──────────────────────────────────────────────
        html.Div([

            # Canvas
            dcc.Graph(
                id="canvas",
                figure=blank_figure(),
                style={"flex": "1", "height": "100%"},
                config={"scrollZoom": True, "displaylogo": False,
                        "modeBarButtonsToRemove": ["toImage"]},
            ),

            # Polygon list (right of canvas)
            html.Div([
                html.P("Polygons", style={
                    "fontWeight": "bold", "fontSize": "12px",
                    "color": "#555", "margin": "0 0 6px",
                }),
                html.Div(id="polygon-list",
                         style={"overflowY": "auto", "flex": "1"}),
            ], style={
                "width": "160px",
                "padding": "8px",
                "background": "#f9f9f9",
                "borderLeft": "1px solid #ddd",
                "display": "flex",
                "flexDirection": "column",
            }),

        ], style={"display": "flex", "flex": "1", "overflow": "hidden"}),

    ], style={
        "display": "flex",
        "flex": "1",
        "overflow": "hidden",
    }),

    # ── Status bar ────────────────────────────────────────────────────────────
    html.Div(id="status-bar", children="Ready. Click canvas to draw, or import a map.",
             style={
                 "padding": "4px 12px",
                 "background": "#e8eaf6",
                 "fontSize": "12px",
                 "color": "#333",
                 "borderTop": "1px solid #ddd",
             }),

], style={
    "display": "flex",
    "flexDirection": "column",
    "height": "100vh",
    "fontFamily": "sans-serif",
})
