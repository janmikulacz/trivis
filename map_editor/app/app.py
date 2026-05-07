"""TřiVis Polygonal Map Editor — Dash entry point.

Build the trivis_clipper extension first:
    cd map_editor
    cmake -B build -DCMAKE_BUILD_TYPE=Release
    cmake --build build --parallel

Then run:
    cd map_editor/app
    python app.py
"""

import dash

app = dash.Dash(
    __name__,
    suppress_callback_exceptions=True,
    title="TřiVis Map Editor",
)

from layout import layout      # noqa: E402 — import after app creation
import callbacks               # noqa: E402 — registers callbacks via decorators

app.layout = layout

if __name__ == "__main__":
    app.run(debug=True, port=8050, host="127.0.0.1")
