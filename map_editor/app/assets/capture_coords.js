/**
 * Captures the actual mouse position in data coordinates whenever the user
 * clicks on the canvas.  Plotly's clickData snaps to the nearest data point,
 * so we intercept the raw mousedown event and convert using plotly's internal
 * axis geometry (range + offset + length).
 *
 * The result is stored in window._lastClickCoords = {x, y} and is read by the
 * Dash clientside callback in callbacks.py.
 */
(function () {
    "use strict";

    function svgPixelToData(svgX, svgY, layout) {
        const xa = layout.xaxis;
        const ya = layout.yaxis;
        if (!xa || !ya || !xa._length || !ya._length) return null;

        const dataX = xa.range[0] + (svgX - xa._offset) / xa._length * (xa.range[1] - xa.range[0]);
        // Y axis is inverted in SVG (pixel 0 = top = ya.range[1] in data coords)
        const dataY = ya.range[1] - (svgY - ya._offset) / ya._length * (ya.range[1] - ya.range[0]);

        // Clamp to current axis range
        const xMin = Math.min(xa.range[0], xa.range[1]);
        const xMax = Math.max(xa.range[0], xa.range[1]);
        const yMin = Math.min(ya.range[0], ya.range[1]);
        const yMax = Math.max(ya.range[0], ya.range[1]);

        return {
            x: Math.round(Math.min(Math.max(dataX, xMin), xMax) * 100) / 100,
            y: Math.round(Math.min(Math.max(dataY, yMin), yMax) * 100) / 100
        };
    }

    function attachListener() {
        const graphDiv = document.getElementById("canvas");
        if (!graphDiv) {
            // Retry until Dash has rendered the graph
            setTimeout(attachListener, 150);
            return;
        }

        graphDiv.addEventListener("mousedown", function (event) {
            if (!graphDiv._fullLayout) return;
            const svg = graphDiv.querySelector(".main-svg");
            if (!svg) return;
            const rect = svg.getBoundingClientRect();
            const coords = svgPixelToData(
                event.clientX - rect.left,
                event.clientY - rect.top,
                graphDiv._fullLayout
            );
            if (coords) {
                window._lastClickCoords = coords;
            }
        });
    }

    if (document.readyState === "loading") {
        document.addEventListener("DOMContentLoaded", attachListener);
    } else {
        setTimeout(attachListener, 0);
    }
})();
