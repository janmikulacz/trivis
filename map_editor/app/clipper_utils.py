"""
Friendly Python wrapper around the compiled trivis_clipper C++ extension.

All public functions accept and return polygons as ``list[list[[x, y]]]`` —
i.e. a list of paths, each path being a list of [x, y] coordinate pairs.
The trivis_clipper module's PointD / PathsD types are kept internal.
"""

from __future__ import annotations

import os
import sys
from typing import TYPE_CHECKING

import trivis_clipper as _tc
from trivis_clipper import FillRule, JoinType, EndType  # re-export enums

# ── Type helpers ──────────────────────────────────────────────────────────────

Polygon = list[list[float]]   # [[x0,y0], [x1,y1], ...]
Polygons = list[Polygon]      # list of polygons


def _to_paths(polygons: Polygons) -> list:
    return [[_tc.PointD(float(x), float(y)) for x, y in poly] for poly in polygons]


def _from_paths(paths) -> Polygons:
    return [[[p.x, p.y] for p in path] for path in paths if len(path) >= 3]


# ── Boolean operations ────────────────────────────────────────────────────────

def union(subjects: Polygons, clips: Polygons,
          fill_rule: FillRule = FillRule.NonZero) -> Polygons:
    """Return the union of *subjects* and *clips*."""
    result = _tc.union_paths(_to_paths(subjects), _to_paths(clips), fill_rule)
    return _from_paths(result)


def intersect(subjects: Polygons, clips: Polygons,
              fill_rule: FillRule = FillRule.NonZero) -> Polygons:
    """Return the intersection of *subjects* and *clips*."""
    result = _tc.intersect_paths(_to_paths(subjects), _to_paths(clips), fill_rule)
    return _from_paths(result)


def difference(subjects: Polygons, clips: Polygons,
               fill_rule: FillRule = FillRule.NonZero) -> Polygons:
    """Subtract *clips* from *subjects*."""
    result = _tc.difference_paths(_to_paths(subjects), _to_paths(clips), fill_rule)
    return _from_paths(result)


def xor(subjects: Polygons, clips: Polygons,
        fill_rule: FillRule = FillRule.NonZero) -> Polygons:
    """Return the symmetric difference of *subjects* and *clips*."""
    result = _tc.xor_paths(_to_paths(subjects), _to_paths(clips), fill_rule)
    return _from_paths(result)


# ── Inflate / offset ──────────────────────────────────────────────────────────

def inflate(polygons: Polygons, delta: float,
            join_type: JoinType = JoinType.Round,
            end_type: EndType = EndType.Polygon,
            miter_limit: float = 2.0) -> Polygons:
    """Inflate (positive *delta*) or deflate (negative *delta*) polygons."""
    result = _tc.inflate_paths(
        _to_paths(polygons), delta, join_type, end_type, miter_limit)
    return _from_paths(result)


# ── Simplify ──────────────────────────────────────────────────────────────────

def simplify(polygons: Polygons, epsilon: float, closed: bool = True) -> Polygons:
    """Simplify polygons using the Ramer–Douglas–Peucker algorithm."""
    result = _tc.simplify_paths(_to_paths(polygons), epsilon, closed)
    return _from_paths(result)


# ── TřiVis map file I/O ───────────────────────────────────────────────────────

def save_trivis_map(border: Polygon, holes: list[Polygon], path: str) -> None:
    """Save a polygonal map to a TřiVis .txt file.

    Format (matches trivis_plus SavePolyMap):
        [SCALE]
        1.0

        [BORDER]
        x y
        ...

        [OBSTACLE]
        x y
        ...
    """
    with open(path, "w") as f:
        f.write("[SCALE]\n1.0\n")
        f.write("\n[BORDER]\n")
        for x, y in border:
            f.write(f"{x:.17g} {y:.17g}\n")
        for hole in holes:
            f.write("\n[OBSTACLE]\n")
            for x, y in hole:
                f.write(f"{x:.17g} {y:.17g}\n")


def load_trivis_map(path: str) -> tuple[Polygon, list[Polygon]]:
    """Load a TřiVis .txt map file.

    Returns *(border, holes)* where both are lists of [x, y] pairs.
    """
    border: Polygon = []
    holes: list[Polygon] = []
    current: Polygon | None = None

    with open(path) as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            if line == "[SCALE]":
                current = None
                continue
            if line == "[BORDER]":
                current = border
                continue
            if line == "[OBSTACLE]":
                current = []
                holes.append(current)
                continue
            # Try to parse as coordinate pair
            if current is not None:
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        current.append([float(parts[0]), float(parts[1])])
                    except ValueError:
                        pass  # skip non-numeric tokens like section headers

    return border, holes
