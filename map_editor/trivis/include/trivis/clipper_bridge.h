#pragma once

#include "trivis/types.h"

#include <clipper2/clipper.h>

namespace trivis {

// ── PointInt ↔ Point64 ───────────────────────────────────────────────────────

inline Clipper2Lib::Point64 ToClipper(const PointInt& p) {
    return {p.x, p.y};
}

inline PointInt FromClipper(const Clipper2Lib::Point64& p) {
    return {p.x, p.y};
}

// ── PathInt ↔ Path64 ─────────────────────────────────────────────────────────

inline Clipper2Lib::Path64 ToClipper(const PathInt& path) {
    Clipper2Lib::Path64 result;
    result.reserve(path.size());
    for (const auto& pt : path)
        result.push_back(ToClipper(pt));
    return result;
}

inline PathInt FromClipper(const Clipper2Lib::Path64& path) {
    PathInt result;
    result.reserve(path.size());
    for (const auto& pt : path)
        result.push_back(FromClipper(pt));
    return result;
}

// ── PathsInt ↔ Paths64 ───────────────────────────────────────────────────────

inline Clipper2Lib::Paths64 ToClipper(const PathsInt& paths) {
    Clipper2Lib::Paths64 result;
    result.reserve(paths.size());
    for (const auto& path : paths)
        result.push_back(ToClipper(path));
    return result;
}

inline PathsInt FromClipper(const Clipper2Lib::Paths64& paths) {
    PathsInt result;
    result.reserve(paths.size());
    for (const auto& path : paths)
        result.push_back(FromClipper(path));
    return result;
}

// ── RectInt ↔ Rect64 ─────────────────────────────────────────────────────────

inline Clipper2Lib::Rect64 ToClipper(const RectInt& r) {
    return {r.left, r.top, r.right, r.bottom};
}

inline RectInt FromClipper(const Clipper2Lib::Rect64& r) {
    return {r.left, r.top, r.right, r.bottom};
}

// ── PolyTreeInt ← PolyPath64 ─────────────────────────────────────────────────
// PolyTree64 is a read-only clipper output; there is no ToClipper counterpart.

inline PolyTreeInt FromClipper(const Clipper2Lib::PolyPath64& node) {
    PolyTreeInt result;
    result.polygon = FromClipper(node.Polygon());
    for (const auto& child : node)
        result.children.push_back(
            std::make_unique<PolyTreeInt>(FromClipper(*child)));
    return result;
}

}  // namespace trivis
