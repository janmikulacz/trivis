#pragma once

#include "trivis/types.h"

#include <clipper2/clipper.h>

#include <cstddef>

namespace trivis {

// Layout assertions — caught at compile time if Clipper2 is built with USINGZ
// or if field ordering ever diverges.
static_assert(sizeof(PointInt) == sizeof(Clipper2Lib::Point64));
static_assert(offsetof(PointInt, x) == offsetof(Clipper2Lib::Point64, x));
static_assert(offsetof(PointInt, y) == offsetof(Clipper2Lib::Point64, y));

static_assert(sizeof(RectInt) == sizeof(Clipper2Lib::Rect64));
static_assert(offsetof(RectInt, left)   == offsetof(Clipper2Lib::Rect64, left));
static_assert(offsetof(RectInt, top)    == offsetof(Clipper2Lib::Rect64, top));
static_assert(offsetof(RectInt, right)  == offsetof(Clipper2Lib::Rect64, right));
static_assert(offsetof(RectInt, bottom) == offsetof(Clipper2Lib::Rect64, bottom));

// std::vector layout: same size regardless of element type (pointer + two size_t).
static_assert(sizeof(PathInt)  == sizeof(Clipper2Lib::Path64));
static_assert(sizeof(PathsInt) == sizeof(Clipper2Lib::Paths64));

// ── PointInt ↔ Point64 ───────────────────────────────────────────────────────

inline const Clipper2Lib::Point64& ToClipper(const PointInt& p) {
    return reinterpret_cast<const Clipper2Lib::Point64&>(p);
}

inline const PointInt& FromClipper(const Clipper2Lib::Point64& p) {
    return reinterpret_cast<const PointInt&>(p);
}

// ── PathInt ↔ Path64 ─────────────────────────────────────────────────────────

inline const Clipper2Lib::Path64& ToClipper(const PathInt& path) {
    return reinterpret_cast<const Clipper2Lib::Path64&>(path);
}

inline const PathInt& FromClipper(const Clipper2Lib::Path64& path) {
    return reinterpret_cast<const PathInt&>(path);
}

// ── PathsInt ↔ Paths64 ───────────────────────────────────────────────────────

inline const Clipper2Lib::Paths64& ToClipper(const PathsInt& paths) {
    return reinterpret_cast<const Clipper2Lib::Paths64&>(paths);
}

inline const PathsInt& FromClipper(const Clipper2Lib::Paths64& paths) {
    return reinterpret_cast<const PathsInt&>(paths);
}

// ── RectInt ↔ Rect64 ─────────────────────────────────────────────────────────

inline const Clipper2Lib::Rect64& ToClipper(const RectInt& r) {
    return reinterpret_cast<const Clipper2Lib::Rect64&>(r);
}

inline const RectInt& FromClipper(const Clipper2Lib::Rect64& r) {
    return reinterpret_cast<const RectInt&>(r);
}

// ── PolyTreeInt ← PolyPath64 ─────────────────────────────────────────────────
// PolyPath64 has a base class with parent_ pointer and private fields in a
// different order, so layout reinterpretation is not safe here.
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
