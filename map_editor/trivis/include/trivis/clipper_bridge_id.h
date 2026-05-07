#pragma once

// This header must be compiled in a translation unit that links Clipper2Z,
// which defines USINGZ as a public compile definition. Do NOT include it
// in the same translation unit as clipper_bridge.h.

#include "trivis/types.h"

#include <clipper2/clipper.h>

#include <cstddef>

namespace trivis {

// Layout assertions — the USINGZ build of Point64 must be exactly {x, y, z}.
static_assert(sizeof(PointIntId) == sizeof(Clipper2Lib::Point64));
static_assert(offsetof(PointIntId, x)  == offsetof(Clipper2Lib::Point64, x));
static_assert(offsetof(PointIntId, y)  == offsetof(Clipper2Lib::Point64, y));
static_assert(offsetof(PointIntId, id) == offsetof(Clipper2Lib::Point64, z));

static_assert(sizeof(PathIntId)  == sizeof(Clipper2Lib::Path64));
static_assert(sizeof(PathsIntId) == sizeof(Clipper2Lib::Paths64));

// ── PointIntId ↔ Point64 (USINGZ) ────────────────────────────────────────────

inline const Clipper2Lib::Point64& ToClipper(const PointIntId& p) {
    return reinterpret_cast<const Clipper2Lib::Point64&>(p);
}

inline const PointIntId& FromClipper(const Clipper2Lib::Point64& p) {
    return reinterpret_cast<const PointIntId&>(p);
}

// ── PathIntId ↔ Path64 (USINGZ) ──────────────────────────────────────────────

inline const Clipper2Lib::Path64& ToClipper(const PathIntId& path) {
    return reinterpret_cast<const Clipper2Lib::Path64&>(path);
}

inline const PathIntId& FromClipper(const Clipper2Lib::Path64& path) {
    return reinterpret_cast<const PathIntId&>(path);
}

// ── PathsIntId ↔ Paths64 (USINGZ) ────────────────────────────────────────────

inline const Clipper2Lib::Paths64& ToClipper(const PathsIntId& paths) {
    return reinterpret_cast<const Clipper2Lib::Paths64&>(paths);
}

inline const PathsIntId& FromClipper(const Clipper2Lib::Paths64& paths) {
    return reinterpret_cast<const PathsIntId&>(paths);
}

// ── PolyTreeIntId ← PolyPath64 (USINGZ) ──────────────────────────────────────
// PolyPath64 layout is incompatible regardless of USINGZ; recursive copy required.

inline PolyTreeIntId FromClipper(const Clipper2Lib::PolyPath64& node) {
    PolyTreeIntId result;
    result.polygon = FromClipper(node.Polygon());
    for (const auto& child : node)
        result.children.push_back(
            std::make_unique<PolyTreeIntId>(FromClipper(*child)));
    return result;
}

}  // namespace trivis
