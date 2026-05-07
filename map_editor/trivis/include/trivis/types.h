#pragma once

#include <cstdint>
#include <vector>
#include <memory>

namespace trivis {

// ── Scalar types ──────────────────────────────────────────────────────────────

using Int   = int64_t;
using Float = double;

// ── Integer geometric types ───────────────────────────────────────────────────

struct PointInt {
    Int x;
    Int y;
};

using PathInt  = std::vector<PointInt>;
using PathsInt = std::vector<PathInt>;

struct RectInt {
    Int left;
    Int top;
    Int right;
    Int bottom;
};

struct PolyTreeInt {
    PathInt polygon;
    std::vector<std::unique_ptr<PolyTreeInt>> children;
};

// ── Integer geometric types with vertex IDs (maps to Clipper2Z / USINGZ) ─────

struct PointIntId {
    Int x;
    Int y;
    Int id;
};

using PathIntId  = std::vector<PointIntId>;
using PathsIntId = std::vector<PathIntId>;

struct PolyTreeIntId {
    PathIntId polygon;
    std::vector<std::unique_ptr<PolyTreeIntId>> children;
};

// ── Float geometric types ─────────────────────────────────────────────────────

struct PointFloat {
    Float x;
    Float y;
};

using PathFloat  = std::vector<PointFloat>;
using PathsFloat = std::vector<PathFloat>;

struct RectFloat {
    Float left;
    Float top;
    Float right;
    Float bottom;
};

struct PolyTreeFloat {
    PathFloat polygon;
    std::vector<std::unique_ptr<PolyTreeFloat>> children;
};

}  // namespace trivis
