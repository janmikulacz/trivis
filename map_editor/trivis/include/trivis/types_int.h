#pragma once

#include <cstdint>
#include <vector>
#include <memory>

namespace trivis {

using Int = int64_t;

// ── Integer geometric types ───────────────────────────────────────────────────

struct PointInt { Int x; Int y; };

using PathInt  = std::vector<PointInt>;
using PathsInt = std::vector<PathInt>;

struct RectInt { Int left; Int top; Int right; Int bottom; };

struct PolyTreeInt {
    PathInt polygon;
    std::vector<std::unique_ptr<PolyTreeInt>> children;
};

// ── Integer geometric types with vertex IDs ───────────────────────────────────

struct PointIntId { Int x; Int y; Int id; };

using PathIntId  = std::vector<PointIntId>;
using PathsIntId = std::vector<PathIntId>;

struct PolyTreeIntId {
    PathIntId polygon;
    std::vector<std::unique_ptr<PolyTreeIntId>> children;
};

}  // namespace trivis
