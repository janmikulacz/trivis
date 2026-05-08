#pragma once

#include <cstdint>
#include <vector>
#include <memory>

namespace trivis {

using Float = double;

// ── Float geometric types ─────────────────────────────────────────────────────

struct PointFloat { Float x; Float y; };

using PathFloat  = std::vector<PointFloat>;
using PathsFloat = std::vector<PathFloat>;

struct RectFloat { Float left; Float top; Float right; Float bottom; };

struct PolyTreeFloat {
    PathFloat polygon;
    std::vector<std::unique_ptr<PolyTreeFloat>> children;
};

// ── Float geometric types with vertex IDs ────────────────────────────────────

struct PointFloatId { Float x; Float y; int64_t id; };

using PathFloatId  = std::vector<PointFloatId>;
using PathsFloatId = std::vector<PathFloatId>;

struct PolyTreeFloatId {
    PathFloatId polygon;
    std::vector<std::unique_ptr<PolyTreeFloatId>> children;
};

}  // namespace trivis
