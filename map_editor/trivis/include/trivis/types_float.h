#pragma once

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

}  // namespace trivis
