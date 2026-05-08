#include "trivis/clipper_int_id.h"

#include <clipper2/clipper.h>   // compiled with USINGZ via Clipper2Z linkage

#include <cstddef>

namespace trivis {

// ── Bridge: trivis IntId ↔ Clipper2 types (USINGZ build) ─────────────────────

static_assert(sizeof(PointIntId) == sizeof(Clipper2Lib::Point64));
static_assert(offsetof(PointIntId, x)  == offsetof(Clipper2Lib::Point64, x));
static_assert(offsetof(PointIntId, y)  == offsetof(Clipper2Lib::Point64, y));
static_assert(offsetof(PointIntId, id) == offsetof(Clipper2Lib::Point64, z));

static_assert(sizeof(PathIntId)  == sizeof(Clipper2Lib::Path64));
static_assert(sizeof(PathsIntId) == sizeof(Clipper2Lib::Paths64));

namespace {

const Clipper2Lib::Point64& ToClipper(const PointIntId& p) {
    return reinterpret_cast<const Clipper2Lib::Point64&>(p);
}
const Clipper2Lib::Path64& ToClipper(const PathIntId& p) {
    return reinterpret_cast<const Clipper2Lib::Path64&>(p);
}
const Clipper2Lib::Paths64& ToClipper(const PathsIntId& p) {
    return reinterpret_cast<const Clipper2Lib::Paths64&>(p);
}

const PathsIntId& FromClipper(const Clipper2Lib::Paths64& p) {
    return reinterpret_cast<const PathsIntId&>(p);
}

PolyTreeIntId FromClipper(const Clipper2Lib::PolyPath64& node) {
    PolyTreeIntId result;
    result.polygon = reinterpret_cast<const PathIntId&>(node.Polygon());
    for (const auto& child : node)
        result.children.push_back(
            std::make_unique<PolyTreeIntId>(FromClipper(*child)));
    return result;
}

}  // namespace

// ── Enum casts (same enums as clipper.h, same underlying values) ──────────────

static Clipper2Lib::FillRule C(FillRule v) { return static_cast<Clipper2Lib::FillRule>(v); }
static Clipper2Lib::ClipType C(ClipType v) { return static_cast<Clipper2Lib::ClipType>(v); }
static Clipper2Lib::JoinType C(JoinType v) { return static_cast<Clipper2Lib::JoinType>(v); }
static Clipper2Lib::EndType  C(EndType  v) { return static_cast<Clipper2Lib::EndType>(v); }

// ── Free functions ────────────────────────────────────────────────────────────

PathsIntId Union(const PathsIntId& subjects, const PathsIntId& clips, FillRule fill_rule) {
    PathsIntId result;
    Clipper2Lib::Clipper64 c;
    c.AddSubject(ToClipper(subjects));
    c.AddClip(ToClipper(clips));
    c.Execute(Clipper2Lib::ClipType::Union, C(fill_rule),
              reinterpret_cast<Clipper2Lib::Paths64&>(result));
    return result;
}

PathsIntId Intersect(const PathsIntId& subjects, const PathsIntId& clips, FillRule fill_rule) {
    PathsIntId result;
    Clipper2Lib::Clipper64 c;
    c.AddSubject(ToClipper(subjects));
    c.AddClip(ToClipper(clips));
    c.Execute(Clipper2Lib::ClipType::Intersection, C(fill_rule),
              reinterpret_cast<Clipper2Lib::Paths64&>(result));
    return result;
}

PathsIntId Difference(const PathsIntId& subjects, const PathsIntId& clips, FillRule fill_rule) {
    PathsIntId result;
    Clipper2Lib::Clipper64 c;
    c.AddSubject(ToClipper(subjects));
    c.AddClip(ToClipper(clips));
    c.Execute(Clipper2Lib::ClipType::Difference, C(fill_rule),
              reinterpret_cast<Clipper2Lib::Paths64&>(result));
    return result;
}

PathsIntId Xor(const PathsIntId& subjects, const PathsIntId& clips, FillRule fill_rule) {
    PathsIntId result;
    Clipper2Lib::Clipper64 c;
    c.AddSubject(ToClipper(subjects));
    c.AddClip(ToClipper(clips));
    c.Execute(Clipper2Lib::ClipType::Xor, C(fill_rule),
              reinterpret_cast<Clipper2Lib::Paths64&>(result));
    return result;
}

PathsIntId InflatePaths(const PathsIntId& paths, double delta,
                        JoinType join_type, EndType end_type,
                        double miter_limit, double arc_tolerance) {
    PathsIntId result;
    Clipper2Lib::ClipperOffset co(miter_limit, arc_tolerance);
    co.AddPaths(ToClipper(paths), C(join_type), C(end_type));
    co.Execute(delta, reinterpret_cast<Clipper2Lib::Paths64&>(result));
    return result;
}

PathsIntId SimplifyPaths(const PathsIntId& paths, double epsilon, bool is_open) {
    return FromClipper(Clipper2Lib::SimplifyPaths(ToClipper(paths), epsilon, is_open));
}

double Area     (const PathIntId& path)   { return Clipper2Lib::Area(ToClipper(path)); }
double Area     (const PathsIntId& paths) { return Clipper2Lib::Area(ToClipper(paths)); }
bool   IsPositive(const PathIntId& path)  { return Clipper2Lib::IsPositive(ToClipper(path)); }

PointInPolygonResult PointInPolygon(const PointIntId& pt, const PathIntId& polygon) {
    return static_cast<PointInPolygonResult>(
        Clipper2Lib::PointInPolygon(ToClipper(pt), ToClipper(polygon)));
}

// ── ClipperZ ──────────────────────────────────────────────────────────────────

struct ClipperZ::Impl {
    Clipper2Lib::Clipper64 clipper;
};

ClipperZ::ClipperZ()  : impl_(std::make_unique<Impl>()) {}
ClipperZ::~ClipperZ() = default;
ClipperZ::ClipperZ(ClipperZ&&) noexcept            = default;
ClipperZ& ClipperZ::operator=(ClipperZ&&) noexcept = default;

void ClipperZ::AddSubject    (const PathsIntId& s) { impl_->clipper.AddSubject(ToClipper(s)); }
void ClipperZ::AddOpenSubject(const PathsIntId& s) { impl_->clipper.AddOpenSubject(ToClipper(s)); }
void ClipperZ::AddClip       (const PathsIntId& c) { impl_->clipper.AddClip(ToClipper(c)); }

bool ClipperZ::Execute(ClipType ct, FillRule fr, PathsIntId& solution) {
    return impl_->clipper.Execute(C(ct), C(fr),
                                  reinterpret_cast<Clipper2Lib::Paths64&>(solution));
}

bool ClipperZ::Execute(ClipType ct, FillRule fr, PathsIntId& closed, PathsIntId& open) {
    return impl_->clipper.Execute(C(ct), C(fr),
                                  reinterpret_cast<Clipper2Lib::Paths64&>(closed),
                                  reinterpret_cast<Clipper2Lib::Paths64&>(open));
}

bool ClipperZ::Execute(ClipType ct, FillRule fr, PolyTreeIntId& tree) {
    Clipper2Lib::PolyTree64 polytree;
    bool ok = impl_->clipper.Execute(C(ct), C(fr), polytree);
    if (ok) tree = FromClipper(polytree);
    return ok;
}

bool ClipperZ::Execute(ClipType ct, FillRule fr, PolyTreeIntId& tree, PathsIntId& open) {
    Clipper2Lib::PolyTree64 polytree;
    bool ok = impl_->clipper.Execute(C(ct), C(fr), polytree,
                                     reinterpret_cast<Clipper2Lib::Paths64&>(open));
    if (ok) tree = FromClipper(polytree);
    return ok;
}

void ClipperZ::Clear() { impl_->clipper.Clear(); }

// ── ClipperOffsetZ ────────────────────────────────────────────────────────────

struct ClipperOffsetZ::Impl {
    Clipper2Lib::ClipperOffset offset;
    Impl(double ml, double at, bool pc, bool rs) : offset(ml, at, pc, rs) {}
};

ClipperOffsetZ::ClipperOffsetZ(double ml, double at, bool pc, bool rs)
    : impl_(std::make_unique<Impl>(ml, at, pc, rs)) {}
ClipperOffsetZ::~ClipperOffsetZ() = default;
ClipperOffsetZ::ClipperOffsetZ(ClipperOffsetZ&&) noexcept            = default;
ClipperOffsetZ& ClipperOffsetZ::operator=(ClipperOffsetZ&&) noexcept = default;

void ClipperOffsetZ::AddPath(const PathIntId& path, JoinType jt, EndType et) {
    impl_->offset.AddPath(ToClipper(path), C(jt), C(et));
}

void ClipperOffsetZ::AddPaths(const PathsIntId& paths, JoinType jt, EndType et) {
    impl_->offset.AddPaths(ToClipper(paths), C(jt), C(et));
}

PathsIntId ClipperOffsetZ::Execute(double delta) {
    PathsIntId result;
    impl_->offset.Execute(delta, reinterpret_cast<Clipper2Lib::Paths64&>(result));
    return result;
}

PolyTreeIntId ClipperOffsetZ::ExecuteTree(double delta) {
    Clipper2Lib::PolyTree64 polytree;
    impl_->offset.Execute(delta, polytree);
    return FromClipper(polytree);
}

double ClipperOffsetZ::MiterLimit() const   { return impl_->offset.MiterLimit(); }
void   ClipperOffsetZ::MiterLimit(double v) { impl_->offset.MiterLimit(v); }
double ClipperOffsetZ::ArcTolerance() const   { return impl_->offset.ArcTolerance(); }
void   ClipperOffsetZ::ArcTolerance(double v) { impl_->offset.ArcTolerance(v); }

void ClipperOffsetZ::Clear() { impl_->offset.Clear(); }

}  // namespace trivis
