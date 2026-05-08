#include "trivis/clipper_int.h"

#include <clipper2/clipper.h>

#include <cstddef>

namespace trivis {

// ── Bridge: trivis ↔ Clipper2 types ──────────────────────────────────────────

static_assert(sizeof(PointInt) == sizeof(Clipper2Lib::Point64));
static_assert(offsetof(PointInt, x) == offsetof(Clipper2Lib::Point64, x));
static_assert(offsetof(PointInt, y) == offsetof(Clipper2Lib::Point64, y));

static_assert(sizeof(RectInt) == sizeof(Clipper2Lib::Rect64));
static_assert(offsetof(RectInt, left)   == offsetof(Clipper2Lib::Rect64, left));
static_assert(offsetof(RectInt, top)    == offsetof(Clipper2Lib::Rect64, top));
static_assert(offsetof(RectInt, right)  == offsetof(Clipper2Lib::Rect64, right));
static_assert(offsetof(RectInt, bottom) == offsetof(Clipper2Lib::Rect64, bottom));

static_assert(sizeof(PathInt)  == sizeof(Clipper2Lib::Path64));
static_assert(sizeof(PathsInt) == sizeof(Clipper2Lib::Paths64));

namespace {

const Clipper2Lib::Point64& ToClipper(const PointInt& p) {
    return reinterpret_cast<const Clipper2Lib::Point64&>(p);
}
const Clipper2Lib::Path64& ToClipper(const PathInt& p) {
    return reinterpret_cast<const Clipper2Lib::Path64&>(p);
}
const Clipper2Lib::Paths64& ToClipper(const PathsInt& p) {
    return reinterpret_cast<const Clipper2Lib::Paths64&>(p);
}
const Clipper2Lib::Rect64& ToClipper(const RectInt& r) {
    return reinterpret_cast<const Clipper2Lib::Rect64&>(r);
}

const PathsInt& FromClipper(const Clipper2Lib::Paths64& p) {
    return reinterpret_cast<const PathsInt&>(p);
}

PolyTreeInt FromClipper(const Clipper2Lib::PolyPath64& node) {
    PolyTreeInt result;
    result.polygon = reinterpret_cast<const PathInt&>(node.Polygon());
    for (const auto& child : node)
        result.children.push_back(
            std::make_unique<PolyTreeInt>(FromClipper(*child)));
    return result;
}

}  // namespace

// ── Enum value assertions ─────────────────────────────────────────────────────

static_assert(static_cast<int>(FillRule::EvenOdd)  == static_cast<int>(Clipper2Lib::FillRule::EvenOdd));
static_assert(static_cast<int>(FillRule::NonZero)  == static_cast<int>(Clipper2Lib::FillRule::NonZero));
static_assert(static_cast<int>(FillRule::Positive) == static_cast<int>(Clipper2Lib::FillRule::Positive));
static_assert(static_cast<int>(FillRule::Negative) == static_cast<int>(Clipper2Lib::FillRule::Negative));

static_assert(static_cast<int>(ClipType::None)         == static_cast<int>(Clipper2Lib::ClipType::None));
static_assert(static_cast<int>(ClipType::Intersection) == static_cast<int>(Clipper2Lib::ClipType::Intersection));
static_assert(static_cast<int>(ClipType::Union)        == static_cast<int>(Clipper2Lib::ClipType::Union));
static_assert(static_cast<int>(ClipType::Difference)   == static_cast<int>(Clipper2Lib::ClipType::Difference));
static_assert(static_cast<int>(ClipType::Xor)          == static_cast<int>(Clipper2Lib::ClipType::Xor));

static_assert(static_cast<int>(JoinType::Square) == static_cast<int>(Clipper2Lib::JoinType::Square));
static_assert(static_cast<int>(JoinType::Round)  == static_cast<int>(Clipper2Lib::JoinType::Round));
static_assert(static_cast<int>(JoinType::Miter)  == static_cast<int>(Clipper2Lib::JoinType::Miter));

static_assert(static_cast<int>(EndType::Polygon) == static_cast<int>(Clipper2Lib::EndType::Polygon));
static_assert(static_cast<int>(EndType::Joined)  == static_cast<int>(Clipper2Lib::EndType::Joined));
static_assert(static_cast<int>(EndType::Butt)    == static_cast<int>(Clipper2Lib::EndType::Butt));
static_assert(static_cast<int>(EndType::Square)  == static_cast<int>(Clipper2Lib::EndType::Square));
static_assert(static_cast<int>(EndType::Round)   == static_cast<int>(Clipper2Lib::EndType::Round));

static_assert(static_cast<int>(PointInPolygonResult::IsOn)      == static_cast<int>(Clipper2Lib::PointInPolygonResult::IsOn));
static_assert(static_cast<int>(PointInPolygonResult::IsInside)  == static_cast<int>(Clipper2Lib::PointInPolygonResult::IsInside));
static_assert(static_cast<int>(PointInPolygonResult::IsOutside) == static_cast<int>(Clipper2Lib::PointInPolygonResult::IsOutside));

// ── Enum casts ────────────────────────────────────────────────────────────────

static Clipper2Lib::FillRule C(FillRule v) { return static_cast<Clipper2Lib::FillRule>(v); }
static Clipper2Lib::ClipType C(ClipType v) { return static_cast<Clipper2Lib::ClipType>(v); }
static Clipper2Lib::JoinType C(JoinType v) { return static_cast<Clipper2Lib::JoinType>(v); }
static Clipper2Lib::EndType  C(EndType  v) { return static_cast<Clipper2Lib::EndType>(v); }

// ── Free functions ────────────────────────────────────────────────────────────

PathsInt Union(const PathsInt& subjects, const PathsInt& clips, FillRule fill_rule) {
    PathsInt result;
    Clipper2Lib::Clipper64 c;
    c.AddSubject(ToClipper(subjects));
    c.AddClip(ToClipper(clips));
    c.Execute(Clipper2Lib::ClipType::Union, C(fill_rule),
              reinterpret_cast<Clipper2Lib::Paths64&>(result));
    return result;
}

PathsInt Intersect(const PathsInt& subjects, const PathsInt& clips, FillRule fill_rule) {
    PathsInt result;
    Clipper2Lib::Clipper64 c;
    c.AddSubject(ToClipper(subjects));
    c.AddClip(ToClipper(clips));
    c.Execute(Clipper2Lib::ClipType::Intersection, C(fill_rule),
              reinterpret_cast<Clipper2Lib::Paths64&>(result));
    return result;
}

PathsInt Difference(const PathsInt& subjects, const PathsInt& clips, FillRule fill_rule) {
    PathsInt result;
    Clipper2Lib::Clipper64 c;
    c.AddSubject(ToClipper(subjects));
    c.AddClip(ToClipper(clips));
    c.Execute(Clipper2Lib::ClipType::Difference, C(fill_rule),
              reinterpret_cast<Clipper2Lib::Paths64&>(result));
    return result;
}

PathsInt Xor(const PathsInt& subjects, const PathsInt& clips, FillRule fill_rule) {
    PathsInt result;
    Clipper2Lib::Clipper64 c;
    c.AddSubject(ToClipper(subjects));
    c.AddClip(ToClipper(clips));
    c.Execute(Clipper2Lib::ClipType::Xor, C(fill_rule),
              reinterpret_cast<Clipper2Lib::Paths64&>(result));
    return result;
}

PathsInt InflatePaths(const PathsInt& paths, double delta,
                      JoinType join_type, EndType end_type,
                      double miter_limit, double arc_tolerance) {
    PathsInt result;
    Clipper2Lib::ClipperOffset co(miter_limit, arc_tolerance);
    co.AddPaths(ToClipper(paths), C(join_type), C(end_type));
    co.Execute(delta, reinterpret_cast<Clipper2Lib::Paths64&>(result));
    return result;
}

PathsInt SimplifyPaths(const PathsInt& paths, double epsilon, bool is_open) {
    return FromClipper(Clipper2Lib::SimplifyPaths(ToClipper(paths), epsilon, is_open));
}

double Area(const PathInt& path)       { return Clipper2Lib::Area(ToClipper(path)); }
double Area(const PathsInt& paths)     { return Clipper2Lib::Area(ToClipper(paths)); }
bool   IsPositive(const PathInt& path) { return Clipper2Lib::IsPositive(ToClipper(path)); }

PointInPolygonResult PointInPolygon(const PointInt& pt, const PathInt& polygon) {
    return static_cast<PointInPolygonResult>(
        Clipper2Lib::PointInPolygon(ToClipper(pt), ToClipper(polygon)));
}

// ── ClipperInt ───────────────────────────────────────────────────────────────────

struct ClipperInt::Impl {
    Clipper2Lib::Clipper64 clipper;
};

ClipperInt::ClipperInt()  : impl_(std::make_unique<Impl>()) {}
ClipperInt::~ClipperInt() = default;
ClipperInt::ClipperInt(ClipperInt&&) noexcept            = default;
ClipperInt& ClipperInt::operator=(ClipperInt&&) noexcept = default;

void ClipperInt::AddSubject    (const PathsInt& s) { impl_->clipper.AddSubject(ToClipper(s)); }
void ClipperInt::AddOpenSubject(const PathsInt& s) { impl_->clipper.AddOpenSubject(ToClipper(s)); }
void ClipperInt::AddClip       (const PathsInt& c) { impl_->clipper.AddClip(ToClipper(c)); }

bool ClipperInt::Execute(ClipType ct, FillRule fr, PathsInt& solution) {
    return impl_->clipper.Execute(C(ct), C(fr),
                                  reinterpret_cast<Clipper2Lib::Paths64&>(solution));
}

bool ClipperInt::Execute(ClipType ct, FillRule fr, PathsInt& closed, PathsInt& open) {
    return impl_->clipper.Execute(C(ct), C(fr),
                                  reinterpret_cast<Clipper2Lib::Paths64&>(closed),
                                  reinterpret_cast<Clipper2Lib::Paths64&>(open));
}

bool ClipperInt::Execute(ClipType ct, FillRule fr, PolyTreeInt& tree) {
    Clipper2Lib::PolyTree64 polytree;
    bool ok = impl_->clipper.Execute(C(ct), C(fr), polytree);
    if (ok) tree = FromClipper(polytree);
    return ok;
}

bool ClipperInt::Execute(ClipType ct, FillRule fr, PolyTreeInt& tree, PathsInt& open) {
    Clipper2Lib::PolyTree64 polytree;
    bool ok = impl_->clipper.Execute(C(ct), C(fr), polytree,
                                     reinterpret_cast<Clipper2Lib::Paths64&>(open));
    if (ok) tree = FromClipper(polytree);
    return ok;
}

void ClipperInt::Clear() { impl_->clipper.Clear(); }

// ── ClipperOffsetInt ─────────────────────────────────────────────────────────────

struct ClipperOffsetInt::Impl {
    Clipper2Lib::ClipperOffset offset;
    Impl(double ml, double at, bool pc, bool rs) : offset(ml, at, pc, rs) {}
};

ClipperOffsetInt::ClipperOffsetInt(double ml, double at, bool pc, bool rs)
    : impl_(std::make_unique<Impl>(ml, at, pc, rs)) {}
ClipperOffsetInt::~ClipperOffsetInt() = default;
ClipperOffsetInt::ClipperOffsetInt(ClipperOffsetInt&&) noexcept            = default;
ClipperOffsetInt& ClipperOffsetInt::operator=(ClipperOffsetInt&&) noexcept = default;

void ClipperOffsetInt::AddPath(const PathInt& path, JoinType jt, EndType et) {
    impl_->offset.AddPath(ToClipper(path), C(jt), C(et));
}

void ClipperOffsetInt::AddPaths(const PathsInt& paths, JoinType jt, EndType et) {
    impl_->offset.AddPaths(ToClipper(paths), C(jt), C(et));
}

PathsInt ClipperOffsetInt::Execute(double delta) {
    PathsInt result;
    impl_->offset.Execute(delta, reinterpret_cast<Clipper2Lib::Paths64&>(result));
    return result;
}

PolyTreeInt ClipperOffsetInt::ExecuteTree(double delta) {
    Clipper2Lib::PolyTree64 polytree;
    impl_->offset.Execute(delta, polytree);
    return FromClipper(polytree);
}

double ClipperOffsetInt::MiterLimit() const   { return impl_->offset.MiterLimit(); }
void   ClipperOffsetInt::MiterLimit(double v) { impl_->offset.MiterLimit(v); }
double ClipperOffsetInt::ArcTolerance() const   { return impl_->offset.ArcTolerance(); }
void   ClipperOffsetInt::ArcTolerance(double v) { impl_->offset.ArcTolerance(v); }

void ClipperOffsetInt::Clear() { impl_->offset.Clear(); }

}  // namespace trivis
