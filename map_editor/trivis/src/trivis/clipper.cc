#include "trivis/clipper.h"
#include "trivis/clipper_bridge.h"

#include <clipper2/clipper.h>

namespace trivis {

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
    PathsInt result;
    const auto& simplified = Clipper2Lib::SimplifyPaths(ToClipper(paths), epsilon, is_open);
    result = FromClipper(simplified);
    return result;
}

double Area(const PathInt& path)   { return Clipper2Lib::Area(ToClipper(path)); }
double Area(const PathsInt& paths) { return Clipper2Lib::Area(ToClipper(paths)); }
bool   IsPositive(const PathInt& path) { return Clipper2Lib::IsPositive(ToClipper(path)); }

PointInPolygonResult PointInPolygon(const PointInt& pt, const PathInt& polygon) {
    return static_cast<PointInPolygonResult>(
        Clipper2Lib::PointInPolygon(ToClipper(pt), ToClipper(polygon)));
}

// ── Clipper ───────────────────────────────────────────────────────────────────

struct Clipper::Impl {
    Clipper2Lib::Clipper64 clipper;
};

Clipper::Clipper()  : impl_(std::make_unique<Impl>()) {}
Clipper::~Clipper() = default;
Clipper::Clipper(Clipper&&) noexcept            = default;
Clipper& Clipper::operator=(Clipper&&) noexcept = default;

void Clipper::AddSubject    (const PathsInt& s) { impl_->clipper.AddSubject(ToClipper(s)); }
void Clipper::AddOpenSubject(const PathsInt& s) { impl_->clipper.AddOpenSubject(ToClipper(s)); }
void Clipper::AddClip       (const PathsInt& c) { impl_->clipper.AddClip(ToClipper(c)); }

bool Clipper::Execute(ClipType ct, FillRule fr, PathsInt& solution) {
    return impl_->clipper.Execute(C(ct), C(fr),
                                  reinterpret_cast<Clipper2Lib::Paths64&>(solution));
}

bool Clipper::Execute(ClipType ct, FillRule fr, PathsInt& closed, PathsInt& open) {
    return impl_->clipper.Execute(C(ct), C(fr),
                                  reinterpret_cast<Clipper2Lib::Paths64&>(closed),
                                  reinterpret_cast<Clipper2Lib::Paths64&>(open));
}

bool Clipper::Execute(ClipType ct, FillRule fr, PolyTreeInt& tree) {
    Clipper2Lib::PolyTree64 polytree;
    bool ok = impl_->clipper.Execute(C(ct), C(fr), polytree);
    if (ok) tree = FromClipper(polytree);
    return ok;
}

bool Clipper::Execute(ClipType ct, FillRule fr, PolyTreeInt& tree, PathsInt& open) {
    Clipper2Lib::PolyTree64 polytree;
    bool ok = impl_->clipper.Execute(C(ct), C(fr), polytree,
                                     reinterpret_cast<Clipper2Lib::Paths64&>(open));
    if (ok) tree = FromClipper(polytree);
    return ok;
}

void Clipper::Clear() { impl_->clipper.Clear(); }

// ── ClipperOffset ─────────────────────────────────────────────────────────────

struct ClipperOffset::Impl {
    Clipper2Lib::ClipperOffset offset;
    Impl(double ml, double at, bool pc, bool rs) : offset(ml, at, pc, rs) {}
};

ClipperOffset::ClipperOffset(double ml, double at, bool pc, bool rs)
    : impl_(std::make_unique<Impl>(ml, at, pc, rs)) {}
ClipperOffset::~ClipperOffset() = default;
ClipperOffset::ClipperOffset(ClipperOffset&&) noexcept            = default;
ClipperOffset& ClipperOffset::operator=(ClipperOffset&&) noexcept = default;

void ClipperOffset::AddPath(const PathInt& path, JoinType jt, EndType et) {
    impl_->offset.AddPath(ToClipper(path), C(jt), C(et));
}

void ClipperOffset::AddPaths(const PathsInt& paths, JoinType jt, EndType et) {
    impl_->offset.AddPaths(ToClipper(paths), C(jt), C(et));
}

PathsInt ClipperOffset::Execute(double delta) {
    PathsInt result;
    impl_->offset.Execute(delta, reinterpret_cast<Clipper2Lib::Paths64&>(result));
    return result;
}

PolyTreeInt ClipperOffset::ExecuteTree(double delta) {
    Clipper2Lib::PolyTree64 polytree;
    impl_->offset.Execute(delta, polytree);
    return FromClipper(polytree);
}

double ClipperOffset::MiterLimit() const  { return impl_->offset.MiterLimit(); }
void   ClipperOffset::MiterLimit(double v){ impl_->offset.MiterLimit(v); }
double ClipperOffset::ArcTolerance() const  { return impl_->offset.ArcTolerance(); }
void   ClipperOffset::ArcTolerance(double v){ impl_->offset.ArcTolerance(v); }

void ClipperOffset::Clear() { impl_->offset.Clear(); }

}  // namespace trivis
