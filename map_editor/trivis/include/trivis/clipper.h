#pragma once

#include "trivis/types.h"
#include "trivis/clipper_bridge.h"

#include <clipper2/clipper.h>

namespace trivis {

// ── Enum aliases ──────────────────────────────────────────────────────────────

using FillRule            = Clipper2Lib::FillRule;
using ClipType            = Clipper2Lib::ClipType;
using JoinType            = Clipper2Lib::JoinType;
using EndType             = Clipper2Lib::EndType;
using PointInPolygonResult = Clipper2Lib::PointInPolygonResult;

// ── Boolean operations (free functions) ──────────────────────────────────────

inline PathsInt Union(const PathsInt& subjects, const PathsInt& clips,
                      FillRule fill_rule = FillRule::NonZero) {
    PathsInt result;
    Clipper2Lib::Clipper64 c;
    c.AddSubject(ToClipper(subjects));
    c.AddClip(ToClipper(clips));
    c.Execute(ClipType::Union, fill_rule,
              reinterpret_cast<Clipper2Lib::Paths64&>(result));
    return result;
}

inline PathsInt Intersect(const PathsInt& subjects, const PathsInt& clips,
                           FillRule fill_rule = FillRule::NonZero) {
    PathsInt result;
    Clipper2Lib::Clipper64 c;
    c.AddSubject(ToClipper(subjects));
    c.AddClip(ToClipper(clips));
    c.Execute(ClipType::Intersection, fill_rule,
              reinterpret_cast<Clipper2Lib::Paths64&>(result));
    return result;
}

inline PathsInt Difference(const PathsInt& subjects, const PathsInt& clips,
                            FillRule fill_rule = FillRule::NonZero) {
    PathsInt result;
    Clipper2Lib::Clipper64 c;
    c.AddSubject(ToClipper(subjects));
    c.AddClip(ToClipper(clips));
    c.Execute(ClipType::Difference, fill_rule,
              reinterpret_cast<Clipper2Lib::Paths64&>(result));
    return result;
}

inline PathsInt Xor(const PathsInt& subjects, const PathsInt& clips,
                    FillRule fill_rule = FillRule::NonZero) {
    PathsInt result;
    Clipper2Lib::Clipper64 c;
    c.AddSubject(ToClipper(subjects));
    c.AddClip(ToClipper(clips));
    c.Execute(ClipType::Xor, fill_rule,
              reinterpret_cast<Clipper2Lib::Paths64&>(result));
    return result;
}

inline PathsInt InflatePaths(const PathsInt& paths, double delta,
                              JoinType join_type, EndType end_type,
                              double miter_limit = 2.0,
                              double arc_tolerance = 0.0) {
    PathsInt result;
    Clipper2Lib::ClipperOffset co(miter_limit, arc_tolerance);
    co.AddPaths(ToClipper(paths), join_type, end_type);
    co.Execute(delta, reinterpret_cast<Clipper2Lib::Paths64&>(result));
    return result;
}

inline PathsInt SimplifyPaths(const PathsInt& paths, double epsilon,
                               bool is_open = false) {
    return FromClipper(Clipper2Lib::SimplifyPaths(ToClipper(paths), epsilon, is_open));
}

// ── Geometry utilities ────────────────────────────────────────────────────────

inline double Area(const PathInt& path) {
    return Clipper2Lib::Area(ToClipper(path));
}

inline double Area(const PathsInt& paths) {
    return Clipper2Lib::Area(ToClipper(paths));
}

inline bool IsPositive(const PathInt& path) {
    return Clipper2Lib::IsPositive(ToClipper(path));
}

inline PointInPolygonResult PointInPolygon(const PointInt& pt,
                                           const PathInt& polygon) {
    return Clipper2Lib::PointInPolygon(ToClipper(pt), ToClipper(polygon));
}

// ── Clipper class (stateful clipping) ────────────────────────────────────────

class Clipper {
    Clipper2Lib::Clipper64 impl_;
public:
    void AddSubject(const PathsInt& subjects) {
        impl_.AddSubject(ToClipper(subjects));
    }

    void AddOpenSubject(const PathsInt& subjects) {
        impl_.AddOpenSubject(ToClipper(subjects));
    }

    void AddClip(const PathsInt& clips) {
        impl_.AddClip(ToClipper(clips));
    }

    bool Execute(ClipType clip_type, FillRule fill_rule, PathsInt& solution) {
        return impl_.Execute(clip_type, fill_rule,
                             reinterpret_cast<Clipper2Lib::Paths64&>(solution));
    }

    bool Execute(ClipType clip_type, FillRule fill_rule,
                 PathsInt& closed_paths, PathsInt& open_paths) {
        return impl_.Execute(clip_type, fill_rule,
                             reinterpret_cast<Clipper2Lib::Paths64&>(closed_paths),
                             reinterpret_cast<Clipper2Lib::Paths64&>(open_paths));
    }

    bool Execute(ClipType clip_type, FillRule fill_rule, PolyTreeInt& tree) {
        Clipper2Lib::PolyTree64 polytree;
        bool ok = impl_.Execute(clip_type, fill_rule, polytree);
        if (ok) tree = FromClipper(polytree);
        return ok;
    }

    bool Execute(ClipType clip_type, FillRule fill_rule,
                 PolyTreeInt& tree, PathsInt& open_paths) {
        Clipper2Lib::PolyTree64 polytree;
        bool ok = impl_.Execute(clip_type, fill_rule, polytree,
                                reinterpret_cast<Clipper2Lib::Paths64&>(open_paths));
        if (ok) tree = FromClipper(polytree);
        return ok;
    }

    void Clear() { impl_.Clear(); }
};

// ── ClipperOffset class (inflate / deflate) ───────────────────────────────────

class ClipperOffset {
    Clipper2Lib::ClipperOffset impl_;
public:
    explicit ClipperOffset(double miter_limit = 2.0, double arc_tolerance = 0.0,
                           bool preserve_collinear = false,
                           bool reverse_solution = false)
        : impl_(miter_limit, arc_tolerance, preserve_collinear, reverse_solution) {}

    void AddPath(const PathInt& path, JoinType join_type, EndType end_type) {
        impl_.AddPath(ToClipper(path), join_type, end_type);
    }

    void AddPaths(const PathsInt& paths, JoinType join_type, EndType end_type) {
        impl_.AddPaths(ToClipper(paths), join_type, end_type);
    }

    PathsInt Execute(double delta) {
        PathsInt result;
        impl_.Execute(delta, reinterpret_cast<Clipper2Lib::Paths64&>(result));
        return result;
    }

    PolyTreeInt ExecuteTree(double delta) {
        Clipper2Lib::PolyTree64 polytree;
        impl_.Execute(delta, polytree);
        return FromClipper(polytree);
    }

    double MiterLimit() const    { return impl_.MiterLimit(); }
    void MiterLimit(double v)    { impl_.MiterLimit(v); }
    double ArcTolerance() const  { return impl_.ArcTolerance(); }
    void ArcTolerance(double v)  { impl_.ArcTolerance(v); }

    void Clear() { impl_.Clear(); }
};

}  // namespace trivis
