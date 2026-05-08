#pragma once

#include "trivis/types.h"

#include <memory>

namespace trivis {

// ── Enums ─────────────────────────────────────────────────────────────────────

enum class FillRule { EvenOdd, NonZero, Positive, Negative };
enum class ClipType { None, Intersection, Union, Difference, Xor };
enum class JoinType { Square, Round, Miter };
enum class EndType  { Polygon, Joined, Butt, Square, Round };
enum class PointInPolygonResult { IsOn, IsInside, IsOutside };

// ── Free functions ────────────────────────────────────────────────────────────

PathsInt Union     (const PathsInt& subjects, const PathsInt& clips,
                    FillRule fill_rule = FillRule::NonZero);
PathsInt Intersect (const PathsInt& subjects, const PathsInt& clips,
                    FillRule fill_rule = FillRule::NonZero);
PathsInt Difference(const PathsInt& subjects, const PathsInt& clips,
                    FillRule fill_rule = FillRule::NonZero);
PathsInt Xor       (const PathsInt& subjects, const PathsInt& clips,
                    FillRule fill_rule = FillRule::NonZero);

PathsInt InflatePaths(const PathsInt& paths, double delta,
                      JoinType join_type, EndType end_type,
                      double miter_limit = 2.0, double arc_tolerance = 0.0);

PathsInt SimplifyPaths(const PathsInt& paths, double epsilon,
                       bool is_open = false);

double Area(const PathInt& path);
double Area(const PathsInt& paths);
bool   IsPositive(const PathInt& path);

PointInPolygonResult PointInPolygon(const PointInt& pt, const PathInt& polygon);

// ── Clipper (stateful clipping) ───────────────────────────────────────────────

class Clipper {
    struct Impl;
    std::unique_ptr<Impl> impl_;
public:
    Clipper();
    ~Clipper();
    Clipper(Clipper&&) noexcept;
    Clipper& operator=(Clipper&&) noexcept;

    void AddSubject    (const PathsInt& subjects);
    void AddOpenSubject(const PathsInt& subjects);
    void AddClip       (const PathsInt& clips);

    bool Execute(ClipType, FillRule, PathsInt& solution);
    bool Execute(ClipType, FillRule, PathsInt& closed, PathsInt& open);
    bool Execute(ClipType, FillRule, PolyTreeInt&);
    bool Execute(ClipType, FillRule, PolyTreeInt&, PathsInt& open);

    void Clear();
};

// ── ClipperOffset (inflate / deflate) ────────────────────────────────────────

class ClipperOffset {
    struct Impl;
    std::unique_ptr<Impl> impl_;
public:
    explicit ClipperOffset(double miter_limit = 2.0, double arc_tolerance = 0.0,
                           bool preserve_collinear = false,
                           bool reverse_solution  = false);
    ~ClipperOffset();
    ClipperOffset(ClipperOffset&&) noexcept;
    ClipperOffset& operator=(ClipperOffset&&) noexcept;

    void AddPath (const PathInt&  path,  JoinType, EndType);
    void AddPaths(const PathsInt& paths, JoinType, EndType);

    PathsInt    Execute    (double delta);
    PolyTreeInt ExecuteTree(double delta);

    double MiterLimit() const;
    void   MiterLimit(double v);
    double ArcTolerance() const;
    void   ArcTolerance(double v);

    void Clear();
};

}  // namespace trivis
