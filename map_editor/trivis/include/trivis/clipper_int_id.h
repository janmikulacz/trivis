#pragma once

#include "trivis/clipper_int.h"   // shared enums: FillRule, ClipType, JoinType, EndType
#include "trivis/types_int.h"

#include <memory>

namespace trivis {

// ── Free functions (IntId variants — require Clipper2Z / USINGZ) ──────────────

PathsIntId Union     (const PathsIntId& subjects, const PathsIntId& clips,
                      FillRule fill_rule = FillRule::NonZero);
PathsIntId Intersect (const PathsIntId& subjects, const PathsIntId& clips,
                      FillRule fill_rule = FillRule::NonZero);
PathsIntId Difference(const PathsIntId& subjects, const PathsIntId& clips,
                      FillRule fill_rule = FillRule::NonZero);
PathsIntId Xor       (const PathsIntId& subjects, const PathsIntId& clips,
                      FillRule fill_rule = FillRule::NonZero);

PathsIntId InflatePaths(const PathsIntId& paths, double delta,
                        JoinType join_type, EndType end_type,
                        double miter_limit = 2.0, double arc_tolerance = 0.0);

PathsIntId SimplifyPaths(const PathsIntId& paths, double epsilon,
                         bool is_open = false);

double Area     (const PathIntId& path);
double Area     (const PathsIntId& paths);
bool   IsPositive(const PathIntId& path);

PointInPolygonResult PointInPolygon(const PointIntId& pt, const PathIntId& polygon);

// ── ClipperZ (stateful clipping with vertex IDs) ──────────────────────────────

class ClipperZ {
    struct Impl;
    std::unique_ptr<Impl> impl_;
public:
    ClipperZ();
    ~ClipperZ();
    ClipperZ(ClipperZ&&) noexcept;
    ClipperZ& operator=(ClipperZ&&) noexcept;

    void AddSubject    (const PathsIntId& subjects);
    void AddOpenSubject(const PathsIntId& subjects);
    void AddClip       (const PathsIntId& clips);

    bool Execute(ClipType, FillRule, PathsIntId& solution);
    bool Execute(ClipType, FillRule, PathsIntId& closed, PathsIntId& open);
    bool Execute(ClipType, FillRule, PolyTreeIntId&);
    bool Execute(ClipType, FillRule, PolyTreeIntId&, PathsIntId& open);

    void Clear();
};

// ── ClipperOffsetZ (inflate / deflate with vertex IDs) ───────────────────────

class ClipperOffsetZ {
    struct Impl;
    std::unique_ptr<Impl> impl_;
public:
    explicit ClipperOffsetZ(double miter_limit = 2.0, double arc_tolerance = 0.0,
                            bool preserve_collinear = false,
                            bool reverse_solution   = false);
    ~ClipperOffsetZ();
    ClipperOffsetZ(ClipperOffsetZ&&) noexcept;
    ClipperOffsetZ& operator=(ClipperOffsetZ&&) noexcept;

    void AddPath (const PathIntId&  path,  JoinType, EndType);
    void AddPaths(const PathsIntId& paths, JoinType, EndType);

    PathsIntId    Execute    (double delta);
    PolyTreeIntId ExecuteTree(double delta);

    double MiterLimit() const;
    void   MiterLimit(double v);
    double ArcTolerance() const;
    void   ArcTolerance(double v);

    void Clear();
};

}  // namespace trivis
