#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "clipper2/clipper.h"

namespace py = pybind11;
using namespace Clipper2Lib;

PYBIND11_MODULE(trivis_clipper, m) {
    m.doc() = "Python bindings for the Clipper2 polygon clipping library.";

    // ── PointD ────────────────────────────────────────────────────────────────
    py::class_<PointD>(m, "PointD")
        .def(py::init<double, double>(), py::arg("x"), py::arg("y"))
        .def_readwrite("x", &PointD::x)
        .def_readwrite("y", &PointD::y)
        .def("__repr__", [](const PointD& p) {
            return "PointD(" + std::to_string(p.x) + ", " + std::to_string(p.y) + ")";
        })
        .def("__iter__", [](const PointD& p) {
            return py::iter(py::make_tuple(p.x, p.y));
        })
        .def("__eq__", [](const PointD& a, const PointD& b) {
            return a.x == b.x && a.y == b.y;
        });

    // ── Enums ─────────────────────────────────────────────────────────────────
    py::enum_<FillRule>(m, "FillRule")
        .value("EvenOdd",  FillRule::EvenOdd)
        .value("NonZero",  FillRule::NonZero)
        .value("Positive", FillRule::Positive)
        .value("Negative", FillRule::Negative)
        .export_values();

    py::enum_<JoinType>(m, "JoinType")
        .value("Square", JoinType::Square)
        .value("Round",  JoinType::Round)
        .value("Miter",  JoinType::Miter)
        .export_values();

    py::enum_<EndType>(m, "EndType")
        .value("Polygon", EndType::Polygon)
        .value("Joined",  EndType::Joined)
        .value("Butt",    EndType::Butt)
        .value("Square",  EndType::Square)
        .value("Round",   EndType::Round)
        .export_values();

    // ── Boolean operations ────────────────────────────────────────────────────
    m.def("union_paths",
        [](const PathsD& subjects, const PathsD& clips, FillRule fill_rule, int precision) {
            return Union(subjects, clips, fill_rule, precision);
        },
        py::arg("subjects"), py::arg("clips"),
        py::arg("fill_rule") = FillRule::NonZero, py::arg("precision") = 2,
        "Compute the union of two sets of polygons.");

    m.def("intersect_paths",
        [](const PathsD& subjects, const PathsD& clips, FillRule fill_rule, int precision) {
            return Intersect(subjects, clips, fill_rule, precision);
        },
        py::arg("subjects"), py::arg("clips"),
        py::arg("fill_rule") = FillRule::NonZero, py::arg("precision") = 2,
        "Compute the intersection of two sets of polygons.");

    m.def("difference_paths",
        [](const PathsD& subjects, const PathsD& clips, FillRule fill_rule, int precision) {
            return Difference(subjects, clips, fill_rule, precision);
        },
        py::arg("subjects"), py::arg("clips"),
        py::arg("fill_rule") = FillRule::NonZero, py::arg("precision") = 2,
        "Subtract clip polygons from subject polygons.");

    m.def("xor_paths",
        [](const PathsD& subjects, const PathsD& clips, FillRule fill_rule, int precision) {
            return Xor(subjects, clips, fill_rule, precision);
        },
        py::arg("subjects"), py::arg("clips"),
        py::arg("fill_rule") = FillRule::NonZero, py::arg("precision") = 2,
        "Compute the symmetric difference of two sets of polygons.");

    // ── Inflate (offset/buffer) ───────────────────────────────────────────────
    m.def("inflate_paths",
        [](const PathsD& paths, double delta,
           JoinType join_type, EndType end_type,
           double miter_limit, int precision, double arc_tolerance) {
            return InflatePaths(paths, delta, join_type, end_type,
                                miter_limit, precision, arc_tolerance);
        },
        py::arg("paths"), py::arg("delta"),
        py::arg("join_type")      = JoinType::Round,
        py::arg("end_type")       = EndType::Polygon,
        py::arg("miter_limit")    = 2.0,
        py::arg("precision")      = 2,
        py::arg("arc_tolerance")  = 0.0,
        "Inflate (positive delta) or deflate (negative delta) a set of polygons.");

    // ── Simplify (Ramer–Douglas–Peucker) ─────────────────────────────────────
    m.def("simplify_paths",
        [](const PathsD& paths, double epsilon, bool closed) {
            return SimplifyPaths(paths, epsilon, closed);
        },
        py::arg("paths"), py::arg("epsilon"), py::arg("closed") = true,
        "Simplify paths using the Ramer–Douglas–Peucker algorithm.");
}
