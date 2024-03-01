/**
 * File:    ptree_json.cc
 *
 * Date:    17.11.2020
 * Author:  Jan Mikula
 * E-mail:  jan.mikula@cvut.cz
 *
 */

#include "trivis_plus/utils/ptree_json.h"

#include <boost/property_tree/json_parser.hpp>

using namespace trivis;
using namespace trivis::geom;

using namespace trivis_plus;
using namespace trivis_plus::utils;

void utils::PTreeJSON::Read(const std::string &file_name) noexcept(false) {
  boost::property_tree::read_json(file_name, pt_);
}

void utils::PTreeJSON::Write(const std::string &file_name) const noexcept(false) {
  boost::property_tree::write_json(file_name, pt_);
}

FPoint utils::PTreeJSON::GetPoint() const noexcept(false) {
  return geom::MakePoint(Get<double>("x"), Get<double>("y"));
}

FPoints utils::PTreeJSON::GetPoints() const noexcept(false) {
  auto n = Get<int>("n");
  auto x = GetVector<double>("x");
  auto y = GetVector<double>("y");
  FPoints points;
  points.reserve(n);
  for (int kI = 0; kI < n; ++kI) {
    points.emplace_back(x[kI], y[kI]);
  }
  return points;
}

FPoints utils::PTreeJSON::GetPolygons() const noexcept(false) {
  // TODO
  std::cerr << "utils::PropertyTreeJson::GetPolygons() not yet implemented!\n";
  return FPoints();
}

FPosition utils::PTreeJSON::GetPosition() const noexcept(false) {
  return geom::MakePosition(Get<double>("x"), Get<double>("y"), Get<double>("a"));
}

FPositions utils::PTreeJSON::GetPositions() const noexcept(false) {
  auto n = Get<int>("n");
  auto x = GetVector<double>("x");
  auto y = GetVector<double>("y");
  auto a = GetVector<double>("a");
  FPositions positions;
  positions.reserve(n);
  for (int kI = 0; kI < n; ++kI) {
    positions.emplace_back(x[kI], y[kI], a[kI]);
  }
  return positions;
}

void utils::PTreeJSON::PushBackGeom(const FPoint &point) {
  Put("x", point.x);
  Put("y", point.y);
}

void utils::PTreeJSON::PushBackGeom(const FPoints &points) {
  std::vector<double> x(points.size());
  std::vector<double> y(points.size());
  for (int kI = 0; kI < points.size(); ++kI) {
    x[kI] = points[kI].x;
    y[kI] = points[kI].y;
  }
  auto pt_x = VecToPTree(x);
  auto pt_y = VecToPTree(y);
  Put("n", points.size());
  PutChild("x", pt_x);
  PutChild("y", pt_y);
}

void utils::PTreeJSON::PushBackGeom(const FPolygons &polygons) {
  Put("n", polygons.size());
  PTreeJSON pt_polygons;
  for (const auto &polygon : polygons) {
    PTreeJSON pt_polygon;
    pt_polygon.PushBackGeom(polygon);
    pt_polygons.PushBack(std::make_pair("", pt_polygon.pt_));
  }
  PutChild("polygons", pt_polygons);
}

void utils::PTreeJSON::PushBackGeom(const FPosition &position) {
  Put("x", position.x);
  Put("y", position.y);
  Put("a", position.a);
}

void utils::PTreeJSON::PushBackGeom(const FPositions &positions) {
  std::vector<double> x(positions.size());
  std::vector<double> y(positions.size());
  std::vector<double> a(positions.size());
  for (int kI = 0; kI < positions.size(); ++kI) {
    x[kI] = positions[kI].x;
    y[kI] = positions[kI].y;
    a[kI] = positions[kI].a;
  }
  auto pt_x = VecToPTree(x);
  auto pt_y = VecToPTree(y);
  auto pt_a = VecToPTree(a);
  Put("n", positions.size());
  PutChild("x", pt_x);
  PutChild("y", pt_y);
  PutChild("a", pt_a);
}
