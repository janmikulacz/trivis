/**
 * File:    ptree_json.h
 *
 * Date:    17.11.2020
 * Author:  Jan Mikula
 * E-mail:  jan.mikula@cvut.cz
 *
 */

#ifndef TRIVIS_PLUS_UTILS_PTREE_JSON_H_
#define TRIVIS_PLUS_UTILS_PTREE_JSON_H_

#include <boost/property_tree/ptree.hpp>

#include "trivis/geom/geom_types.h"

#include <iostream>

namespace trivis_plus::utils {

using PTree = boost::property_tree::ptree;

class PTreeJSON {

 public:

  PTreeJSON() = default;

  explicit PTreeJSON(const PTree &pt) : pt_(pt) {}

  virtual ~PTreeJSON() = default;

  PTreeJSON(const PTreeJSON &other) = default;

  PTreeJSON(PTreeJSON &&other) = default;

  PTreeJSON &operator=(const PTreeJSON &other) = default;

  PTreeJSON &operator=(PTreeJSON &&other) = default;

  template<typename T>
  [[nodiscard]] auto Get(const std::string &path) const noexcept(false) {
    return pt_.get<T>(path);
  }

  template<typename T>
  void Get(const std::string &path, T &var) const noexcept(false) {
    var = pt_.get<T>(path);
  }

  [[nodiscard]] auto GetChild(const std::string &path) const noexcept(false) {
    return PTreeJSON(pt_.get_child(path));
  }

  template<typename T>
  [[nodiscard]] std::vector<T> GetVector(const std::string &path) const noexcept(false) {
    return PTreeToVec<T>(GetChild(path).pt_);
  }

  [[nodiscard]] trivis::geom::FPoint GetPoint(const std::string &path) const noexcept(false) {
    return GetChild(path).GetPoint();
  }

  [[nodiscard]] trivis::geom::FPoints GetPoints(const std::string &path) const noexcept(false) {
    return GetChild(path).GetPoints();
  }

  [[nodiscard]] trivis::geom::FPosition GetPosition(const std::string &path) const noexcept(false) {
    return GetChild(path).GetPosition();
  }

  [[nodiscard]] trivis::geom::Positions<double> GetPositions(const std::string &path) const noexcept(false) {
    return GetChild(path).GetPositions();
  }

  template<typename T>
  void PushBack(const T &val) {
    pt_.push_back(val);
  }

  template<typename T>
  void PushBack(const std::vector<T> &vec) {
    push_back(std::make_pair("", VecToPTree(vec)));
  }

  template<typename T>
  void Put(const std::string &path, const T &val) {
    pt_.put(path, val);
  }

  template<typename T>
  void PutChild(const std::string &path, const T &child) {
    pt_.put_child(path, child);
  }

  void PutChild(const std::string &path, const PTreeJSON &child) {
    pt_.put_child(path, child.pt_);
  }

  template<typename T>
  void PutGeom(const std::string &path, const T &val) {
    PTreeJSON pt_temp;
    pt_temp.PushBackGeom(val);
    PutChild(path, pt_temp.pt_);
  }

  template<typename T>
  void PutVector(const std::string &path, const std::vector<T> &vec) {
    PutChild(path, VecToPTree(vec));
  }

  template<typename T>
  void PutMatrix(const std::string &path, const std::vector<T> &mat) {
    PutChild(path, MatToPTree(mat));
  }

  void Read(const std::string &file_name) noexcept(false);

  void Write(const std::string &file_name) const noexcept(false);

 private:

  PTree pt_;

  [[nodiscard]] trivis::geom::FPoint GetPoint() const noexcept(false);

  [[nodiscard]] trivis::geom::FPoints GetPoints() const noexcept(false);

  [[nodiscard]] trivis::geom::FPoints GetPolygons() const noexcept(false);

  [[nodiscard]] trivis::geom::FPosition GetPosition() const noexcept(false);

  [[nodiscard]] trivis::geom::Positions<double> GetPositions() const noexcept(false);

  void PushBackGeom(const trivis::geom::FPoint &point);

  void PushBackGeom(const trivis::geom::FPoints &points);

  void PushBackGeom(const trivis::geom::Polygons<double> &polygons);

  void PushBackGeom(const trivis::geom::FPosition &position);

  void PushBackGeom(const trivis::geom::Positions<double> &positions);

  template<typename T>
  static PTree VecToPTree(const std::vector<T> &vec) {
    PTree pt;
    for (const auto &val : vec) {
      PTree child;
      child.put("", val);
      pt.push_back(std::make_pair("", child));
    }
    return pt;
  }

  template<typename T>
  static PTree MatToPTree(const std::vector<std::vector<T>> &mat) {
    PTree pt;
    for (const auto &vec : mat) {
      pt.push_back(std::make_pair("", VecToPTree(vec)));
    }
    return pt;
  }

  template<typename T>
  static std::vector<T> PTreeToVec(const PTree &pt) {
    std::vector<T> vec;
    vec.reserve(pt.size());
    for (const auto &val : pt) {
      T val_t;
      std::stringstream(val.second.data()) >> val_t;
      vec.push_back(val_t);
    }
    return vec;
  }

};

}

#endif //TRIVIS_PLUS_UTILS_PTREE_JSON_H_
