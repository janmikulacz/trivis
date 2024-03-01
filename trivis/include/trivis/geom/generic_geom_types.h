/**
 * File:    generic_geom_types.h
 *
 * Date:    18.08.2020
 * Author:  Jan Mikula
 * E-mail:  jan.mikula@cvut.cz
 *
 */

#ifndef TRIVIS_GEOM_GENERIC_GEOM_TYPES_H_
#define TRIVIS_GEOM_GENERIC_GEOM_TYPES_H_

#include <cmath>
#include <vector>
#include <string>
#include <ostream>

namespace trivis::geom {

template<typename T>
struct Point {

    T x = T(0);
    T y = T(0);

    Point() = default;

    Point(T x, T y) : x(x), y(y) {}

    Point(const Point &) = default;

    Point(Point &&) noexcept = default;

    Point &operator=(const Point &) = default;

    Point &operator=(Point &&) noexcept = default;

    virtual ~Point() = default;

    Point operator+(const Point &p) const { return Point(x + p.x, y + p.y); }

    Point operator-(const Point &p) const { return Point(x - p.x, y - p.y); }

    Point operator-() const { return Point(-x, -y); }

    Point operator*(T f) const { return Point(f * x, f * y); }

    Point operator/(T f) const { return Point(x / f, y / f); }

    bool operator==(const Point &p) const { return x == p.x && y == p.y; }

    bool operator!=(const Point &p) const { return !(x == p.x && y == p.y); }

    T NormSquared() const { return x * x + y * y; }

    template<typename TT = T>
    TT Norm() const { return std::sqrt(static_cast<TT>(NormSquared())); }

    T SquaredDistanceTo(Point p) const { return (*this - p).NormSquared(); }

    template<typename TT = T>
    TT DistanceTo(Point p) const { return (*this - p).template Norm<TT>(); }

    template<typename TT = T>
    Point<TT> Copy() const { return Point(static_cast<TT>(x), static_cast<TT>(y)); }

    template<typename TT = T>
    Point<TT> CopyNormalized() const { return Point(static_cast<TT>(x), static_cast<TT>(y)) / Norm<TT>(); }

    template<typename TT = T>
    Point<TT> CopySwappedXY() const { return Point(static_cast<TT>(y), static_cast<TT>(x)); }

    void SwapXY() { std::swap(x, y); }

    [[nodiscard]] std::string ToString(const std::string &left_str = "{", const std::string &delimiter = ", ", const std::string &right_str = "}") const {
        return left_str + std::to_string(x) + delimiter + std::to_string(y) + right_str;
    }

    friend std::ostream &operator<<(std::ostream &os, const Point &p) {
        return os << "{" << p.x << ", " << p.y << "}";
    }

};

template<typename T>
inline auto MakePoint(const T &x, const T &y) { return Point<T>(x, y); }

template<typename T>
using Points = std::vector<Point<T>>;

template<typename T>
using Polygon = Points<T>;

template<typename T>
using Polygons = std::vector<Points<T>>;

template<typename T>
struct Position {

    T x = T(0);
    T y = T(0);
    T a = T(0);

    Position() = default;

    Position(const T &x, const T &y) : x(x), y(y) {}

    Position(const T &x, const T &y, const T &a) : x(x), y(y), a(a) {}

    Position(const Position &) = default;

    Position(Position &&) noexcept = default;

    Position &operator=(const Position &) = default;

    Position &operator=(Position &&) noexcept = default;

    virtual ~Position() = default;

    Position operator+(const Position &p) const { return Position(x + p.x, y + p.y, a + p.a); }

    Position operator-(const Position &p) const { return Position(x - p.x, y - p.y, a + p.a); }

    bool operator==(const Position &p) const { return x == p.x && y == p.y && a == p.a; }

    bool operator!=(const Position &p) const { return !(x == p.x && y == p.y && a == p.a); }

    Point<T> ToPoint() const { return Point<T>(x, y); }

    [[nodiscard]] std::string ToString(const std::string &left_str = "{", const std::string &delimiter = ", ", const std::string &right_str = "}") const {
        return left_str + std::to_string(x) + delimiter + std::to_string(y) + delimiter + std::to_string(a) + right_str;
    }

    friend std::ostream &operator<<(std::ostream &os, const Position &p) {
        return os << "{" << std::to_string(p.x) << ", " << std::to_string(p.y) << ", " << std::to_string(p.a) << "}";
    }

};

template<typename T>
inline auto MakePosition(const T &x, const T &y) { return Position<T>(x, y); }

template<typename T>
inline auto MakePosition(const T &x, const T &y, const T &a) { return Position<T>(x, y, a); }

template<typename T>
using Positions = std::vector<Position<T>>;

}

#endif //TRIVIS_GEOM_GENERIC_GEOM_TYPES_H_
