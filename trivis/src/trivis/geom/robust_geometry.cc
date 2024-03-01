/**
 * File:   robust_geometry.cc
 *
 * Date:   23.03.2022
 * Author: Jan Mikula
 * E-mail: jan.mikula@cvut.cz
 *
 */

#include "trivis/geom/robust_geometry.h"

#ifndef TRIVIS_USE_INEXACT_PREDICATES
#include "robust-predicate/geompred.hpp"
#endif

#include <cassert>

using namespace trivis;
using namespace trivis::geom;

Orientation geom::Orient(const FPoint &a, const FPoint &b, const FPoint &c) {
#ifdef TRIVIS_USE_INEXACT_PREDICATES
    double acx = a.x - c.x;
    double bcx = b.x - c.x;
    double acy = a.y - c.y;
    double bcy = b.y - c.y;
    double result = acx * bcy - acy * bcx;
#else
    auto result = geompred::orient2d(&a.x, &b.x, &c.x);
#endif
    if (result > 0.0) {
        return Orientation::kLeftTurn;
    }
    if (result < 0.0) {
        return Orientation::kRightTurn;
    }
    return Orientation::kCollinear;
}

Orientation geom::OrientWithEps(const FPoint &a, const FPoint &b, const FPoint &c, double eps) {
#ifdef TRIVIS_USE_INEXACT_PREDICATES
    double acx = a.x - c.x;
    double bcx = b.x - c.x;
    double acy = a.y - c.y;
    double bcy = b.y - c.y;
    double result = acx * bcy - acy * bcx;
#else
    auto result = geompred::orient2d(&a.x, &b.x, &c.x);
#endif
    if (result > eps) {
        return Orientation::kLeftTurn;
    }
    if (result < -eps) {
        return Orientation::kRightTurn;
    }
    return Orientation::kCollinear;
}

bool geom::TurnsLeft(const FPoint &a, const FPoint &b, const FPoint &c) {
    return Orient(a, b, c) == Orientation::kLeftTurn;
}

bool geom::TurnsRight(const FPoint &a, const FPoint &b, const FPoint &c) {
    return Orient(a, b, c) == Orientation::kRightTurn;
}

bool geom::Collinear(const FPoint &a, const FPoint &b, const FPoint &c) {
    return Orient(a, b, c) == Orientation::kCollinear;
}

bool geom::TurnsLeftWithEps(const FPoint &a, const FPoint &b, const FPoint &c, double eps) {
    return OrientWithEps(a, b, c, eps) == Orientation::kLeftTurn;
}

bool geom::TurnsRightWithEps(const FPoint &a, const FPoint &b, const FPoint &c, double eps) {
    return OrientWithEps(a, b, c, eps) == Orientation::kRightTurn;
}

bool geom::CollinearWithEps(const FPoint &a, const FPoint &b, const FPoint &c, double eps) {
    return OrientWithEps(a, b, c, eps) == Orientation::kCollinear;
}

char geom::PointTriangleRelation(const FPoint &q, const FPoint &a, const FPoint &b, const FPoint &c) {
    Orientation o_abq = Orient(a, b, q);
    if (o_abq == Orientation::kRightTurn) {
        return '0';
    }
    Orientation o_bcq = Orient(b, c, q);
    if (o_bcq == Orientation::kRightTurn) {
        return '0';
    }
    Orientation o_caq = Orient(c, a, q);
    if (o_caq == Orientation::kRightTurn) {
        return '0';
    }
    if (o_abq == Orientation::kCollinear) {
        if (o_bcq == Orientation::kCollinear) {
            return 'B';
        }
        if (o_caq == Orientation::kCollinear) {
            return 'A';
        }
        return 'a';
    }
    if (o_bcq == Orientation::kCollinear) {
        if (o_caq == Orientation::kCollinear) {
            return 'C';
        }
        return 'b';
    }
    if (o_caq == Orientation::kCollinear) {
        return 'c';
    }
    return '1';
}

char geom::PointTriangleRelationWithEps(const FPoint &q, const FPoint &a, const FPoint &b, const FPoint &c, double eps) {
    Orientation o_abq = OrientWithEps(a, b, q, eps);
    if (o_abq == Orientation::kRightTurn) {
        return '0';
    }
    Orientation o_bcq = OrientWithEps(b, c, q, eps);
    if (o_bcq == Orientation::kRightTurn) {
        return '0';
    }
    Orientation o_caq = OrientWithEps(c, a, q, eps);
    if (o_caq == Orientation::kRightTurn) {
        return '0';
    }
    if (o_abq == Orientation::kCollinear) {
        if (o_bcq == Orientation::kCollinear) {
            return 'B';
        }
        if (o_caq == Orientation::kCollinear) {
            return 'A';
        }
        return 'a';
    }
    if (o_bcq == Orientation::kCollinear) {
        if (o_caq == Orientation::kCollinear) {
            return 'C';
        }
        return 'b';
    }
    if (o_caq == Orientation::kCollinear) {
        return 'c';
    }
    return '1';
}

bool geom::IsPointInCone(const FPoint &q, const FPoint &a, const FPoint &b, const FPoint &c) {
    if (!TurnsRight(a, b, c)) {
        // convex at b
        return !TurnsRight(q, a, b) && !TurnsRight(q, b, c);
    }
    // reflex at b
    return !TurnsRight(q, a, b) || !TurnsRight(q, b, c);
}

FPolygon geom::RemoveDuplicatePoints(const FPolygon &polygon) {
    geom::FPolygon ret;
    int n = static_cast<int>(polygon.size());
    int i = 0;
    ret.reserve(n);
    while (true) {
        ret.push_back(polygon[i]);
        int j = i;
        bool same;
        while (j < n && (same = (polygon[++j] == polygon[i])));
        if (j == n - 1) {
            if (!same) {
                i = j;
            }
            if (polygon[i] != polygon[0]) {
                ret.push_back(polygon[i]);
            }
            break;
        }
        i = j;
    }
    return ret;
}

FPolygon geom::RemoveCollinearPoints(const FPolygon &polygon) {
    geom::FPolygon ret;
    int n = static_cast<int>(polygon.size());
    int u = n - 1, v = 0, w = 1;
    while (v < n) {
        if (!Collinear(polygon[u], polygon[v], polygon[w])) {
            ret.push_back(polygon[v]);
        }
        ++u, ++v, ++w;
        if (u == n) u = 0;
        if (w == n) w = 0;
    }
    return ret;
}

FPoint geom::FindAnyPointInNonConvexPolygon(const FPolygon &polygon, bool clockwise) {
    assert(!polygon.empty());

    int n = static_cast<int>(polygon.size());

    // Deal with a point, segment, and triangle.
    if (n <= 3) {
        FPoint p = polygon[0];
        for (int i = 1; i < polygon.size(); ++i) {
            p = p + polygon[i];
        }
        p = p / n;
        return p;
    }

    // 1) Identify a convex vertex v
    int idx_u = n - 1, idx_v = 0, idx_w = 1;
    while (idx_v < n) {
        if (clockwise && TurnsRight(polygon[idx_u], polygon[idx_v], polygon[idx_w])) {
            break;
        }
        if (!clockwise && TurnsLeft(polygon[idx_u], polygon[idx_v], polygon[idx_w])) {
            break;
        }
        ++idx_u, ++idx_v, ++idx_w;
        if (idx_u == n) idx_u = 0;
        if (idx_w == n) idx_w = 0;
    }
    if (clockwise) {
        std::swap(idx_u, idx_w);
    }
    const auto &u = polygon[idx_u];
    const auto &v = polygon[idx_v];
    const auto &w = polygon[idx_w];

    // 2) For each other vertex p do
    int idx_q = -1;
    double min_sq_distance_to_v = std::numeric_limits<double>::max();
    for (int idx_p = 0; idx_p < polygon.size(); ++idx_p) {
        const auto &p = polygon[idx_p];
        // 2a) if p is inside triangle uvw,
        if (idx_p != idx_u && idx_p != idx_v && idx_p != idx_w && TurnsLeft(u, v, p) && TurnsLeft(v, w, p) && !TurnsRight(w, u, p)) {
            // compute squared distance to p.
            auto sq_distance_to_v = p.SquaredDistanceTo(v);
            // 2b) Save point p as point q if distance is a new min.
            if (sq_distance_to_v < min_sq_distance_to_v) {
                min_sq_distance_to_v = sq_distance_to_v;
                idx_q = idx_p;
            }
        }
    }

    if (idx_q == -1) {
        // 3) If no point is inside, return centroid of uvw.
        return geom::MakePoint((u.x + v.x + w.x) / 3.0, (u.y + v.y + w.y) / 3.0);
    } else {
        // 4) Else if some point inside, qv is internal: return its midpoint.
        return geom::MakePoint((polygon[idx_q].x + v.x) / 2.0, (polygon[idx_q].y + v.y) / 2.0);
    }
}

int geom::PointInPolygon(
    const FPoint &p,
    const FPolygon &polygon
) {

    // FIXME: this function is not tested!

    // based on: https://github.com/mikolalysenko/robust-point-in-polygon

    double x = p.x;
    double y = p.y;
    int n = static_cast<int>(polygon.size());
    int inside = 1;
    int lim = n;
    for (int i = 0, j = n - 1; i < lim; j = i++) {
        const auto &a = polygon[i];
        const auto &b = polygon[j];
        double yi = a.y;
        double yj = b.y;
        if (yj < yi) {
            if (yj < y && y < yi) {
                Orientation s = Orient(a, b, p);
                if (s == Orientation::kCollinear) {
                    return 0;
                } else {
                    inside ^= static_cast<int>(s == Orientation::kLeftTurn) | 0;
                }
            } else if (y == yi) {
                const auto &c = polygon[(i + 1) % n];
                double yk = c.y;
                if (yi < yk) {
                    Orientation s = Orient(a, b, p);
                    if (s == Orientation::kCollinear) {
                        return 0;
                    } else {
                        inside ^= static_cast<int>(s == Orientation::kLeftTurn) | 0;
                    }
                }
            }
        } else if (yi < yj) {
            if (yi < y && y < yj) {
                Orientation s = Orient(a, b, p);
                if (s == Orientation::kCollinear) {
                    return 0;
                } else {
                    inside ^= static_cast<int>(s == Orientation::kLeftTurn) | 0;
                }
            } else if (y == yi) {
                const auto &c = polygon[(i + 1) % n];
                double yk = c.y;
                if (yk < yi) {
                    Orientation s = Orient(a, b, p);
                    if (s == Orientation::kCollinear) {
                        return 0;
                    } else {
                        inside ^= static_cast<int>(s == Orientation::kLeftTurn) | 0;
                    }
                }
            }
        } else if (y == yi) {
            double x0 = std::min(a.x, b.x);
            double x1 = std::max(a.x, b.x);
            if (i == 0) {
                while (j > 0) {
                    int k = (j + n - 1) % n;
                    const auto &pp = polygon[k];
                    if (pp.y != y) {
                        break;
                    }
                    double px = pp.x;
                    x0 = std::min(x0, px);
                    x1 = std::max(x1, px);
                    j = k;
                }
                if (j == 0) {
                    if (x0 <= x && x <= x1) {
                        return 0;
                    }
                    return 1;
                }
                lim = j + 1;
            }
            double y0 = polygon[(j + n - 1) % n].y;
            while (i + 1 < lim) {
                const auto &pp = polygon[i + 1];
                if (pp.y != y) {
                    break;
                }
                double px = pp.x;
                x0 = std::min(x0, px);
                x1 = std::max(x1, px);
                i += 1;
            }
            if (x0 <= x && x <= x1) {
                return 0;
            }
            double y1 = polygon[(i + 1) % n].y;
            if (x < x0 && (y0 < y != y1 < y)) {
                inside ^= 1;
            }
        }
    }
    return 2 * inside - 1;
}
