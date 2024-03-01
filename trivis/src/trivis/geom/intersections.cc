/**
 * File:   intersections.cc
 *
 * Date:   14.04.2022
 * Author: Jan Mikula
 * E-mail: jan.mikula@cvut.cz
 *
 */

#include "trivis/geom/intersections.h"

#include "trivis/geom/robust_geometry.h"

using namespace trivis;
using namespace trivis::geom;

/**
 * Assumption: 'abc' is not collinear AND 'abd' is not collinear.
 */
bool geom::LineLineIntersectionNotCollinear(const FPoint &a, const FPoint &b, const FPoint &c, const FPoint &d, FPoint &p) {
    double den = a.x * (d.y - c.y) + b.x * (c.y - d.y) + d.x * (b.y - a.y) + c.x * (a.y - b.y);
    // den should never be zero if the assumption is satisfied
    if (den == 0.0) {
        // just to be sure (should not happen!)
        return false;
    }
    double num = a.x * (d.y - c.y) + c.x * (a.y - d.y) + d.x * (c.y - a.y);
    double t = num / den;
    p.x = a.x + t * (b.x - a.x);
    p.y = a.y + t * (b.y - a.y);
    return true;
}

/**
 * Assumption: 'abc' is collinear.
 */
inline char SegmentPointIntersectionCollinear(const FPoint &a, const FPoint &b, const FPoint &c, FPoint &p) {
    if (IsBetween(a, b, c)) {
        p = c;
        return 'i';
    }
    if (c == a || c == b) {
        p = c;
        return 'v';
    }
    return 'c';
}

/**
 * Assumption: 'c' is left from 'ab' AND 'd' is right from 'ab'.
 */
char SegmentSegmentIntersectionGeneral(const FPoint &a, const FPoint &b, const FPoint &c, const FPoint &d, FPoint &p) {
    Orientation o_cda = Orient(c, d, a);
    if (o_cda == Orientation::kLeftTurn) {
        return '0';
    }
    if (o_cda == Orientation::kRightTurn) {
        Orientation o_cdb = Orient(c, d, b);
        if (o_cdb == Orientation::kLeftTurn) {
            // The intersection must exist!
            if (LineLineIntersectionNotCollinear(a, b, c, d, p)) {
                return '1';
            }
            // Should not happen!
            return '0';
        }
        if (o_cdb == Orientation::kRightTurn) {
            // No intersection.
            return '0';
        }
        // The intersection must be in 'b'.
        p = b;
        return 'v';
    }
    // The intersection must be in 'a'.
    p = a;
    return 'v';
}

/**
 * Assumption: 'abc' is collinear AND 'abd' is collinear.
 */
char SegmentSegmentIntersectionCollinear(const FPoint &a, const FPoint &b, const FPoint &c, const FPoint &d, FPoint &p, FPoint &q) {
    if (Extends(a, b, c)) {
        if (Extends(a, b, d)) {
            // The intersection either is in 'b' or does not exist.
            if (b == c || b == d) {
                p = b;
                return 'V';
            }
            // No intersection.
            return 'c';
        }
        if (Extends(b, a, d)) {
            // The intersection is the whole segment 'ab'.
            p = a;
            q = b;
            return 's';
        }
        // The intersection is segment 'db'.
        p = d;
        q = b;
        return 's';
    }
    if (Extends(b, a, c)) {
        if (Extends(b, a, d)) {
            // The intersection either is in 'a' or does not exist.
            if (a == c || a == d) {
                p = a;
                return 'V';
            }
            // No intersection.
            return 'c';
        }
        if (Extends(a, b, d)) {
            // The intersection is the whole segment 'ab'.
            p = a;
            q = b;
            return 's';
        }
        // The intersection is segment 'ad'.
        p = a;
        q = d;
        return 's';
    }
    // Point 'c' is between 'a' and 'b'.
    if (Extends(a, b, d)) {
        // The intersection is segment 'cb'.
        p = c;
        q = b;
        return 's';
    }
    if (Extends(b, a, d)) {
        // The intersection is segment 'ac'.
        p = a;
        q = c;
        return 's';
    }
    // The intersection is the whole segment 'cd'.
    p = c;
    q = d;
    return 's';
}

char geom::SegmentSegmentIntersection(const FPoint &a, const FPoint &b, const FPoint &c, const FPoint &d, FPoint &p, FPoint &q) {
    Orientation o_abc = Orient(a, b, c);
    Orientation o_abd = Orient(a, b, d);
    if (o_abc == Orientation::kLeftTurn) {
        if (o_abd == Orientation::kLeftTurn) {
            // === LL ===
            // No intersection.
            return '0';
        }
        if (o_abd == Orientation::kRightTurn) {
            // === LR ===
            // The segments are in general position, the intersection can also be in either 'a' or 'b'.
            return SegmentSegmentIntersectionGeneral(a, b, c, d, p);
        }
        // === LC ===
        // The intersection either is in 'd' or does not exist.
        return SegmentPointIntersectionCollinear(a, b, d, p);
    }
    if (o_abc == Orientation::kRightTurn) {
        if (o_abd == Orientation::kLeftTurn) {
            // === RL ===
            // The segments are in general position, the intersection can also be in either 'a' or 'b'.
            return SegmentSegmentIntersectionGeneral(a, b, d, c, p);
        }
        if (o_abd == Orientation::kRightTurn) {
            // === RR ===
            // No intersection.
            return '0';
        }
        // === RC ===
        // The intersection either is in 'd' or does not exist.
        return SegmentPointIntersectionCollinear(a, b, d, p);
    }
    if (o_abd == Orientation::kLeftTurn) {
        // === CL ===
        // The intersection either is in 'c' or does not exist.
        return SegmentPointIntersectionCollinear(a, b, c, p);
    }
    if (o_abd == Orientation::kRightTurn) {
        // === CR ===
        // The intersection either is in 'c' or does not exist.
        return SegmentPointIntersectionCollinear(a, b, c, p);
    }
    // === CC ===
    // Segments are collinear.
    return SegmentSegmentIntersectionCollinear(a, b, c, d, p, q);
}

/**
 * Assumption: 'abc' is collinear.
 */
inline char RayPointIntersectionCollinear(const FPoint &a, const FPoint &b, const FPoint &c, FPoint &p) {
    if (Extends(b, a, c)) {
        if (c == a) {
            p = c;
            return 'v';
        }
        return 'c';
    }
    p = c;
    return 'i';
}

/**
 * Assumption: 'c' is left from 'ab' AND 'd' is right from 'ab'.
 */
char RaySegmentIntersectionGeneral(const FPoint &a, const FPoint &b, const FPoint &c, const FPoint &d, FPoint &p) {
    Orientation o_cda = Orient(c, d, a);
    if (o_cda == Orientation::kLeftTurn) {
        return '0';
    }
    if (o_cda == Orientation::kRightTurn) {
        // The intersection must exist!
        if (LineLineIntersectionNotCollinear(a, b, c, d, p)) {
            return '1';
        }
        // Should not happen!
        return '0';
    }
    // The intersection must be in 'a'.
    p = a;
    return 'v';
}

/**
 * Assumption: 'abc' is collinear AND 'abd' is collinear.
 */
char RaySegmentIntersectionCollinear(const FPoint &a, const FPoint &b, const FPoint &c, const FPoint &d, FPoint &p, FPoint &q) {
    if (Extends(b, a, c)) {
        if (Extends(b, a, d)) {
            // The intersection either is in 'a' or does not exist.
            if (a == c || a == d) {
                p = a;
                return 'V';
            }
            // No intersection.
            return 'c';
        }
        // The intersection is segment 'ad'.
        p = a;
        q = d;
        return 's';
    }
    // Point 'c' is right of 'a'.
    if (Extends(b, a, d)) {
        // The intersection is segment 'ac'.
        p = a;
        q = c;
        return 's';
    }
    // The intersection is the whole segment 'cd'.
    p = c;
    q = d;
    return 's';
}

char geom::RaySegmentIntersection(const FPoint &a, const FPoint &b, const FPoint &c, const FPoint &d, FPoint &p, FPoint &q) {
    Orientation o_abc = Orient(a, b, c);
    Orientation o_abd = Orient(a, b, d);
    if (o_abc == Orientation::kLeftTurn) {
        if (o_abd == Orientation::kLeftTurn) {
            // === LL ===
            // No intersection.
            return '0';
        }
        if (o_abd == Orientation::kRightTurn) {
            // === LR ===
            // The ray and segment are in general position, the intersection can also be in either 'a' or 'b'.
            return RaySegmentIntersectionGeneral(a, b, c, d, p);
        }
        // === LC ===
        // The intersection either is in 'd' or does not exist.
        return RayPointIntersectionCollinear(a, b, d, p);
    }
    if (o_abc == Orientation::kRightTurn) {
        if (o_abd == Orientation::kLeftTurn) {
            // === RL ===
            // The ray and segment are in general position, the intersection can also be in either 'a' or 'b'.
            return RaySegmentIntersectionGeneral(a, b, d, c, p);
        }
        if (o_abd == Orientation::kRightTurn) {
            // === RR ===
            // No intersection.
            return '0';
        }
        // === RC ===
        // The intersection either is in 'd' or does not exist.
        return RayPointIntersectionCollinear(a, b, d, p);
    }
    if (o_abd == Orientation::kLeftTurn) {
        // === CL ===
        // The intersection either is in 'c' or does not exist.
        return RayPointIntersectionCollinear(a, b, c, p);
    }
    if (o_abd == Orientation::kRightTurn) {
        // === CR ===
        // The intersection either is in 'c' or does not exist.
        return RayPointIntersectionCollinear(a, b, c, p);
    }
    // === CC ===
    // Segments are collinear.
    return RaySegmentIntersectionCollinear(a, b, c, d, p, q);
}

inline double sq(double x) {
    return x * x;
}

geom::Points<double> geom::LineCircleIntersections(const FPoint &p1, const FPoint &p2, const FPoint &cp, double r, bool segment, double eps) {
    geom::Points<double> res;
    auto x0 = cp.x;
    auto y0 = cp.y;
    auto x1 = p1.x;
    auto y1 = p1.y;
    auto x2 = p2.x;
    auto y2 = p2.y;
    auto A = y2 - y1;
    auto B = x1 - x2;
    auto C = x2 * y1 - x1 * y2;
    auto a = sq(A) + sq(B);
    double b, c;
    bool bnz = true;
    if (std::abs(B) >= eps) {
        b = 2 * (A * C + A * B * y0 - sq(B) * x0);
        c = sq(C) + 2 * B * C * y0 - sq(B) * (sq(r) - sq(x0) - sq(y0));
    } else {
        b = 2 * (B * C + A * B * x0 - sq(A) * y0);
        c = sq(C) + 2 * A * C * x0 - sq(A) * (sq(r) - sq(x0) - sq(y0));
        bnz = false;
    }
    auto d = sq(b) - 4 * a * c; // discriminant
    if (d < 0) {
        return res;
    }

    // checks whether a point is within a segment
    auto within = [x1, y1, x2, y2, eps](double x, double y) {
        auto d1 = std::sqrt(sq(x2 - x1) + sq(y2 - y1));  // distance between end-points
        auto d2 = std::sqrt(sq(x - x1) + sq(y - y1));    // distance from point to one end
        auto d3 = std::sqrt(sq(x2 - x) + sq(y2 - y));    // distance from point to other end
        auto delta = d1 - d2 - d3;
        return std::abs(delta) < eps;                    // true if delta is less than a small tolerance
    };

    auto fx = [A, B, C](double x) {
        return -(A * x + C) / B;
    };

    auto fy = [A, B, C](double y) {
        return -(B * y + C) / A;
    };

    auto rxy = [segment, &res, within](double x, double y) {
        if (!segment || within(x, y)) {
            res.emplace_back(x, y);
        }
    };

    double x, y;
    if (d == 0.0) {
        // line is tangent to circle, so just one intersect at most
        if (bnz) {
            x = -b / (2 * a);
            y = fx(x);
            rxy(x, y);
        } else {
            y = -b / (2 * a);
            x = fy(y);
            rxy(x, y);
        }
    } else {
        // two intersects at most
        d = std::sqrt(d);
        if (bnz) {
            x = (-b + d) / (2 * a);
            y = fx(x);
            rxy(x, y);
            x = (-b - d) / (2 * a);
            y = fx(x);
            rxy(x, y);
        } else {
            y = (-b + d) / (2 * a);
            x = fy(y);
            rxy(x, y);
            y = (-b - d) / (2 * a);
            x = fy(y);
            rxy(x, y);
        }
    }

    return res;
}

char geom::PointSegmentClosenessRelation(const FPoint &q, const FPoint &a, const FPoint &b) {
    const double r_num = (q.x - a.x) * (b.x - a.x) + (q.y - a.y) * (b.y - a.y);
    const double r_den = (b.x - a.x) * (b.x - a.x) + (b.y - a.y) * (b.y - a.y);
    const double r = r_num / r_den;
    if (0.0 <= r && r <= 1.0) {
        return 'l';
    }
    // segment distance is the distance from c to a or from c to b
    const double dist1 = q.SquaredDistanceTo(a);
    const double dist2 = q.SquaredDistanceTo(b);
    if (dist1 < dist2) {
        return 'a';
    } else {
        return 'b';
    }
}

double geom::PointSegmentSquaredDistance(const FPoint &q, const FPoint &a, const FPoint &b, char *flag) {
    const double r_num = (q.x - a.x) * (b.x - a.x) + (q.y - a.y) * (b.y - a.y);
    const double r_den = (b.x - a.x) * (b.x - a.x) + (b.y - a.y) * (b.y - a.y);
    const double r = r_num / r_den;
    if (0.0 <= r && r <= 1.0) {
        // segment distance is line distance
        const double s = (a.y - q.y) * (b.x - a.x) - (a.x - q.x) * (b.y - a.y);
        if (flag) *flag = 'l';
        return s * s / r_den;
    }
    // segment distance is the distance from c to a or from c to b
    const double dist1 = q.SquaredDistanceTo(a);
    const double dist2 = q.SquaredDistanceTo(b);
    if (dist1 < dist2) {
        if (flag) *flag = 'a';
        return dist1;
    } else {
        if (flag) *flag = 'b';
        return dist2;
    }
}
