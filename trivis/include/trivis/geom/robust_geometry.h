/**
 * File:   robust_geometry.h
 *
 * Date:   22.03.2022
 * Author: Jan Mikula
 * E-mail: jan.mikula@cvut.cz
 *
 */

#ifndef TRIVIS_GEOM_ROBUST_GEOMETRY_H_
#define TRIVIS_GEOM_ROBUST_GEOMETRY_H_

#include "trivis/geom/geom_types.h"

namespace trivis::geom {

enum class Orientation : int {
    kRightTurn = -1,
    kCollinear = 0,
    kLeftTurn = 1
};

Orientation Orient(const FPoint &a, const FPoint &b, const FPoint &c);

Orientation OrientWithEps(const FPoint &a, const FPoint &b, const FPoint &c, double eps);

bool TurnsLeft(const FPoint &a, const FPoint &b, const FPoint &c);

bool TurnsRight(const FPoint &a, const FPoint &b, const FPoint &c);

bool Collinear(const FPoint &a, const FPoint &b, const FPoint &c);

bool TurnsLeftWithEps(const FPoint &a, const FPoint &b, const FPoint &c, double eps);

bool TurnsRightWithEps(const FPoint &a, const FPoint &b, const FPoint &c, double eps);

bool CollinearWithEps(const FPoint &a, const FPoint &b, const FPoint &c, double eps);

/**
 *
 * Returns true if c is on the extension of the ray starting in b in the direction b-->a, i.e., if (b-a)*(c-b) >= 0, and false otherwise.
 *  <---------a-------b
 *      FALSE | FALSE | TRUE
 * Exact implementation, works only if a, b, c are EXACTLY collinear.
 */
inline bool Extends(const FPoint &a, const FPoint &b, const FPoint &c) {
    return (a.x == b.x) ? ((a.y <= b.y) ? (b.y <= c.y) : (b.y >= c.y)) : ((a.x <= b.x) ? (b.x <= c.x) : (b.x >= c.x));
}

/**
 *  <---------a------b--------->
 *      FALSE | TRUE | FALSE
 *
 * Exact implementation, works only if a, b, c are EXACTLY collinear.
 */
inline bool IsBetween(const FPoint &a, const FPoint &b, const FPoint &c) {
    return !Extends(a, b, c) && !Extends(b, a, c);
}

/**
 *
 * @param q IN: The query point.
 * @param a IN: 1st point defining the counter-clockwise triangle a-b-c.
 * @param b IN: 2nd point defining the counter-clockwise triangle a-b-c.
 * @param c IN: 3rd point defining the counter-clockwise triangle a-b-c.
 * @return Code (char), see bellow:
 *      Code 'a': Point q lies on segment a-b.
 *      Code 'b': Point q lies on segment b-c.
 *      Code 'c': Point q lies on segment c-a.
 *      Code 'A': Point q is equal to point a.
 *      Code 'B': Point q is equal to point b
 *      Code 'C': Point q is equal to point c.
 *      Code '1': Point q is inside triangle a-b-c but not on its border.
 *      Code '0': Point is outside of triangle a-b-c.
 *  Point q is strictly inside OR on the border of a-b-c: '1', 'a', 'b', 'c', 'A', 'B', 'C'
 *  Point q lies on the border of a-b-c: 'a', 'b', 'c', 'A', 'B', 'C'.
 *  Point q is equal to one of the triangle vertices:  'A', 'B', 'C'.
 */
char PointTriangleRelation(const FPoint &q, const FPoint &a, const FPoint &b, const FPoint &c);

char PointTriangleRelationWithEps(const FPoint &q, const FPoint &a, const FPoint &b, const FPoint &c, double eps);

// FIXME: this function is not tested!
int PointInPolygon(const FPoint &p, const FPolygon &polygon);

bool IsPointInCone(const FPoint &q, const FPoint &a, const FPoint &b, const FPoint &c);

FPolygon RemoveDuplicatePoints(const FPolygon &polygon);

FPolygon RemoveCollinearPoints(const FPolygon &polygon);

FPoint FindAnyPointInNonConvexPolygon(const FPolygon &polygon, bool clockwise);

}

#endif //TRIVIS_GEOM_ROBUST_GEOMETRY_H_
