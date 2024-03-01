/**
 * File:   intersections.h
 *
 * Date:   14.04.2022
 * Author: Jan Mikula
 * E-mail: jan.mikula@cvut.cz
 *
 */

#ifndef TRIVIS_GEOM_INTERSECTIONS_H_
#define TRIVIS_GEOM_INTERSECTIONS_H_

#include "trivis/geom/geom_types.h"

namespace trivis::geom {

bool LineLineIntersectionNotCollinear(const FPoint &a, const FPoint &b, const FPoint &c, const FPoint &d, FPoint &p);

/**
 *
 * @param a IN: 1st endpoint of segment a-b.
 * @param b IN: 2nd endpoint of segment a-b.
 * @param c IN: 1st endpoint of segment c-d.
 * @param d IN: 2nd endpoint of segment c-d.
 * @param p OUT: 1st endpoint of the intersection segment, if a-b and c-d intersect in a segment (see Code 's'); else intersection of a-b and c-d if it exists.
 * @param q OUT: 2nd endpoint of the intersection segment, if a-b and c-d intersect in a segment (see Code 's').
 * @return Code (char), see bellow:
 *      Code 'c': Segments are collinear and do not intersect.
 *      Code 'V': Segments are collinear and share a single endpoint p.
 *      Code 's': Segments are collinear and intersect in segment p-q.
 *      Code '0': Segments are NOT collinear and do not intersect.
 *      Code 'v': Segments are NOT collinear and share a single endpoint p.
 *      Code 'i': Segments are NOT collinear and intersect in p, which is an endpoint of either a-b or c-d.
 *      Code '1': Segments are NOT collinear and have a PROPER intersection in p.
 *      Code 'e': Segments are NOT collinear and SHOULD intersect but the intersection was not found (ERROR STATE).
 *  Segments are collinear: 'c', 'V', 's'.
 *  Segments are NOT collinear: '0', 'v', 'i', '1', 'e'.
 *  Segments intersect (intersection found): 'V', 's', 'v', 'i', '1'.
 *  Segments intersect (intersection NOT found): 'e'.
 *  Segments DO NOT intersect: 'c', '0'.
 *  Intersection NOT found: 'c', '0', 'e'.
 */
char SegmentSegmentIntersection(const FPoint &a, const FPoint &b, const FPoint &c, const FPoint &d, FPoint &p, FPoint &q);

/**
 *
 * @param a IN: An endpoint of ray a-b->.
 * @param b IN: A point defining the direction of a-b->.
 * @param c IN: 1st endpoint of segment c-d.
 * @param d IN: 2nd endpoint of segment c-d.
 * @param p OUT: 1st endpoint of the intersection segment, if a-b-> and c-d intersect in a segment (see Code 's'); else intersection of a-b-> and c-d if it exists.
 * @param q OUT: 2nd endpoint of the intersection segment, if a-b-> and c-d intersect in a segment (see Code 's').
 * @return Code (char), see bellow:
 *      Code 'c': The ray and segment are collinear and do not intersect.
 *      Code 'V': The ray and segment are collinear and share a single endpoint p.
 *      Code 's': The ray and segment are collinear and intersect in segment p-q.
 *      Code '0': The ray and segment are NOT collinear and do not intersect.
 *      Code 'v': The ray and segment are NOT collinear and share a single endpoint p.
 *      Code 'i': The ray and segment are NOT collinear and intersect in p, which is either a or an endpoint of c-d.
 *      Code '1': The ray and segment are NOT collinear and have a PROPER intersection in p.
 *      Code 'e': The ray and segment are NOT collinear and SHOULD intersect but the intersection was not found (ERROR STATE).
 *  The ray and segment are collinear: 'c', 'V', 's'.
 *  The ray and segment are NOT collinear: '0', 'v', 'i', '1', 'e'.
 *  The ray and segment intersect (intersection found): 'V', 's', 'v', 'i', '1'.
 *  The ray and segment intersect (intersection NOT found): 'e'.
 *  The ray and segment DO NOT intersect: 'c', '0'.
 *  Intersection NOT found: 'c', '0', 'e'.
 */
char RaySegmentIntersection(const FPoint &a, const FPoint &b, const FPoint &c, const FPoint &d, FPoint &p, FPoint &q);

geom::Points<double> LineCircleIntersections(const FPoint &p1, const FPoint &p2, const FPoint &cp, double r, bool segment, double eps = 1e-6);

char PointSegmentClosenessRelation(const FPoint &q, const FPoint &a, const FPoint &b);

double PointSegmentSquaredDistance(const FPoint &q, const FPoint &a, const FPoint &b, char *flag = nullptr);

}

#endif //TRIVIS_GEOM_INTERSECTIONS_H_
