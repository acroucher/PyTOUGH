"""Geometry calculations.

Copyright 2011 University of Auckland.

This file is part of PyTOUGH.

PyTOUGH is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

PyTOUGH is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with PyTOUGH.  If not, see <http://www.gnu.org/licenses/>."""

from __future__ import print_function

try:
    import numpy as np
    from numpy import float64
    from numpy import int8 as int8
    from numpy.linalg import norm, solve, LinAlgError
except ImportError: # try importing Numeric on old installs
    import Numeric as np
    from Numeric import Float64 as float64
    from Numeric import Int8 as int8
    def norm(a): # Numeric doesn't have a norm function
        from math import sqrt
        return sqrt(np.dot(a, a))
    from LinearAlgebra import solve_linear_equations as solve
    from LinearAlgebra import LinAlgError

# check Python version:
import sys
vinfo = sys.version_info
if (vinfo[0] <= 2) and (vinfo[1] < 5): # < 2.5:
    def any(a):
        for x in a:
            if x: return True
        return False
    def all(a):
        for x in a:
            if not x: return False
        return True

def in_polygon(pos, polygon):
    """Tests if the (2-D) point a lies within a given polygon."""
    tolerance = 1.e-6
    numcrossings = 0
    ref = polygon[0]
    v = pos - ref
    for i in range(len(polygon)):
        p1 = polygon[i] - ref
        i2 = (i+1) % len(polygon)
        p2 = polygon[i2] - ref
        if p1[1] <= v[1] < p2[1] or p2[1] <= v[1] < p1[1]:
            d = p2 - p1
            if abs(d[1]) > tolerance:
                x = p1[0] + (v[1] - p1[1]) * d[0] / d[1]
                if v[0] < x: numcrossings += 1
    return (numcrossings % 2)

def in_rectangle(pos, rect):
    """Tests if the 2-D point lies in an axis-aligned rectangle, defined
    as a two-element list of arrays [bottom left, top right].
    """
    return all([rect[0][i] <= pos[i] <= rect[1][i] for i in range(2)])

def rectangles_intersect(rect1, rect2):
    """Returns True if two rectangles intersect."""
    return all([(rect1[1][i] >= rect2[0][i]) and
                (rect2[1][i] >= rect1[0][i]) for i in range(2)])

def sub_rectangles(rect):
    """Returns the sub-rectangles formed by subdividing the given rectangle evenly in four."""
    centre = 0.5 * (rect[0] + rect[1])
    r0 = [rect[0], centre]
    r1 = [np.array([centre[0], rect[0][1]]), np.array([rect[1][0], centre[1]])]
    r2 = [np.array([rect[0][0], centre[1]]), np.array([centre[0], rect[1][1]])]
    r3 = [centre, rect[1]]
    return [r0, r1, r2, r3]

def bounds_of_points(points):
    """Returns bounding box around the specified 2D points."""
    bottomleft = np.array([min([pos[i] for pos in points]) for i in range(2)])
    topright = np.array([max([pos[i] for pos in points]) for i in range(2)])
    return [bottomleft, topright]

def rect_to_poly(rect):
    """Converts a rectangle to a polygon."""
    return [rect[0], np.array(rect[1][0], rect[0][1]),
            rect[1], np.array(rect[0][0], rect[1][1])]

def polygon_area(polygon):
    """Calculates the area of an arbitrary polygon."""
    area = 0.0
    n = len(polygon)
    if n > 0:
        polygon -= polygon[0]
        for j, p1 in enumerate(polygon):
            p2 = polygon[(j+1) % n]
            area += p1[0] * p2[1] - p2[0] * p1[1]
    return 0.5 * area

def polygon_centroid(polygon):
    """Calculates the centroid of an arbitrary polygon."""
    c, area = np.zeros(2), 0.
    n = len(polygon)
    shift = polygon[0]
    polygon -= shift # shift to reduce roundoff for large coordinates
    if n < 3: return sum(polygon) / n + shift
    else:
        for j, p1 in enumerate(polygon):
            p2 = polygon[(j+1) % n]
            t = p1[0] * p2[1] - p2[0] * p1[1]
            area += t
            c += (p1 + p2) * t
        area *= 0.5
        return c / (6. * area) + shift

def line_polygon_intersections(polygon, line, bound_line = (True,True),
                               indices = False):
    """Returns a list of the intersection points at which a line crosses a
    polygon.  The list is sorted by distance from the start of the
    line.  The parameter bound_line controls whether to limit
    intersections between the line's start and end points.  If indices
    is True, also return polygon side indices of intersections.
    """
    crossings = []
    ref = polygon[0]
    l1, l2 = line[0] - ref, line[1] - ref
    tol = 1.e-12
    ind = {}
    def in_unit(x): return -tol <= x <= 1.0 + tol
    for i, p in enumerate(polygon):
        p1 = p - ref
        p2 = polygon[(i+1)%len(polygon)] - ref
        dp = p2 - p1
        A, b = np.column_stack((dp, l1 - l2)), l1 - p1
        try:
            xi = solve(A,b)
            inline = True
            if bound_line[0]: inline = inline and (-tol <= xi[1])
            if bound_line[1]: inline = inline and (xi[1] <= 1.0 + tol)
            if in_unit(xi[0]) and inline :
                c = tuple(ref + p1 + xi[0] * dp)
                ind[c] = i
        except LinAlgError: continue
    crossings = [np.array(c) for c, i in ind.items()]
    # Sort by distance from start of line:
    sortindex = np.argsort([norm(c - line[0]) for c in crossings])
    if indices: return [crossings[i] for i in sortindex], \
       [ind[tuple(crossings[i])] for i in sortindex]
    else: return [crossings[i] for i in sortindex]

def polyline_polygon_intersections(polygon, polyline):
    """Returns a list of intersection points at which a polyline (list of
    2-D points) crosses a polygon."""
    N = len(polyline)
    intersections = [
        line_polygon_intersections(polygon, [pt, polyline[(i+1) % N]])
        for i, pt in enumerate(polyline)]
    from itertools import chain # flatten list of lists
    return list(chain.from_iterable(intersections))

def simplify_polygon(polygon, tolerance = 1.e-6):
    """Simplifies a polygon by deleting colinear points.  The tolerance
    for detecting colinearity of points can optionally be
    specified.
    """
    s = []
    N = len(polygon)
    for i, p in enumerate(polygon[:]):
        l, n = polygon[i-1], polygon[(i+1)%N]
        dl, dn = p - l, n - p
        if np.dot(dl, dn) < (1. - tolerance) * norm(dl) * norm(dn):
            s.append(p)
    return s

def polygon_boundary(this, other, polygon):
    """Returns point on a line between vector this and other and also on
    the boundary of the polygon"""
    big = 1.e10
    ref = polygon[0]
    a = this - ref
    b = other - ref
    v = None
    dmin = big
    for i in range(len(polygon)):
        c = polygon[i] - ref
        i2 = (i + 1) % len(polygon)
        d = polygon[i2] - ref
        M = np.transpose(np.array([b - a, c - d]))
        try:
            r = solve(M, c - a)
            if r[0] >= 0.0 and 0.0 <= r[1] <= 1.0:
                bdy = c * (1.0 - r[1]) + d * r[1]
                dist = norm(bdy - a)
                if dist < dmin:
                    dmin = dist
                    v = bdy + ref
        except LinAlgError: continue
    return v

def line_projection(a, line, return_xi = False):
    """Finds projection of point a onto a line (defined by two vectors).  Optionally
    return the non-dimensional distance xi between the line start and end."""
    d = line[1] - line[0]
    try:
        xi = np.dot(a - line[0], d) / np.dot(d, d)
        p = line[0] + d * xi
    except ZeroDivisionError: # line ill-defined
        p, xi = None, None
    if return_xi: return p, xi
    else: return p

def point_line_distance(a, line):
    """Finds distance between point a and a line."""
    return np.linalg.norm(a - line_projection(a, line))

def polyline_line_distance(polyline, line):
    """Returns minimum distance between a polyline and a line."""
    dists = []
    for i, pt in enumerate(polyline):
        pline = [pt, polyline[(i + 1) % len(polyline)]]
        if line_polygon_intersections(line, pline): return 0.0
        else: dists.append(min([point_line_distance(p, line) for p in pline]))
    return min(dists)

def vector_heading(p):
    """Returns heading angle of a 2-D vector p, in radians clockwise from
    the y-axis ('north')."""
    from math import asin
    theta = asin(p[0] / norm(p))
    if p[1] < 0: theta = np.pi - theta
    if theta < 0: theta += 2. * np.pi
    elif theta > 2. * np.pi: theta -= 2. * np.pi
    return theta

class linear_trans2(object):
    """Class for 2D linear transformation (Ax+b)."""
    def __init__(self, A = None, b = None):
        if A is None: A = np.identity(2)
        if b is None: b = np.zeros(2)
        self.A = A
        self.b = b
    def __repr__(self): return repr(self.A) + ' ' + repr(self.b)
    def __call__(self, x):
        """Transforms a 2-D point x, or another linear transformation."""
        if isinstance(x, linear_trans2):
            result = linear_trans2()
            result.A = np.dot(self.A, x.A)
            result.b = self.b + np.dot(self.A, x.b)
            return result
        else: return np.dot(self.A, x) + self.b
    def get_inverse(self):
        """Returns the inverse of a transformation, if it exists."""
        result = linear_trans2()
        try:
            result.A = np.linalg.inv(self.A)
            result.b= -np.dot(result.A, self.b)
            return result
        except np.linalg.LinAlgError: return None
    inverse = property(get_inverse)
    def between_rects(self, rect1, rect2):
        """Returns a linear transformation mapping one rectangle to another.
        Each rectangle should be a list or tuple of two 2-D numpy
        arrays (bottom left and top right corners).
        """
        tolerance = 1.e-6
        d1 = rect1[1] - rect1[0]
        if abs(d1[0]) > tolerance and abs(d1[1]) > tolerance:
            d2 = rect2[1] - rect2[0]
            s = [l2 / l1 for l1, l2 in zip(d1, d2)]
            self.A = np.array([[s[0], 0.0], [0.0, s[1]]])
            self.b = rect2[0] - rect1[0]
            self.b = np.array([r2 - s * r1 for r1, r2, s in zip(rect1[0], rect2[0], s)])
        else:
            self.A=None
        return self
    def between_points(self, points1, points2):
        """Returns a least-squares best fit linear transformation between two
        sets of points (e.g. the same set of points in two different
        coordinate systems).  The parameters points1 and points2
        should contain lists of the corresponding points in the first
        and second coordinate systems respectively.  At least three
        points are needed to find a best-fit linear transformation.
        """
        npts1, npts2 = len(points1), len(points2)
        if npts1 == npts2:
            npts = npts1
            if npts >= 3:
                # first normalise the coordinates:
                r1, r2 = bounds_of_points(points1), bounds_of_points(points2)
                dx1 = np.array([r1[1][0] - r1[0][0], r1[1][1] - r1[0][1]])
                dx2 = np.array([r2[1][0] - r2[0][0], r2[1][1] - r2[0][1]])
                N1 = linear_trans2(np.array([[1. / dx1[0], 0.], [0., 1. / dx1[1]]]),
                                   -np.array([r1[0][0] / dx1[0], r1[0][1] / dx1[1]]))
                N2 = linear_trans2(np.array([[1. / dx2[0], 0.], [0., 1. / dx2[1]]]),
                                   -np.array([r2[0][0] / dx2[0], r2[0][1] / dx2[1]]))
                n1, n2 = [N1(p) for p in points1], [N2(p) for p in points2]
                # assemble the least squares fitting system for the normalised points:
                M, r = np.zeros((6,6)), np.zeros(6)
                for p1, p2 in zip(n1, n2):
                    p00, p01, p11 = p1[0] * p1[0], p1[0] * p1[1], p1[1] * p1[1]
                    M[0,0] += p00; M[0,1] += p01; M[0,4] += p1[0]
                    M[1,0] += p01; M[1,1] += p11; M[1,4] += p1[1]
                    M[2,2] += p00; M[2,3] += p01; M[2,5] += p1[0]
                    M[3,2] += p01; M[3,3] += p11; M[3,5] += p1[1]
                    M[4,0] += p1[0]; M[4,1] += p1[1]; M[4,4] += 1.0
                    M[5,2] += p1[0]; M[5,3] += p1[1]; M[5,5] += 1.0
                    r += np.array([p1[0] * p2[0], p1[1] * p2[0], p1[0] * p2[1],
                                   p1[1] * p2[1], p2[0], p2[1]])
                try:
                    a = np.linalg.solve(M, r)
                    L = linear_trans2(np.array([[a[0], a[1]], [a[2], a[3]]]),
                                      np.array([a[4], a[5]]))
                    return N2.inverse(L(N1))
                except np.linalg.LinAlgError:
                    print('Could not solve least squares fitting system.')
                    return None
            else:
                print('At least three points are needed to find the best fit linear transformation.')
                return None
        else:
            print('The two points lists must contain the same number of points.')
            return None
    def rotation(self, angle, centre = None):
        """Returns a linear transformation representing a rotation by the
        specified angle (degrees clockwise), about an optionally
        specified centre.
        """
        if centre is None: centre = np.zeros(2)
        T = linear_trans2(b = -centre)
        from math import radians, sin, cos
        angleradians = radians(angle)
        cosangle, sinangle = cos(angleradians), sin(angleradians)
        R = linear_trans2(np.array([[cosangle, sinangle], [-sinangle, cosangle]]))
        result = R(T)
        result.b += centre  # i.e. return T.inverse(R(T))
        return result
