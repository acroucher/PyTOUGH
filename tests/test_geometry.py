import unittest
from geometry import *

class geometryTestCase(unittest.TestCase):

    def test_in_polygon(self):
        """in_polygon()"""
        poly = [
            np.array([0., 0.]),
            np.array([1., 0.]),
            np.array([0., 1.]),
        ]
        p = np.array([0.8, 0.9])
        self.assertFalse(in_polygon(p, poly))
        p = np.array([0.4, 0.3])
        self.assertTrue(in_polygon(p, poly))

    def test_in_rectangle(self):
        """in_rectangle()"""
        rect = [
            np.array([5., 4.]),
            np.array([8., 11.]),
        ]
        p = np.array([0.8, 0.9])
        self.assertFalse(in_rectangle(p, rect))
        p = np.array([7., 10.5])
        self.assertTrue(in_rectangle(p, rect))

    def test_rect_intersect(self):
        """rectangles_intersect()"""
        r1 = [
            np.array([5., 4.]),
            np.array([8., 11.]),
        ]
        r2 = [
            np.array([7., 5.]),
            np.array([13., 10.]),
        ]
        self.assertTrue(rectangles_intersect(r1, r2))
        self.assertTrue(rectangles_intersect(r2, r1))
        r2 = [
            np.array([1., 2.]),
            np.array([4.6, 10.]),
        ]
        self.assertFalse(rectangles_intersect(r1, r2))
        self.assertFalse(rectangles_intersect(r2, r1))

    def test_bounds_of_points(self):
        """bounds_of_points()"""
        pts = [
            np.array([-1.5, 2.1]),
            np.array([1., -1.1]),
            np.array([4., 3.])
        ]
        bds = bounds_of_points(pts)
        self.assertTrue(np.allclose(np.array([-1.5, -1.1]), bds[0]))
        self.assertTrue(np.allclose(np.array([4., 3.]), bds[1]))

    def test_polygon_area(self):
        """poygon_area()"""
        poly = [
            np.array([0., 0.]),
            np.array([1., 0.]),
            np.array([0., 1.]),
        ]
        self.assertAlmostEqual(0.5, polygon_area(poly))
        poly = [
            np.array([1., 2.]),
            np.array([3., 2.]),
            np.array([3., 4.]),
            np.array([1., 4.]),
        ]
        self.assertAlmostEqual(4., polygon_area(poly))
        t = linear_trans2().rotation(30.)
        poly = [t(p) for p in poly]
        self.assertAlmostEqual(4., polygon_area(poly))

    def test_polygon_centroid(self):
        """polygon_centroid()"""
        poly = [
            np.array([1., 2.]),
            np.array([3., 2.]),
            np.array([3., 4.]),
            np.array([1., 4.]),
        ]
        c = polygon_centroid(poly)
        self.assertTrue(np.allclose(np.array([2., 3.]), c))
        poly = [
            np.array([0., 0.]),
            np.array([15., 0.]),
            np.array([15., 5.]),
            np.array([10., 5.]),
            np.array([10., 10.]),
            np.array([0., 10.])
        ]
        c = polygon_centroid(poly)
        self.assertTrue(np.allclose(np.array([6.5, 4.5]), c))

    def test_line_polygon_intersect(self):
        """line_polygon_intersections()"""
        poly = [
            np.array([0., 0.]),
            np.array([15., 0.]),
            np.array([15., 5.]),
            np.array([10., 5.]),
            np.array([10., 10.]),
            np.array([0., 10.])
        ]
        line = [
            np.array([0., 0.]),
            np.array([15., 4.])
        ]
        pts = line_polygon_intersections(poly, line)
        self.assertEqual(2, len(pts))
        self.assertTrue(np.allclose(np.array([0., 0.]), pts[0]))
        self.assertTrue(np.allclose(np.array([15., 4.]), pts[1]))
        line = [
            np.array([0., 17.5]),
            np.array([17.5, 0.])
        ]
        pts = line_polygon_intersections(poly, line)
        self.assertEqual(4, len(pts))
        self.assertTrue(np.allclose(np.array([7.5, 10.]), pts[0]))
        self.assertTrue(np.allclose(np.array([10., 7.5]), pts[1]))
        self.assertTrue(np.allclose(np.array([12.5, 5.]), pts[2]))
        self.assertTrue(np.allclose(np.array([15., 2.5]), pts[3]))

    def test_simplify_polygon(self):
        """simplify_polygon()"""
        poly = [
            np.array([0., 0.]),
            np.array([1., 0.]),
            np.array([1., 1.]),
            np.array([0.5, 1.]),
            np.array([0., 1.])
        ]
        s = simplify_polygon(poly)
        self.assertEqual(4, len(s))
        poly = [
            np.array([0., 0.]),
            np.array([10., 0.]),
            np.array([15., 0.]),
            np.array([15., 5.]),
            np.array([11., 5.]),
            np.array([10., 5.]),
            np.array([10., 10.]),
            np.array([2., 10.]),
            np.array([0., 10.])
        ]
        s = simplify_polygon(poly)
        self.assertEqual(6, len(s))
        
    def test_polygon_boundary(self):
        """polygon_boundary()"""
        poly = [
            np.array([0., 0.]),
            np.array([15., 0.]),
            np.array([15., 5.]),
            np.array([10., 5.]),
            np.array([10., 10.]),
            np.array([0., 10.])
        ]
        p1 = np.array([5., -1.])
        p2 = np.array([5., 12.])
        b = polygon_boundary(p1, p2, poly)
        self.assertTrue(np.allclose(np.array([5., 0.]), b))
        p1 = np.array([5., 2.5])
        p2 = np.array([20., 2.5])
        b = polygon_boundary(p1, p2, poly)
        self.assertTrue(np.allclose(np.array([15., 2.5]), b))
        b = polygon_boundary(p2, p1, poly)
        self.assertTrue(np.allclose(np.array([15., 2.5]), b))

    def test_line_projection_distance(self):
        """line_projection() and point_line_distance()"""
        from math import sqrt
        line = [
            np.array([0., 0.]),
            np.array([1., 1.])
        ]
        a = np.array([1., 0.])
        p, xi = line_projection(a, line, True)
        self.assertTrue(np.allclose(np.array([0.5, 0.5]), p))
        self.assertAlmostEqual(0.5, xi)
        d = point_line_distance(a, line)
        self.assertAlmostEqual(sqrt(0.5), d)
        a = np.array([2., 1.])
        p, xi = line_projection(a, line, True)
        self.assertTrue(np.allclose(np.array([1.5, 1.5]), p))
        self.assertAlmostEqual(1.5, xi)
        d = point_line_distance(a, line)
        self.assertAlmostEqual(sqrt(0.5), d)

    def test_line_intersects_rectangle(self):
        """line_intersects_rectangle()"""
        r = [np.array([1., 1.]), np.array([5., 3.])]
        line = [np.array([0., 0.]), np.array([1., 1.])]
        self.assertTrue(line_intersects_rectangle(r, line))
        line = r
        self.assertTrue(line_intersects_rectangle(r, line))
        line = r[::-1]
        self.assertTrue(line_intersects_rectangle(r, line))
        line = [np.array([1., 1.]), np.array([5., 1.])]
        self.assertTrue(line_intersects_rectangle(r, line))
        line = [np.array([3., 1.]), np.array([4., 3.])]
        self.assertTrue(line_intersects_rectangle(r, line))
        line = [np.array([3., 0.]), np.array([4., 6.])]
        self.assertTrue(line_intersects_rectangle(r, line))
        line = [np.array([2., 2.]), np.array([4., 2.5])]
        self.assertTrue(line_intersects_rectangle(r, line))
        line = [np.array([6., 6.]), np.array([8., 0.])]
        self.assertFalse(line_intersects_rectangle(r, line))

    def test_vector_heading(self):
        """vector_heading()"""
        from math import atan
        p = np.array([0., 10.])
        h = vector_heading(p)
        self.assertAlmostEqual(0., h)

        p = np.array([1., 1.])
        h = vector_heading(p)
        self.assertAlmostEqual(0.25 * np.pi, h)

        p = np.array([-2., -2.])
        h = vector_heading(p)
        self.assertAlmostEqual(5. * np.pi / 4., h)

        p = np.array([-3., 4.])
        h = vector_heading(p)
        self.assertAlmostEqual(1.5 * np.pi + atan(4./ 3.), h)

    def test_linear_trans(self):
        """Linear transformations"""
        from math import sqrt
        r1 = [
            np.array([5., 4.]),
            np.array([8., 11.]),
        ]
        r2 = [
            np.array([7., 5.]),
            np.array([13., 10.]),
        ]
        r1c = sum(r1) / len(r1)
        r2c = sum(r2) / len(r2)

        t = linear_trans2().between_rects(r1, r1)
        self.assertTrue(np.allclose(r1c, t(r1c)))

        t = linear_trans2().between_rects(r1, r2)
        self.assertTrue(np.allclose(r2c, t(r1c)))

        ti = t.inverse
        self.assertTrue(np.allclose(r1c, ti(r2c)))

        pts1 = [
            np.array([5., 4.]),
            np.array([8., 4.]),
            np.array([5., 11.])
        ]
        pts2 = [
            np.array([7., 5.]),
            np.array([13., 5.]),
            np.array([7., 10.])
        ]
        t = linear_trans2().between_points(pts1, pts2)
        self.assertTrue(np.allclose(r2c, t(r1c)))

        t = linear_trans2().rotation(45., np.ones(2))
        p = np.array([1., 0.])
        sq2 = sqrt(2.)
        a = (sq2 - 1.) / sq2
        p1 = np.array([a, a])
        self.assertTrue(np.allclose(p1, t(p)))
        
        
if __name__ == '__main__':

    suite = unittest.TestLoader().loadTestsFromTestCase(geometryTestCase)
    unittest.TextTestRunner(verbosity = 1).run(suite)
