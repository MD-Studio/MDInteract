# -*- coding: utf-8 -*-

"""
file: module_geometry_test.py

Unit tests for the geometric calculation functions
"""

from interact.core.geometry import *
from tests.module.unittest_baseclass import UnittestPythonCompatibility


class GeometryTests(UnittestPythonCompatibility):

    def test_geometry_distance(self):
        """
        Test calculation of Euclidean distance between two coordinates
        """

        coor1 = numpy.array([-11.719, 89.719, 35.174])
        coor2 = numpy.array([-11.260, 89.101, 37.633])

        self.assertAlmostEqual(distance(coor1, coor2), 2.577, 3)
        self.assertEqual(distance(coor1, coor1), 0.0)

    def test_geometry_angle(self):
        """
        Test calculation of regular angle between three coordinates
        """

        coor1 = numpy.array([-11.066, 88.196, 36.842])
        coor2 = numpy.array([-12.189, 87.877, 35.866])
        coor3 = numpy.array([-12.025, 86.720, 34.856])

        self.assertAlmostEqual(angle(coor1, coor2, coor3), 119.87, 2)
        self.assertAlmostEqual(angle(coor1, coor2, coor3, deg=False), 2.09, 2)

    def test_geometry_angle_exception(self):
        """
        Test exception for angle calculation between two identical vectors
        """

        coor1 = numpy.array([-11.066, 88.196, 36.842])
        coor2 = numpy.array([-12.189, 87.877, 35.866])

        self.assertRaises(ArithmeticError, angle, coor1, coor2, coor1)

    def test_geometry_dihedral(self):
        """
        Test calculation of regular dihedral between four coordinates
        """

        coor1 = numpy.array([-11.260, 89.101, 37.633])
        coor2 = numpy.array([-11.066, 88.196, 36.842])
        coor3 = numpy.array([-9.920, 87.515, 36.768])
        coor4 = numpy.array([-9.840, 86.829, 36.031])

        self.assertAlmostEqual(dihedral(coor1, coor2, coor3, coor4), -175.64, 2)
        self.assertAlmostEqual(dihedral(coor4, coor3, coor2, coor1), -175.64, 2)

        coor1 = numpy.array([-12.025, 86.720, 34.856])
        coor2 = numpy.array([-12.189, 87.877, 35.866])
        coor3 = numpy.array([-11.066, 88.196, 36.842])
        coor4 = numpy.array([-9.920, 87.515, 36.768])

        self.assertAlmostEqual(dihedral(coor1, coor2, coor3, coor4), 3.47, 2)

    def test_geometry_is_planar(self):
        """
        Test evaluation of coordinate (ring/plane) planarity
        """

        planar = numpy.array([[4.207, 88.596, 46.587],
                              [5.469, 89.006, 47.056],
                              [3.390, 89.502, 45.897],
                              [5.934, 90.303, 46.813],
                              [3.861, 90.803, 45.650],
                              [5.131, 91.205, 46.099]])

        self.assertTrue(is_planar(planar))

        notplanar = numpy.array([[5.207, 69.821, 64.683],
                                 [3.833, 69.634, 64.995],
                                 [5.776, 69.459, 63.289],
                                 [4.736, 69.363, 62.145],
                                 [3.515, 68.540, 62.597],
                                 [2.860, 69.177, 63.848]])

        self.assertFalse(is_planar(notplanar))
