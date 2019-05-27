# -*- coding: utf-8 -*-

"""
file: unittest_baseclass.py

Python 2/3 unittest compatibility class
"""

import unittest
import sys

version = sys.version_info
MAJOR_PY_VERSION = sys.version_info.major
PY_VERSION = '{0}.{1}'.format(version.major, version.minor)

# Unicode test
UNICODE_TYPE = str
PY_PRIMITIVES = (int, float, bool, str)
if MAJOR_PY_VERSION == 2:
    UNICODE_TYPE = unicode
    PY_PRIMITIVES = (int, float, bool, long, str, unicode)


class UnittestPythonCompatibility(unittest.TestCase):

    def assertItemsEqual(self, expected_seq, actual_seq, msg=None):
        """
        Universal assertItemsEqual method.

        Python 2.x has assertItemsEqual but it is assertCountEqual in 3.x.
        """

        if MAJOR_PY_VERSION == 2:
            return super(UnittestPythonCompatibility, self).assertItemsEqual(expected_seq, actual_seq, msg=msg)
        return super(UnittestPythonCompatibility, self).assertCountEqual(expected_seq, actual_seq, msg=msg)

    @classmethod
    def assertViewEqual(self, expected_seq, actual_seq, msg=None):
        """
        Test equality in items even if they are 'view' based
        """

        return all([t in expected_seq for t in actual_seq]) and all([t in actual_seq for t in expected_seq])

    def assertDictAlmostEqual(self, expected_seq, actual_seq, places=None, msg=None):
        """
        Test if dict values are almost equal for numberic values
        """

        self.assertEqual(expected_seq.keys(), actual_seq.keys())
        return all([self.assertAlmostEqual(expected_seq[key], value, places=places) for key, value in
                    actual_seq.items()])
