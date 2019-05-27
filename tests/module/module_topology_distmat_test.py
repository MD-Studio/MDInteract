# -*- coding: utf-8 -*-

"""
file: module_topology_distmat_test.py

Unit tests for distance matrix computations
"""

import os

from pandas import DataFrame

from interact.md_system import System
from tests.module.unittest_baseclass import UnittestPythonCompatibility


class DistanceMatrixTests(UnittestPythonCompatibility):
    currpath = os.path.dirname(__file__)
    pdb_file = os.path.abspath(os.path.join(currpath, '../files/1acj.pdb'))
    mol_file = os.path.abspath(os.path.join(currpath, '../files/1acj.mol2'))

    def setUp(self):
        """
        Prepare TopologyDataFrame once for every test
        """

        self.top = System(self.pdb_file, mol2file=self.mol_file).topology

    def test_distmat_overflow_exception(self):
        """
        Test OverflowError exception for (too) large distance matrix
        """

        # Unable to compute distance matrix > max_distmat_size
        self.assertRaises(OverflowError, self.top.distances, max_distmat_size=10000)

    def test_distmat_attribute_exception(self):
        """
        Test AttributeError on missing or incomplete coordinates
        """

        # No coordinates
        self.top._coordinates = None
        self.assertRaises(AttributeError, self.top.distances)

    def test_distmat_square(self):
        """
        Test computation of default square matrix
        """

        distmat = self.top.distances()

        self.assertIsInstance(distmat, DataFrame)
        self.assertEqual(distmat.shape[0], distmat.shape[1])
        self.assertEqual(list(distmat.columns), list(distmat.index))

    def test_distmat_target(self):
        """
        Test computation of matrix with custom source and target selection
        """

        source = self.top[self.top['resSeq'] == 999]
        target = self.top[self.top['resName'] == 'HIS']
        distmat = source.distances(target=target)

        self.assertIsInstance(distmat, DataFrame)
        self.assertEqual(distmat.shape, (17, 138))
        self.assertEqual(len(distmat.columns), len(target))
        self.assertEqual(len(distmat.index), len(source))

    def test_distmat_empty_selection(self):
        """
        Test compuation of matrix when (one of) the input selections is empty
        """

        source = self.top[self.top['resSeq'] == 9999]
        target = self.top[self.top['resName'] == 'HIS']

        self.assertTrue(source.distances().empty)
        self.assertTrue(source.distances(target=target).empty)