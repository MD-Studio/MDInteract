# -*- coding: utf-8 -*-

"""
file: module_topology_series_test.py

Unit tests for the TopologySeries class
"""

import os

from numpy import ndarray

from interact.core.topology_series import TopologySeries
from interact.core.topology_dataframe import TopologyDataFrame
from interact.md_system import System

from tests.module.unittest_baseclass import UnittestPythonCompatibility, PY_PRIMITIVES


class TopologySeriesTests(UnittestPythonCompatibility):
    currpath = os.path.dirname(__file__)
    pdb_file = os.path.abspath(os.path.join(currpath, '../files/1acj.pdb'))
    mol_file = os.path.abspath(os.path.join(currpath, '../files/1acj.mol2'))

    @classmethod
    def setUpClass(cls):
        """
        Prepare TopologyDataFrame once for every test
        """

        cls.top = System(cls.pdb_file, mol2file=cls.mol_file).topology

    def test_topseries_column_select(self):
        """
        Selecting a single DataFrame column should return a TopologySeries
        object
        """

        col = self.top['resName']
        self.assertIsInstance(col, TopologySeries)

    def test_topseries_row_select(self):
        """
        Selecting a single DataFrame row using the custom TopologyDataFrame
        iterrows method should return a TopologySeries.
        """

        for i, a in self.top.iterrows():
            self.assertIsInstance(a, TopologySeries)
            break

    def test_topseries_coord(self):
        """
        Single coordinates or coordinate slices should be available in a
        TopologySeries object depending if it is a single row or single
        column series.
        """

        # Single column TopologySeries
        col = self.top.loc[self.top['resName'] == 'HIS', ['serial']]
        self.assertEqual(len(col.coord), len(col))

        # Single row TopologySeries means one coordinate
        for i, a in self.top[self.top['resName'] == 'HIS'].iterrows():
            self.assertIsInstance(a.coord, ndarray)
            self.assertTrue(a.coord.shape, (3,))

    def test_topseries_attr_inheritance(self):
        """
        The _coordinates, _parent and _distance_matrix attributes should be
        inherited from the parent TopologyDataFrame
        """

        col = self.top.loc[self.top['resName'] == 'HIS', ['resName']]
        row = self.top[self.top['serial'] == 133].squeeze()

        for obj in (self.top, col, row):
            self.assertItemsEqual(obj._metadata, ['_parent', '_coordinates', '_distance_matrix',
                                                  'unitcell_vectors', 'unitcell_lengths', 'unitcell_angles', 'time'])

            # _parent has pointer to toplevel dataframe
            self.assertTrue(hasattr(obj, '_parent'))
            self.assertIsInstance(obj._parent, TopologyDataFrame)
            self.assertEqual(id(obj._parent), id(self.top))

            # _coordinate frame should equal length of parent and have pointer
            # to it
            self.assertTrue(hasattr(obj, '_coordinates'))
            self.assertEqual(id(obj._coordinates), id(self.top._coordinates))
            self.assertEqual(len(obj._coordinates), len(self.top))

            # _distance_matrix not set by default. Should be None
            self.assertTrue(hasattr(obj, '_distance_matrix'))
            self.assertIsNone(obj._distance_matrix)

    def test_topseries_attr_access(self):
        """
        Every data element in a row based TopologySeries should be accessible
        as an attribute
        """

        # Single row attribute access
        for i, a in self.top.iterrows():
            for label in a.axes[0].tolist():
                self.assertIsInstance(getattr(a, label), PY_PRIMITIVES)
            break

    def test_topseries_labels(self):
        """
        Quick access to the labels of a series (axis labels)
        """

        row = self.top[self.top['serial'] == 133].squeeze()
        self.assertItemsEqual(row.labels(), row.axes[0].tolist())

    def test_topseries_neighbours(self):
        """
        Get neighbour selection from single atom
        """

        # Build distance matrix first
        self.top.distances()

        # Get neighbours from 999-THA-N to rest of system within 0.6 nm
        row = self.top[self.top['serial'] == 5194].squeeze()
        ng = row.neighbours()

        self.assertEqual(list(ng['resSeq'].unique()), [84, 199, 330, 439, 440, 441, 442, 634, 999])

    def test_topseries_extend(self):
        """
        Extend atom selection to residue
        """

        row = self.top[self.top['serial'] == 5194].squeeze()
        res = row.extend()
        self.assertItemsEqual(res['serial'], self.top.loc[self.top['resSeq'] == 999, 'serial'])
