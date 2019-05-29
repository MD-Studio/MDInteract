# -*- coding: utf-8 -*-

"""
file: module_load_topology_test.py

Unit tests for the interact topology loader methods
"""

import os

from pandas import DataFrame

from interact.md_system import System
from interact.core.topology_dataframe import TopologyDataFrame

from tests.module.unittest_baseclass import UnittestPythonCompatibility


class LoadTopologyTests(UnittestPythonCompatibility):
    currpath = os.path.dirname(__file__)
    pdb_file = os.path.abspath(os.path.join(currpath, '../files/1acj.pdb'))
    mol_file = os.path.abspath(os.path.join(currpath, '../files/1acj.mol2'))

    def test_load_topology_exceptions(self):
        """
        No file, unknown file or wrong path raises TypeError
        """

        self.assertRaises(TypeError, System, 'file_not_exist')
        self.assertRaises(TypeError, System, os.path.abspath(os.path.join(self.currpath, '../files/tt.csv')))

    def test_load_topology_pdb(self):
        """
        Load PDB file
        """

        top = System(self.pdb_file).topology

        self.assertEqual(len(top), 5202)
        self.assertEqual(list(top.columns), ['serial', 'name', 'element', 'resSeq', 'resName',
                                             'chainID', 'segmentID'])

    def test_load_topology_sybyl(self):
        """
        Load topology including Tripos SYBYL atom types from mol2
        """

        top = System(self.pdb_file, mol2file=self.mol_file).topology
        self.assertTrue('attype' in top.columns)

    def test_load_topology_instance(self):
        """
        Check if returned object is instance of DataFrame, TopologyDataFrame
        """

        top = System(self.pdb_file).topology

        self.assertIsInstance(top, DataFrame)
        self.assertIsInstance(top, TopologyDataFrame)
