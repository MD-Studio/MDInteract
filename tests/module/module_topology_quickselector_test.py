# -*- coding: utf-8 -*-

"""
file: module_topology_quickselector_test.py

Unit tests for quick atom selection methods
"""

import os

from interact import System
from interact.core.topology_series import TopologySeries
from tests.module.unittest_baseclass import UnittestPythonCompatibility

currpath = os.path.dirname(__file__)


class TopologyDataframeTests(UnittestPythonCompatibility):
    pdb_file = os.path.abspath(os.path.join(currpath, '../files/dnmt.pdb'))
    mol_file = os.path.abspath(os.path.join(currpath, '../files/dnmt.mol2'))

    @classmethod
    def setUpClass(cls):
        """
        Prepare TopologyDataFrame once for every test
        """

        top = System(cls.pdb_file, mol2file=cls.mol_file)
        cls.pdb = top.topology

    def test_quickselect_amino_acid(self):
        """
        Quick select all amino-acids
        """

        self.assertIsInstance(self.pdb.is_amino_acid(), TopologySeries)
        self.assertItemsEqual(set(self.pdb[self.pdb.is_amino_acid()]['resSeq']), range(1, 328))

    def test_quickselect_amino_acid_backbone(self):
        """
        Quick select the backbone atoms of all amino-acids
        """

        self.assertIsInstance(self.pdb.is_amino_acid_backbone(), TopologySeries)

        backbone = self.pdb[self.pdb.is_amino_acid_backbone()]
        self.assertItemsEqual(set(backbone['resSeq']), range(1, 328))
        self.assertItemsEqual(set(backbone['name']), ('C', 'CA', 'CB', 'N', 'O', 'H', 'HA'))

    def test_quickselect_nucleic_acid(self):
        """
        Quick select all nucleic acids
        """

        nucleic_acids = list(range(401, 414)) + list(range(421, 434))

        self.assertIsInstance(self.pdb.is_nucleic_acid(), TopologySeries)
        self.assertItemsEqual(set(self.pdb[self.pdb.is_nucleic_acid()]['resSeq']), nucleic_acids)

    def test_quickselect_nucleic_acid_backbone(self):
        """
        Quick select the backbone atoms of all nucleic-acids
        """

        self.assertIsInstance(self.pdb.is_nucleic_acid_backbone(), TopologySeries)

        nucleic_acids = list(range(401, 414)) + list(range(421, 434))
        backbone = self.pdb[self.pdb.is_nucleic_acid_backbone()]
        self.assertItemsEqual(set(backbone['resSeq']), nucleic_acids)
        self.assertItemsEqual(set(backbone['name']), ('P', 'O1P', 'O2P', "O5'", "C5'", "1H5'", "2H5'", "C4'", "H4'",
                                                      "C3'", "O3'", "H3'", "C2'", "1H2'", "2H2'", "C1'", "H1'", "O4'"))

    def test_quickselect_ligand(self):
        """
        Quick select ligands being all but amino-acids, nucleic-acid,
        solvent or ions.
        That will be only the SAM residue here.
        """

        self.assertIsInstance(self.pdb.is_ligand(), TopologySeries)
        self.assertTrue(self.pdb.loc[self.pdb.is_ligand(), 'resSeq'].unique() == [434])
        self.assertNotIn('SAM', self.pdb[~self.pdb.is_ligand()]['resName'].unique())

    def test_quickselect_is_ring(self):
        """
        Quick select for rings in the system. Try on system without waters
        to reduce size of required distance matrix
        """

        nonhoh = self.pdb[self.pdb['resName'] != 'HOH']
        nonhoh.distances()

        self.assertItemsEqual(nonhoh.loc[nonhoh.is_ring(), 'resName'].unique(),
                              ['PHE', 'TRP', 'HIS', 'TYR', 'DA', 'DG', 'DC', 'DT', 'SAM'])
