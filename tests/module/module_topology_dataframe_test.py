# -*- coding: utf-8 -*-

"""
file: module_load_topology_test.py

Unit tests for the interact topology loader methods
"""

import os

from numpy import ndarray

from interact.md_system import System

from tests.module.unittest_baseclass import UnittestPythonCompatibility
from tests.module.test_data import *

currpath = os.path.dirname(__file__)


class TopologyDataframeTests(UnittestPythonCompatibility):
    pdb_file = os.path.abspath(os.path.join(currpath, '../files/1acj.pdb'))
    mol_file = os.path.abspath(os.path.join(currpath, '../files/1acj.mol2'))

    @classmethod
    def setUpClass(cls):
        """
        Prepare TopologyDataFrame once for every test
        """

        cls.top = System(cls.pdb_file, mol2file=cls.mol_file).topology

        # Build distance matrix
        cls.top.distances()

    def test_topdataframe_parent_pointer(self):
        """
        Pointer to original (full) topology DataFrame should persist
        over selections.
        """

        # Parent of initial frame equals self
        self.assertEqual(self.top._parent.shape, self.top.shape)

        # Parent in selection should equal initial frame
        sel = self.top[self.top['resName'] == 'GLU']
        self.assertEqual(sel._parent.shape, self.top.shape)

        # Selection from a selection
        sel2 = sel[sel['resSeq'].isin([462, 463, 484])]
        self.assertEqual(sel2._parent.shape, self.top.shape)

    def test_topdataframe_contains(self):
        """
        TopologyDataFrame contains method which acts as the __contains__
        magic method. We cannot overload __contains__ as it is used by
        Pandas internally.
        """

        sel1 = self.top[self.top['resName'] == 'ALA']
        sel2 = self.top[self.top['resSeq'] == 36] # ALA
        sel3 = self.top[self.top['resSeq'] == 37] # GLU

        self.assertTrue(sel1.contains(sel2))
        self.assertFalse(sel1.contains(sel3))

    def test_topdataframe_set_coords_exceptions(self):
        """
        Set coordinate matrix should match parent topology in number of atoms.
        Raise TypeError if not.
        Set coordinates should have three columns. Raise AssertionError if not
        """

        sel = self.top[self.top.is_amino_acid()]
        self.assertRaises(TypeError, self.top.set_coord, sel.coord)
        self.assertRaises(AssertionError, self.top.set_coord, [[1, 2], [3, 4]])

    def test_topdataframe_coord(self):
        """
        Get coordinates should return a Pandas DataFrame with three columns
        matching the selection
        """

        orig_coord = self.top.coord
        self.assertIsInstance(orig_coord, ndarray)
        self.assertEqual(len(orig_coord), len(self.top))
        self.assertEqual(orig_coord.shape[1], 3)

        sel = self.top[self.top['resName'] == 'THA']
        sel_coord = sel.coord
        self.assertIsInstance(sel_coord, ndarray)
        self.assertEqual(len(sel_coord), len(sel))
        self.assertEqual(sel_coord.shape[1], 3)

    def test_topdataframe_neighbours_exceptions(self):
        """
        Neighbours method exceptions
        """

        # Method requires pairwise distance matrix
        top = System(self.pdb_file, mol2file=self.mol_file).topology
        sel = top[top['resName'] == 'THA']
        self.assertRaises(AttributeError, sel.neighbours)

        # Target selection should be contained in parent
        sel.distances()
        target = top[top['resSeq'].isin([72, 80, 81, 84, 85, 117, 118, 119, 121, 122])]
        parent = top[top['resSeq'].isin([84, 85, 117, 118, 119, 121, 122])]
        sel._parent = parent
        self.assertRaises(TypeError, sel.neighbours, target)

    def test_topdataframe_neighbours(self):
        """
        Get neighbour selection
        """

        # Get neighbours from THA selection to rest of system within 0.6 nm
        sel = self.top[self.top['resName'] == 'THA']
        ng = sel.neighbours()
        self.assertEqual(list(ng['resSeq'].unique()), [72, 80, 81, 84, 85, 117, 118, 119, 121, 122, 130, 199, 200,
                                                       201, 330, 334, 432, 436, 439, 440, 441, 442, 444, 604, 607,
                                                       609, 616, 624, 634, 643])

    def test_topdataframe_neighbours_toself(self):
        """
        Get neighbours using self as target should yield nothing
        """

        sel = self.top[self.top['resName'] == 'THA']
        self.assertTrue(sel.neighbours(target=sel).empty)

    def test_topdataframe_neighbours_covalent(self):
        """
        Neighbours over atom based selections within covalent distance
        should yield bonds
        """

        sel = self.top[self.top['serial'] == 4106]
        self.assertEqual(list(sel.neighbours(cutoff=0.15)['serial']), [4105, 4107, 4108])

    def test_topdataframe_extend(self):
        """
        Test extending selections selections
        """

        # Default resSeq extend or resName should yield the full ligand
        sel1 = self.top[self.top['serial'].isin([4105, 4104])]
        self.assertTrue(len(sel1.extend()), 15)
        self.assertTrue(len(sel1.extend(mode='resName')), 15)

        # Select TRP 84, NE1
        sel2 = self.top[self.top['serial'] == 617]
        sel3 = sel2.extend()
        self.assertEqual(len(sel3), 14) # Extend to full residue
        self.assertEqual(len(sel2.extend(mode='resName')['resSeq'].unique()), 17) # Extend to all TRP residues
        self.assertEqual(len(sel3.extend(mode='chainID')), 4489)  # Full chain


    def test_topdataframe_contacts(self):
        """
        Build a new contact DataFrame for the source (THA residue) to the
        target which is the full system here.
        """

        source = self.top[self.top['resName'] == 'THA']
        cdf = source.contacts(target=self.top)

        self.assertItemsEqual(list(cdf.columns.levels[0]), [u'source', u'target', u'contact'])
        self.assertItemsEqual(list(cdf.columns.levels[1]), [u'index', u'serial', u'name', u'element', u'resSeq',
                                                            u'resName', u'chainID', u'segmentID', u'attype', u'charge',
                                                            u'distance', u'angle', u''])

        self.assertEqual(len(cdf), 88145)
        self.assertEqual(list(cdf['source', 'resSeq'].unique()), [999])
        self.assertEqual(list(cdf.contact.unique()), ['nd'])

        # Residue numbers in contacts 'target' equal full system without source (THA)
        self.assertEqual(len(cdf['target', 'resSeq'].unique()), len(self.top['resSeq'].unique()) - 1)

    def test_topdataframe_contacts_target(self):
        """
        Build a new contact DataFrame from the source (THA residue) to the
        target which is a selection of residues
        """

        source = self.top[self.top['resName'] == 'THA']
        target = self.top[self.top['resSeq'].isin([84, 85, 117, 118, 119, 121, 122, 130])]
        cdf = source.contacts(target=target)

        self.assertEqual(len(cdf), 1326)
        self.assertEqual(list(cdf['source', 'resSeq'].unique()), [999])
        self.assertItemsEqual(list(cdf['target', 'resSeq'].unique()), [84, 85, 117, 118, 119, 121, 122, 130])

    def test_topdataframe_labels(self):
        """
        Quick access to the labels of a dataframe (columns labels)
        """

        self.assertItemsEqual(self.top.labels(), self.top.columns.tolist())

    def test_topdataframe_residues(self):
        """
        Test residue iterator
        """

        res_sel = self.top[self.top['resSeq'].isin(range(500, 505))]
        for residue in res_sel.residues():
            self.assertEqual(residue.shape, self.top[self.top['resSeq'].isin(residue['resSeq'])].shape)


class RingFinderTests(UnittestPythonCompatibility):

    def load_topology(self, pdb_id):

        pdb_file = os.path.abspath(os.path.join(currpath, '../files/{0}.pdb'.format(pdb_id)))
        mol_file = os.path.abspath(os.path.join(currpath, '../files/{0}.mol2'.format(pdb_id)))

        top = System(pdb_file, mol2file=mol_file).topology
        top.distances()

        return top

    def test_rings_1acj(self):

        top = self.load_topology(reference_1acj['pdb_id'])
        lig = top[top['resSeq'] == reference_1acj['ligand']]

        rings_found = [set(ring['serial']) for ring in lig.find_rings(aromatic=False, check_planar=False)]
        ref_rings = [set(ring) for ring in reference_1acj['rings']]

        self.assertEqual(len(rings_found), len(ref_rings))
        for ring in ref_rings:
            self.assertTrue(ring in rings_found)

    def test_rings_1aku(self):

        top = self.load_topology(reference_1aku['pdb_id'])
        lig = top[top['resSeq'] == reference_1aku['ligand']]

        rings_found = [set(ring['serial']) for ring in lig.find_rings(aromatic=False, check_planar=False)]
        ref_rings = [set(ring) for ring in reference_1aku['rings']]

        self.assertEqual(len(rings_found), len(ref_rings))
        for ring in ref_rings:
            self.assertTrue(ring in rings_found)

    def test_rings_1ay8(self):

        top = self.load_topology(reference_1ay8['pdb_id'])
        lig = top[(top['resSeq'] == reference_1ay8['ligand']) & (top['chainID'] == reference_1ay8['chain'])]

        rings_found = [set(ring['serial']) for ring in lig.find_rings(aromatic=False, check_planar=False)]
        ref_rings = [set(ring) for ring in reference_1ay8['rings']]

        self.assertEqual(len(rings_found), len(ref_rings))
        for ring in ref_rings:
            self.assertTrue(ring in rings_found)

    def test_rings_1bju(self):

        top = self.load_topology(reference_1bju['pdb_id'])
        lig = top[top['resSeq'] == reference_1bju['ligand']]

        rings_found = [set(ring['serial']) for ring in lig.find_rings(aromatic=False, check_planar=False)]
        ref_rings = [set(ring) for ring in reference_1bju['rings']]

        self.assertEqual(len(rings_found), len(ref_rings))
        for ring in ref_rings:
            self.assertTrue(ring in rings_found)

    def test_rings_1bma(self):

        top = self.load_topology(reference_1bma['pdb_id'])
        lig = top[top['resSeq'] == reference_1bma['ligand']]

        rings_found = [set(ring['serial']) for ring in lig.find_rings(aromatic=False, check_planar=False)]
        ref_rings = [set(ring) for ring in reference_1bma['rings']]

        self.assertEqual(len(rings_found), len(ref_rings))
        for ring in ref_rings:
            self.assertTrue(ring in rings_found)

    def test_rings_1eve(self):

        top = self.load_topology(reference_1eve['pdb_id'])
        lig = top[top['resSeq'] == reference_1eve['ligand']]

        rings_found = [set(ring['serial']) for ring in lig.find_rings(aromatic=False, check_planar=False)]
        ref_rings = [set(ring) for ring in reference_1eve['rings']]

        self.assertEqual(len(rings_found), len(ref_rings))
        for ring in ref_rings:
            self.assertTrue(ring in rings_found)
