# -*- coding: utf-8 -*-

"""
file: module_evaluators_test.py

Unit tests for the interact profiler evaluators for a selection of
ligand test cases
"""

import os

from interact.md_system import System
from interact.interactions.hydrofobic import eval_hydrophobic_interactions, eval_pistacking
from interact.interactions.hbonds import eval_hbonds, eval_water_bridges
from interact.interactions.charged import eval_saltbridge, eval_pication
from interact.interactions.halogen import eval_halogen_bonds

from tests.module.unittest_baseclass import UnittestPythonCompatibility
from tests.module.test_data import *

currpath = os.path.dirname(__file__)
filepath = os.path.abspath(os.path.join(currpath, '../files'))


class EvaluatorTestBaseClass(object):

    @classmethod
    def setUpClass(cls):
        """
        Prepare TopologyDataFrame once for every test
        """

        pdb_file = '{0}/{1}.pdb'.format(filepath, cls.ref['pdb_id'])
        mol_file = '{0}/{1}.mol2'.format(filepath, cls.ref['pdb_id'])

        # Load structures (PDB + MOL2) and calculate distance matrix
        trajectory = System(pdb_file, mol2file=mol_file)
        cls.top = trajectory.topology
        cls.top.distances()

        # Select ligand and neighbours, build contact frame
        if 'chain' in cls.ref:
            ligand_selection = cls.top[(cls.top['resSeq'] == cls.ref['ligand']) & (cls.top['chainID'] == cls.ref['chain'])]
        else:
            ligand_selection = cls.top[cls.top['resSeq'] == cls.ref['ligand']]

        cls.neighbour_selection = ligand_selection.neighbours()
        cls.cf = ligand_selection.contacts(cls.neighbour_selection)

    @classmethod
    def tearDownClass(cls):
        """
        Report contact frame
        """

        print(cls.cf[cls.cf['contact'] != 'nd'])

    def test_ligand_neighbours(self):
        """
        Test default distance based contacts (neighbours) between ligand and the
        rest of the system except waters.
        """

        neighbours = self.neighbour_selection[self.neighbour_selection['resName'] != 'HOH']
        self.assertItemsEqual(set(neighbours['resSeq']), self.ref['neighbours'])

    def test_ligand_hydrophobic_interactions(self):
        """
        Test overall hydrophobic contacts based on target residue number.
        Not corrected for pi-stacking but that may happen in the pistacking
        test.
        """

        self.cf = eval_hydrophobic_interactions(self.cf, self.top)
        self.assertItemsEqual(list(self.cf.loc[self.cf['contact'] == 'hf', ('target', 'resSeq')].unique()),
                              self.ref['hydrophobic'])

    def test_ligand_pistacking_interactions(self):
        """
        Test pi- T-stacking interactions based on residue number
        """

        self.cf = eval_pistacking(self.cf, self.top)
        self.assertItemsEqual(list(self.cf.loc[self.cf['contact'] == 'ps', ('target', 'resSeq')].unique()),
            self.ref['ps_stacking'])
        self.assertItemsEqual(list(self.cf.loc[self.cf['contact'] == 'ts', ('target', 'resSeq')].unique()),
            self.ref['ts_stacking'])

    def test_ligand_pication_interactions(self):
        """
        Test pi-cation interations based on residue number
        """

        self.cf = eval_pication(self.cf, self.top)
        self.assertItemsEqual(list(self.cf.loc[self.cf['contact'] == 'pc', ('target', 'resSeq')].unique()),
            self.ref['pi_cation'])

    def test_ligand_hbond_interactions(self):
        """
        Test hydrogen bonded interactions based on residue number
        """

        self.cf = eval_hbonds(self.cf, self.top)

        hb = []
        for i, c in self.cf[self.cf['contact'] != 'nd'].iterrows():
            ct = c.contact.values[0]
            if 'hb-ad' in ct or 'hb-da' in ct:
                hb.append(c['target', 'resSeq'])
        self.assertItemsEqual(set(self.ref['hbond']).intersection(hb), self.ref['hbond'])

    def test_ligand_saltbridge_interactions(self):
        """
        Test salt-bridge interactions based on residue number
        """

        self.cf = eval_saltbridge(self.cf, self.top)

        sb = []
        for i, c in self.cf[self.cf['contact'] != 'nd'].iterrows():
            ct = c.contact.values[0]
            if 'sb-np' in ct or 'sb-pn' in ct:
                sb.append(c['target', 'resSeq'])
        self.assertItemsEqual(set(self.ref['salt_bridge']).intersection(set(sb)), self.ref['salt_bridge'])

    def test_ligand_halogen_interactions(self):
        """
        Test halogen interactions based on residue number
        """

        self.cf = eval_halogen_bonds(self.cf, self.top)


    def test_ligand_water_bridges(self):
        """
        Test hydrogen bonded water bridges between ligand and protein using
        water as intermediate
        """

        self.cf = eval_water_bridges(self.cf, self.top)

        wb = []
        for i, c in self.cf[self.cf['contact'] != 'nd'].iterrows():
            ct = c.contact.values[0]
            if 'wb-ad' in ct or 'wb-da' in ct:
                wb.append(i)

        wbf = self.cf.iloc[wb]
        wbsel = []
        for water in wbf['target', 'resSeq'].unique():
            for res in wbf.loc[(wbf['target', 'resSeq'] == water), ('source', 'resSeq')]:
                if res != self.ref['ligand']:
                    wbsel.append({self.ref['ligand'], water, res})

        wbsel = [s for s in wbsel if s in self.ref['water_bridge']]
        self.assertItemsEqual(wbsel, self.ref['water_bridge'])


class Evaluator1ACJTests(EvaluatorTestBaseClass, UnittestPythonCompatibility):
    ref = reference_1acj


class Evaluator1AKUTests(EvaluatorTestBaseClass, UnittestPythonCompatibility):
    ref = reference_1aku


class Evaluator1AY8Tests(EvaluatorTestBaseClass, UnittestPythonCompatibility):
    ref = reference_1ay8


class Evaluator1BJUTests(EvaluatorTestBaseClass, UnittestPythonCompatibility):
    ref = reference_1bju


class Evaluator1BMATests(EvaluatorTestBaseClass, UnittestPythonCompatibility):
    ref = reference_1bma


class Evaluator1EVETests(EvaluatorTestBaseClass, UnittestPythonCompatibility):
    ref = reference_1eve


class Evaluator2REGTests(EvaluatorTestBaseClass, UnittestPythonCompatibility):
    ref = reference_2reg


class Evaluator2W0STests(EvaluatorTestBaseClass, UnittestPythonCompatibility):
    ref = reference_2w0s
