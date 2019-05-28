# -*- coding: utf-8 -*-

import logging
import numpy

from pandas import DataFrame, Series, concat
from scipy.spatial.distance import cdist

from interact import constants, reference_data, __module__
from interact.core.topology_base import TopologyBaseClass
from interact.core.topology_series import TopologySeries
from interact.core.sssr import sssr

logger = logging.getLogger(__module__)


class TopologyDataFrame(TopologyBaseClass, DataFrame):
    """
    TopologyDataFrame class

    An extended Pandas DataFrame for working with MDTraj molecular structure
    topologies and associated atom coordinates.

    A TopologyDataFrame object is initiated from a pandas topology DataFrame
    obtained using `mdtraj.System.to_dataframe` method. An internal
    coordinate representation that is persistent over dataframe selections
    is initiated using the `set_coordinates` method.
    """

    def __init__(self, *args, **kwargs):

        super(TopologyDataFrame, self).__init__(*args, **kwargs)

        # Set parent the first time
        if not hasattr(self, '_parent'):
            self._parent = self
        if not hasattr(self, '_coordinates'):
            self._coordinates = None
        if not hasattr(self, '_distance_matrix'):
            self._distance_matrix = None

        self.unitcell_vectors = None
        self.unitcell_lengths = None
        self.unitcell_angles = None
        self.time = None

    def contains(self, other):
        """
        Check if the atom selection in other is contained in self

        Implementing `contains` as magic method `__contains__` is not possible
        as the native __contains__ is used internally by Pandas.

        :param other: Other topology DataFrame
        :type other:  :interact:TopologyDataFrame

        :rtype:       :py:bool
        """

        return set(other.get_index()).issubset(set(self.get_index()))

    def __getitem__(self, item):

        return self._wrapped_pandas_method('__getitem__', item)

    @property
    def _constructor(self):

        return TopologyDataFrame

    @property
    def _constructor_sliced(self):

        return TopologySeries

    def _wrapped_pandas_method(self, mtd, *args, **kwargs):
        """
        Wrap a generic pandas method to ensure it returns a TopologySeries
        """

        result = getattr(super(TopologyDataFrame, self), mtd)(*args, **kwargs)
        if isinstance(result, TopologySeries):
            for name in self._metadata:
                object.__setattr__(result, name, getattr(self, name, None))

        return result

    def squeeze(self, *args, **kwargs):

        return self._wrapped_pandas_method('squeeze', *args, **kwargs)

    def iterrows(self):

        columns = self.columns
        klass = self._constructor_sliced
        for k, v in zip(self.index, self.values):
            s = klass(v, index=columns, name=k)

            for name in self._metadata:
                setattr(s, name, getattr(self, name, None))

            yield k, s

    def labels(self):
        """
        Unified method to return TopologyDataFrame column labels or
        TopologySeries axis labels as list.

        :return:    data labels
        :rtype:     :py:list
        """

        return self.columns.tolist()

    def get_index(self):

        return list(self.index)

    def charge(self, partial=False):
        """
        Get the partial or formal charge of the selection

        Computes charge as the sum of partial charges in the the 'charge'
        column. Rounds the the nearest integer value unless `partial` is
        True.

        :param partial: Return partial charge of selection
        :type partial:  :py:bool

        :return:        charge of selection
        :rtype:         :py:int or :py:float if partial
        """

        if 'charge' not in self.columns:
            raise AttributeError('No charge column')

        partial_charge = self['charge'].sum()
        if not partial:
            return int(round(partial_charge))
        return partial_charge

    def center(self, mass=False):
        """
        Computes 3D coordinates of the geometrical center or the center of mass
        of atoms in the selection.

        Given the atoms of a ring it will calculate the ring center.

        :param mass: calculate center of mass
        :type mass:  :py:bool

        :return:     coordinate representing center
        :rtype:      :numpy:ndarray
        """

        coords = self.coord

        if mass:
            elements = reference_data['element_data']
            atom_mass = Series(elements.atomicMass.values, index=elements.symbol).to_dict()
            scale = numpy.array([atom_mass.get(element, 12.0) for element in self['element']])
        else:
            scale = numpy.ones((len(coords), 1))

        scaled = coords * scale
        return numpy.mean(scaled, axis=0)

    def distances(self, target=None, max_distmat_size=75000000):
        """
        Compute pairwise distance matrix

        Compute a Euclidean distance matrix between all atom pairs of the
        current source selection and the target selection using the
        `scipy.spatial.distance.cdist` method.
        If the `target` is not set it equals the current source selection
        resulting in a square pairwise distance matrix.

        Once build, the distance matrix is referenced in all TopologyDataFrame
        and TopologySeries selections made from the parent frame until
        `distances` is called again or a new set of coordinates is registered.
        Structure your workflow as such that you build the distance matrix
        from the initial parent frame once and reuse it in subsequent
        selections.

        Note:
            Although building a pairwise distance matrix is pretty fast, the
            memory load increases quickly with every new atom added to the
            selection. A full MD system of a biomolecular structure in a box
            with explicit solvent will probably not fit in memory anymore.

            An upper limit on the matrix size (max_distmat_size) is enforced
            to prevent problems.

        :param target:           target atom selection
        :type target:            :interact:TopologyDataFrame
        :param max_distmat_size: Maximum size of pairwise distance matrix to
                                 prevent memory flooding.
        :type max_distmat_size:  :py:int

        :return:                 pairwise distance matrix with rows for source
                                 selection and columns for the target.
        :rtype:                  :pandas:DataFrame
        :raises:                 OverflowError, matrix size > max_distmat_size
        """

        if target is None:
            target = self

        # Get dataframe/series index for source and target set
        source_atoms = self.get_index()
        target_atoms = target.get_index()

        # Restrict size of pairwise matrix
        matrix_size = len(source_atoms) * len(target_atoms)
        if matrix_size > max_distmat_size:
            raise OverflowError('Pairwise distance matrix size to large {0} > {1}'.format(matrix_size,
                                                                                          max_distmat_size))

        # Calculate pairwise distance matrix and make DataFrame out of it
        # using source and target index.
        distances = cdist(self.coord, target.coord, metric='euclidean')
        self._distance_matrix = DataFrame(distances, index=source_atoms, columns=target_atoms)

        # Set in parent if called from child. Should not be needed and
        # creates one more retained object
        if self._parent._distance_matrix is None:
            self._parent._distance_matrix = self._distance_matrix

        return self._distance_matrix

    def extend(self, mode='resSeq'):
        """
        Extend the current atom selection based on similarity in topology
        column types defined by the `mode` argument.
        Selection is always restricted to a match in chainID to avoid
        duplicate selection in for instance a dimer complex

        for example:
            Default `resSeq` mode will extend the current selection with
            all additional atoms that share the same residue number as
            the ones in the selection.

        :param mode:  extend selection criterium
        :type mode:   :py:str

        :return:      new TopologyDataFrame with extended atom selection
        :rtype:       :interact:TopologyDataFrame
        :raises:      AttributeError
        """

        if mode not in self.labels():
            raise AttributeError('TopologyDataFrame has no column named {0}'.format(mode))

        chainid_restriction = set(self['chainID'])

        return self._parent[(self._parent[mode].isin(self[mode].unique())) &
                            (self._parent['chainID'].isin(chainid_restriction))]

    def contacts(self, target=None, intra=False):
        """
        Get the distance between the atoms in the current selection (source)
        and the target.

        If the target is not specified it is set to the full system (_parent)
        without the source.
        The returned contact DataFrame is a representation of the pairwise
        distance matrix between the source and target containing all atom
        pairs. A distance cutoff is not applied but that can be easily done
        by quering on the resulting contact DataFrame.

        :param target:  system DataFrame congaing a target selection
        :type target:   :interact:TopologyDataFrame
        :param intra:   calculate intra-molecular contacts
        :type intra:    :py:bool

        :return:        all pairwise contacts between source and target
        :rtype:         :pandas:DataFrame
        """

        source_index = set(self.get_index())

        # Get intra-selection contacts (target == self)
        if intra:
            target_index = source_index

        # Inter-selection contacts.
        # Get index of source (current selection) and target (without source)
        else:
            if target is not None:
                target_index = set(target.get_index()).difference(source_index)
            else:
                target_index = set(self.get_index()).difference(source_index)

        # Get slice of contact matrix representing the selection, reformat to row based Dataframe
        contacts = self._distance_matrix.loc[source_index, target_index]
        contacts = contacts.unstack().reset_index()
        contacts.columns = ['target', 'source', 'distance']

        # Get selection for source and target from parent, reindex and concatenate into new DataFrame
        source = self._parent.loc[(self._parent.index.isin(contacts['source'])), :].copy()
        source.insert(0, 'index', source.index)
        source = source.loc[contacts['source'], :]
        source.index = range(len(source))
        target = self._parent.loc[(self._parent.index.isin(contacts['target'])), :].copy()
        target.insert(0, 'index', target.index)
        target = target.loc[contacts['target'], :]
        target.index = range(len(target))

        columns = source.columns.tolist()
        contacts_frame = concat([source, target, contacts['distance']], axis=1)
        multi_index = [(['source'] * len(columns) + ['target'] * (len(columns) + 1)), columns * 2 + ['distance']]
        contacts_frame.columns = multi_index
        contacts_frame.columns.names = [0, 1]

        # Add angle column (for contacts with angle constraints)
        contacts_frame['target', 'angle'] = numpy.nan

        # Add a contact column and fill it with 'nd' (type not determined)
        contacts_frame['contact'] = 'nd'

        return contacts_frame.sort_values(by=('source', 'serial'))

    def covalent_bonds(self, cutoff=constants['max_covalent_bond_dist']):
        """
        Return covalently bonded atoms in the selection as ContactFrame

        :param cutoff:   covalent bond upper distance
        :type cutoff:    :py:float

        :return:         Covalently bonded atoms
        :rtype:          :pandas:DataFrame
        """

        # Get all intra selection distances
        cf = self.contacts(intra=True)

        # Get atom pairs within covalent bonded distance
        return cf[(cf[('target', 'distance')] >= 0.05) & (cf[('target', 'distance')] <= cutoff)]

    def is_amino_acid(self):
        """
        Quick selector for residues of type amino acid (aa) according to their
        three-letter code described in the `residue_data` reference set.

        The returned pandas Series object can be used in additional data
        queries.
        For custom amino acid selection use the pandas `isin` selector.

        :return:    boolean series of the same length as the DataFrame
                    indicating if a row is an amino acid.
        :rtype:     :interact:TopologySeries
        """

        data = reference_data['residue_data']
        aa = data.loc[data['type'] == 'aa', 'three-letter']

        return self['resName'].isin(set(aa))

    def is_amino_acid_backbone(self, backbone=('C', 'CA', 'CB', 'N', 'O', 'H', 'HA')):
        """
        Quick selector for amino acid backbone atoms by first selecting
        all amino acids using `is_amino_acid` followed by selecting the
        backbone atoms defined in the `backbone` attribute

        The returned pandas Series object can be used in additional data
        queries.

        :param backbone:  backbone atom names
        :type backbone:   :py:list

        :return:          boolean series of the same length as the DataFrame
                          indicating if a row is an amino acid backbone atom.
        :rtype:           :interact:TopologySeries
        """

        return self.is_amino_acid() & self['name'].isin(backbone)

    def is_nucleic_acid(self):
        """
        Quick selector for residues of type nucleic acid (na) according to
        their two-letter code described in the `residue_data` reference set.

        The returned pandas Series object can be used in additional data
        queries.
        For custom nucleic acids selection use the pandas `isin` selector.

        :return:    boolean series of the same length as the DataFrame
                    indicating if a row is a nucleic acid.
        :rtype:     :interact:TopologySeries
        """

        data = reference_data['residue_data']
        na = list(data.loc[data['type'] == 'na', 'two-letter'])
        na.extend(list(data.loc[data['type'] == 'na', 'three-letter']))

        return self['resName'].isin(set(na))

    def is_nucleic_acid_backbone(self, backbone=('P', 'O1P', 'O2P', "O5'", "C5'", "1H5'", "2H5'", "C4'", "H4'",
                                                 "C3'", "O3'", "H3'", "C2'", "1H2'", "2H2'", "C1'", "H1'", "O4'")):
        """
        Quick selector for nucleic acid backbone atoms by first selecting
        all nucleic acids using `is_nucleic_acid` followed by selecting the
        backbone atoms defined in the `backbone` attribute

        The returned pandas Series object can be used in additional data
        queries.

        :param backbone:  backbone atom names
        :type backbone:   :py:list

        :return:          boolean series of the same length as the DataFrame
                          indicating if a row is an nucleic acid backbone atom.
        :rtype:           :interact:TopologySeries
        """

        return self.is_nucleic_acid() & self['name'].isin(backbone)

    def is_ligand(self):
        """
        Quick selector for ligand residues identified as those residues not
        part of the amino-acid, nucleic-acid, solvent and ion residue/element
        groups.

        The returned pandas Series object can be used in additional data
        queries.
        For custom ligand selection use the pandas `isin` selector.

        :return:    boolean series of the same length as the DataFrame
                    indicating if a row is a ligand.
        :rtype:     :interact:TopologySeries
        """

        data = reference_data['residue_data']
        known_types = list(data.loc[data['type'].isin(('sol', 'ion')), 'two-letter'])
        known_types.extend(list(data.loc[data['type'].isin(('sol', 'ion')), 'three-letter']))
        known_types.extend(list(reference_data['element_data']['symbol']))

        return ~self.is_nucleic_acid() & ~self.is_amino_acid() & ~self['resName'].isin(set(known_types))

    def is_ring(self, **kwargs):
        """
        Quick selector for rings in the system.

        This method provides an accumulated view on all rings in the system.
        Use the `find_rings` method to obtain individual isolated rings in a
        selection.

        Ring type (aromatic, planar) can be further specified using the keyword
        arguments accepted by the interact.core.sssr method

        :param kwargs: keyword arguments passed along to the `sssr` method

        :return:       boolean series of the same length as the DataFrame
                       indicating if a row is part of a ring.
        :rtype:        :interact:TopologySeries
        """

        serials = []
        for residue in self.residues():
            detected = residue.find_rings(**kwargs)
            for ring in detected:
                serials.extend(list(ring.index))

        return self.index.isin(serials)

    def residues(self, extend=False):
        """
        Residue iterator

        Iterate over all residues in the current TopologyDataFrame by residue
        number and yield new TopologyDataFrame with the residue selection.
        If 'extend', return all atoms of the given residue.

        :param extend: extend to all atoms of the residue in case of subset
        :type extend:  :py:bool

        :return:       TopologyDataFrame
        :rtype:        :interact:TopologyDataFrame
        :raises:       AttributeError
        """

        for residue in self['resSeq'].unique():
            residue_frame = self[self['resSeq'] == residue]
            if extend:
                yield residue_frame.extend(mode='resSeq')
            else:
                yield residue_frame

    def find_rings(self, **kwargs):
        """
        Find rings in the structure

        Uses SSSR implemented in the `interact.core.sssr` method for finding
        the Smallest Subset of Smallest Rings.

        :param kwargs:  keyword arguments passed along to the sssr method

        :return:        list of rings as TopologyDataFrames
        :rtype:         :py:list
        """

        return sssr(self, **kwargs)

    def find_charged_centers(self, negative=True, positive=True,
                             neg_atoms=('O.co2', 'S.O2', 'S.3'), pos_atoms=('N.4', 'N.pl3', 'N.ar')):
        """
        Find charged centers in the current residue selection

        Extend the neighbourhood of the input atoms to include covalently bonded
        atoms that are not of element type 'C' and 'H'.

        :param negative:    include negative charged centers
        :type negative:     :py:bool
        :param positive:    include positive charged centers
        :type positive:     :py:bool
        :param neg_atoms:   sybyl atom types of negative charged atoms to include
        :type neg_atoms:    :py:tuple
        :param pos_atoms:   sysbyl atom types of positive charged atoms to include
        :type pos_atoms:    :py:tuple
        """

        centers = []
        for residue in self.residues(extend=True):
            charge = self.charge()
            if charge <= -1 and negative:
                charge_selection = residue[residue['attype'].isin(neg_atoms)]
            elif charge >= 1 and positive:
                charge_selection = residue[residue['attype'].isin(pos_atoms)]
            else:
                continue

            if not charge_selection.empty:
                # Get direct covalent neighbours not of element type C, H
                resseq = charge_selection['resSeq'].unique()[0]
                parent = charge_selection._parent

                n = list(charge_selection.index)
                sel = charge_selection
                while not sel.empty:
                    neighbours = sel.neighbours(cutoff=constants['max_covalent_bond_dist'])
                    neighbours = neighbours[~neighbours['element'].isin(('C', 'H')) & (neighbours['resSeq'] == resseq)]

                    serials = [i for i in list(neighbours.index) if i not in n]
                    n.extend(serials)
                    sel = parent[parent.index.isin(serials)]

                centers.append((parent[parent.index.isin(n)], charge))

        return centers
