# -*- coding: utf-8 -*-

"""
file: topology_base.py

Baseclass shared by TopologyDataFrame and TopologySeries classes
"""

import collections

from pandas import DataFrame, Series

from interact import constants


class TopologyBaseClass(object):

    _metadata = ['_parent', '_coordinates', '_distance_matrix', 'unitcell_vectors', 'unitcell_lengths',
                 'unitcell_angles', 'time']

    def __finalize__(self, other, method=None):
        """
        Called at the end of DataFrame.__init__ to initiate and propagate
        metadata between the parent DataFrame and DataFrame selections.

        Set a pointer to the _parent DataFrame at first initiation of the
        TopologyDataFrame when the _metadata dictionary is still empty.

        :param other:   parent data frame
        :type other:    :pandas:DataFrame
        """

        for name in self._metadata:
            object.__setattr__(self, name, getattr(other, name, None))

        return self

    @property
    def coord(self):
        """
        Attribute based access to atom 3D coordinate

        Coordinates are selected based on atom index.

        :return:    3D coordinate
        :rtype:     :numpy:ndarray
        """

        if self._coordinates is None:
            raise AttributeError('Parent TopologyDataFrame has no coordinates, use set_coord first')

        coord = self._coordinates.iloc[self.get_index()]
        if len(coord) == 1:
            return coord.values[0]

        if coord._is_view:
            return coord.copy().values
        return coord.values

    def set_coord(self, xyz):
        """
        Set atom coordinates from XYZ numpy array or DataFrame

        Atom coordinates are maintained as seperate DataFrame in the
        TopologyDataFrame object to enable fast switching of coordinate sets
        in MD trajectories for instance.

        The only requirement when loading a new coordinate set is a match in
        the number of atoms in parent topology and the coordinates.
        Existing distance matrix objects are removed.

        :param xyz: atom coordinates
        :type xyz:  :numpy:ndarray
        :raises:    TypeError if topology and coordinates do not share the same
                    number of atoms.
        """

        # Load coordinates in pandas DataFrame. Enforce three columns
        coord_df = DataFrame(xyz, columns=['x','y','z'])

        # Check if coordinate frame matches topology
        if len(self._parent) != len(coord_df):
            raise TypeError('Coordinate frame ({0}) and topology ({1}) do not share the same number of atoms'.format(
                len(coord_df), len(self._parent)))

        # Set the coord_df index to equals atom 'index' column
        coord_df.set_index(self._parent.index.values, inplace=True)
        self._coordinates = coord_df

        # Clear former distance matrix
        if self._distance_matrix is not None:
            self._distance_matrix = None

    def neighbours(self, target=None, covalent=False, cutoff=0.6):
        """
        Return the neighbours of the atoms in the current TopologyDataFrame
        selection (source) with respect to the full system or the `target` if
        defined.

        Selecting covalently linked neighbour atoms can be done by reducing
        the cutoff distance to covalent bond length which is up to 0.17 nm.
        The 'covalent' argument functions as a shortcut method here, changing
        the cutoff the the global max_covalent_bond_dist distance and filtering
        the results on chainID ensuring that neighbouring atoms are part of the
        same chain as the source selection.

        :param target:   Atom selection to use as target for the neighbour
                         search. Needs to be derived from the same parent
                         topology as the source selection.
        :type target:    :interact:TopologyDataFrame
        :param covalent: select covalently linked neighbours only
        :type covalent:  :py:bool
        :param cutoff:   distance cutoff for neighbour selection (nm)
        :type cutoff:    :py:float

        :return:         neighbour selection
        :rtype:          :interact:TopologyDataFrame
        :raises:         AttributeError, when pairwise distance matrix is not
                         set. TypeError, when target selection not contained
                         in self.
        """

        # Neighbours requires a pairwise distance matrix
        if self._distance_matrix is None:
            raise AttributeError('Use set_distance_matrix to build a pairwise distance matrix first')

        # Get index of source (current selection) and target (without source)
        source = self.get_index()
        source = set(source)

        if target is not None:
            if not self.contains(target):
                raise TypeError('Target selection not contained in topology (selection')
            target = set(target.get_index()).difference(source)
        else:
            target = set(self._parent.get_index()).difference(source)

        # Set cutoff in case of covalent but only when default cutoff used
        if covalent and cutoff == 0.6:
            cutoff = constants['max_covalent_bond_dist']

        # Get slice of contact matrix for source to target within cutoff distance
        contacts = self._distance_matrix.loc[target, source]
        contacts = contacts[contacts <= cutoff].dropna(how='all')

        selection = self._parent[self._parent.index.isin(contacts.index)]

        # If covalent, ensure that chainID of selection equals source
        if covalent:
            chainID = self['chainID']
            if not isinstance(chainID, collections.Iterable):
                chainID = [chainID]
            return selection[selection['chainID'].isin(set(chainID))]
        return selection
