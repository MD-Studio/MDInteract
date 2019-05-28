# -*- coding: utf-8 -*-

"""
file: topology_series.py

MDInteract customized version of the Pandas Series class.
"""

from pandas import Series

from interact.core.topology_base import TopologyBaseClass


class TopologySeries(TopologyBaseClass, Series):

    @property
    def _constructor(self):

        return TopologySeries

    def labels(self):
        """
        Unified method to return TopologyDataFrame column labels or
        TopologySeries axis labels as list.

        :return:    data labels
        :rtype:     :py:list
        """

        return self.axes[0].tolist()

    def get_index(self):

        return [self.name]

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

        return self._parent[(self._parent[mode] == self[mode]) &
                            (self._parent['chainID'] == self.chainID)]
