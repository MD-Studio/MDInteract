# -*- coding: utf-8 -*-

"""
file: sssr.py

Implements method for resolving the Smallest Set of Smallest Rings (SSSR) using
graph algorithms.
"""

import logging

from interact import constants, __module__
from interact.core.geometry import is_planar
from interact.core.helpers import is_aromatic

logger = logging.getLogger(__module__)

AROMATIC_RINGS = {'PHE': [('CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ')],
                  'HIS': [('ND1', 'NE2', 'CE1', 'CG', 'CD2')],
                  'TRP': [('CD2', 'CE2', 'CE3', 'CH2', 'CZ2', 'CZ3'), ('CD2', 'CE2', 'CG', 'NE1', 'CD1')],
                  'TYR': [('CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ')],
                  'DC': [('O2', 'N3', 'C4', 'C5', 'C6', 'N1')],
                  'DT': [('N1', 'N3', 'C2', 'C4', 'C5', 'C6')],
                  'DA': [('N1', 'C2', 'C6', 'N3', 'C4', 'C5'), ('C4', 'C5', 'N7', 'N9', 'C8')],
                  'DG': [('N1', 'C2', 'C6', 'N3', 'C4', 'C5'), ('C4', 'C5', 'N7', 'N9', 'C8')]}


def sssr(system, aromatic=True, check_planar=True, planarity_dihedral_max_div=7.5, maxiter=1000):
    """
    Find rings in the structure

    This function uses graph based cycle detection to find the set of
    smallest unique rings in the structure. The graph algorithm will return
    all rings in the system regardless there nature. Optional filters can
    remove rings that are not planar and/or ar not aromatic based on SYBYL
    atom type.

    PROBLEMS: Fused rings that form a 3D box are not well recovered.
    Fused rings like estrogen's always miss out on one ring, the aromatic
    ring is found.

    Algorithm:
    1) Select all heavy atoms in the system. attype != H
    2) Build the adjacency matrix using a default bond length cutoff of
       `max_covalent_bond_dist`. Perform a check to see if the structure
       contains sulphur. Adjust the cutoff to 0.181 nm.
    3) Construct a boolean matrix from the adjacency matrix to identify
       bonded neighbours (1) from non bonded ones (0).
    4) Iteratively remove all atoms having only one neighbour
    5) Build a graph representation from the boolean matrix with the atoms
       as nodes and covalently bonded neighbours as edge list.
    6) Iterate over the nodes of the graph. for each node evaluate
       connectivity using depth-first search (dfs). Minimize number of
       recursions by keeping track of nodes already visited.
    7) Evaluate the spanning tree while it grows to see if there
       can be a path created back to the ancestor (cycle).
    8) For each ring (cycle) found, recreate the visited nodes dictionary
       only with the nodes having two connections. This prevents evaluating
       atoms that are part of a ring already found but still allow atoms
       that may be part of fused rings.
    9) Repeat step 5 till 8 for a graph in which the neighbour list is
       reversed to evaluate connectivity in the opposite direction and
       minimize the changes of a deadlock.
    10) Filter the list of rings: remove duplicate rings, remove fused
        rings for which the individual members have also been found
        (smallest set of smallest rings). If `aromatic` filter rings based
        on aromaticity (rings containing only C.ar or N.ar atom types).
        If check_planar equals true, filter rings that are not planar based
        on any ring dihedral deviating more then planarity_dihedral_max_div
        from planar 180 deg.

    :param system:                     molecular system containing rings
    :type system:                      :interact:TopologyDataFrame
    :param aromatic:                   ring is aromatic
    :type aromatic:                    :py:bool
    :param check_planar:               Check for ring planarity.
    :type check_planar:                :py:bool
    :param planarity_dihedral_max_div: planarity cutoff for aromatic rings
    :type planarity_dihedral_max_div:  :py:float
    :param maxiter:                    maximum number of iterations for the
                                       removal of single bonded atoms.
    :type maxiter:                     :py:int

    :return:                           list of rings as TopologyDataFrames
    :rtype:                            :py:list
    """

    # Get all heavy atoms of the system
    heavyatoms = system[system['attype'] != 'H']

    # Get heavyatoms distance matrix slice within covalent bond distance
    # range max_covalent_bond_dist. Adjust when sulphur atom in selection
    bond_cutoff = constants['max_covalent_bond_dist']
    if 'S' in heavyatoms['element'].unique():
        bond_cutoff = 0.181
        logger.debug("Detected sulphur atom, adjust covalent bond length cutoff to {0:.3f} nm".format(bond_cutoff))

    adjmatr = system._distance_matrix.loc[heavyatoms.index, heavyatoms.index]
    boolmatr = adjmatr[(adjmatr > 0) & (adjmatr < bond_cutoff)].notnull().astype(int)
    logger.debug("{0} covalently linked heavy atoms in structure".format(len(adjmatr)))

    # Iteratively remove all atoms with one neighbour
    singles = True
    itr = 0
    while singles and itr != maxiter:
        onebond = boolmatr.loc[boolmatr.sum(axis=1) == 1]
        if not onebond.empty:
            indexes = set(boolmatr.index.values).difference(set(onebond.index.values))
            boolmatr = boolmatr.loc[indexes, indexes]
            itr += 1
        else:
            singles = False

    logger.debug(
        "Removed {0} terminal, non-cyclic atoms in {1} iterations".format(adjmatr.shape[0] - boolmatr.shape[0],
                                                                          itr))

    def find_cycle_to_ancestor(node, ancestor):
        """
        Find a cycle containing both node and ancestor.
        """

        path = []
        while node != ancestor:
            if node is None:
                return []
            path.append(node)
            node = spanning_tree[node]
        path.append(node)
        path.reverse()

        return path

    def dfs(node):
        """
        Depth-first search subfunction.
        """

        visited[node] = len(graph[node])
        # Explore recursively the connected component
        for each in graph[node]:
            if ring:
                return
            if each not in visited:
                spanning_tree[each] = node
                dfs(each)
            else:
                if spanning_tree[node] != each:
                    ring.extend(find_cycle_to_ancestor(node, each))

    # Efficient ring detection using Depth-first search (dfs).
    # Maintains a list of visited atoms to prevent revisit of atoms
    # having only two connections and previously found to be part of
    # a ring.
    rings = []
    for pathset in ('forward', 'reversed'):

        # Create a graph representation of the structure (atoms are nodes, bonds the edges)
        if pathset == 'forward':
            graph = dict([(ida, list(a[a == 1].index.values)) for ida, a in boolmatr.iterrows()])
        else:
            graph = dict([(ida, list(a[a == 1].index.values[::-1])) for ida, a in boolmatr.iterrows()])

        visited = {}
        spanning_tree = {}
        subrings = []
        ring = []
        for atom in graph.keys():
            # Select a non-visited node
            if atom not in visited:
                spanning_tree[atom] = None
                dfs(atom)  # Explore atom connections
                if ring:
                    subrings.append(ring)
                    visited = dict([(a, 2) for c in subrings for a in c if len(graph[a]) == 2])
                    ring = []  # Reset cycle list

        rings.extend(subrings)

    if not rings:
        logger.debug("No rings found.")
        return []

    # Filter cycles list to remove duplicates but maintain cycle
    # linkage order.
    filtered_rings = []
    indexes = []
    for i, c in enumerate(rings):
        c = set(c)

        superset = False
        for r2 in rings:
            if c.issuperset(set(r2)) and len(c) > len(r2):
                logger.debug("Detected ring {0} superset of smaller rings. Remove".format(list(c)))
                superset = True
                break

        if c not in filtered_rings and not superset:
            filtered_rings.append(c)
            indexes.append(i)

    # Determine ring planarity and check for aromaticity
    ringlist = []
    for r in indexes:

        # Get ring selection
        ring_selection = system[system.index.isin(rings[r])]

        # Check aromaticity
        if aromatic and not is_aromatic(ring_selection, max_div=planarity_dihedral_max_div):
            logger.debug('Ring {0} not aromatic'.format(rings[r]))
            continue

        # Check planarity
        elif check_planar and not is_planar(ring_selection.coord, max_div=planarity_dihedral_max_div):
            logger.debug('Ring {0} not planar'.format(rings[r]))
            continue

        ringlist.append(ring_selection)

    return ringlist
