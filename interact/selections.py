# -*- coding: utf-8 -*-

"""
file: selections.py

Functions for performing advanced selections on molecular systems
"""


def select_interface(topology, cutoff=0.6, include_solvent=False):
    """
    Compute atoms part of an interface

    Select all atoms that are part of an interface between multiple chains
    based on a `cutoff` distance. Chains by default include all molecules
    and ions except the bulk solvent for computational efficiency reasons.
    Using `include_solvent` the solvent molecules present at the identified
    interface will be included in the selection.

    :param topology:        atom selection to compute interface for
    :type topology:         :interact:TopologyDataFrame
    :param cutoff:          cutoff distance in nm.
    :type cutoff:           :py:float
    :param include_solvent: Include interface solvent molecules
    :type include_solvent:  :py:bool

    :return:                new TopologyDataFrame with interface selection
    :rtype:                 :interact:TopologyDataFrame
    :raises:                Exception
    """

    # Exclude waters for now
    selection = topology[topology['resName'] != 'HOH']

    # Should be at least 2 chains
    chains = selection['chainID'].unique()
    if len(chains) < 2:
        raise Exception('Topology has fewer than 2 chains. No interface')

    atom_selection = []
    for chain in chains:

        # Select source chain serials and target as everything but the source
        source_index = set(selection[selection['chainID'] == chain]['serial'])
        target_index = set(selection['serial']).difference(source_index)

        # Get slice of contact matrix representing the selection, reformat to row based Dataframe
        contacts = selection._distance_matrix.loc[source_index, target_index]
        contacts = contacts.unstack().reset_index()
        contacts.columns = ['target', 'source', 'distance']

        # Select all atoms with distance to target below cutoff
        atom_selection.extend(set(contacts.loc[contacts['distance'] <= cutoff, 'target']))

    interface = topology[topology['serial'].isin(set(atom_selection))]

    # Include solvents
    if include_solvent:
        solvent = select_solvent(interface, cutoff=cutoff)
        atom_selection = list(solvent['serial']) + list(interface['serial'])
        interface = topology._parent[topology._parent['serial'].isin(set(atom_selection))]

    return interface


def select_solvent(topology, cutoff=0.6):
    """
    Compute solvent molecules with a maximum of `cutoff` distance away from an
    atom selection

    :param topology:        atom selection to compute solvent content of
    :type topology:         :interact:TopologyDataFrame
    :param cutoff:          cutoff distance in nm.
    :type cutoff:           :py:float

    :return:                new TopologyDataFrame with solvent selection
    :rtype:                 :interact:TopologyDataFrame
    """

    solvent = topology._parent[(topology._parent['resName'] == 'HOH') &
                               (topology._parent['name'] == 'O')]

    sol_distances = topology.distances(target=solvent)
    contacts = sol_distances.loc[topology['serial'], solvent['serial']]
    contacts = contacts.unstack().reset_index()
    contacts.columns = ['target', 'source', 'distance']

    interface_solvents = contacts.loc[contacts['distance'] <= cutoff, 'target']
    interface_solvents = topology._parent[topology._parent['serial'].isin(set(interface_solvents))]

    return interface_solvents.extend()
