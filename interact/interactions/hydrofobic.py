# -*- coding: utf-8 -*-

import logging

from interact import constants, __module__
from interact.core.helpers import remove_contact_type, set_contact_type
from interact.interactions.utils import filter_contacts, is_pistack

logger = logging.getLogger(__module__)


def eval_hydrophobic_interactions(contact_frame, topology, hydroph_dist_max=0.4, min_dist=0.05,
                                  carbon_types=('C.3', 'C.2', 'C.1', 'C.ar')):
    """
    Evaluate hydrophobic interactions between source and target selections

    Contacts are marked hydrophobic ('hf') if both atoms involved are carbons
    in `carbon_types` only having carbon or hydrogen atoms as neighbours and
    their distance is min_dist < distance <  hydroph_dist_max.

    The number of atom pairs involved in hydrophobic contacts can grow quit
    rapidly even surpassing the total number of other identified functional
    contacts. Two filter steps reduce the number of hydrophobic contacts:

    - If one ligand atom contacts multiple target atoms, the one with the
      smallest distance is kept.
    - If atom pairs involved in hydrophobic interactions are also involved in
      pi-stacking interactions they should not be labeled as 'hf' as the
      stacking interaction is also a form of hydrophobic interaction.
      This requires the `eval_pistacking` function to be run following the
      `eval_hydrophobic_interactions` to correct for this.

    :param contact_frame:     contacts DataFrame
    :type contact_frame:      :pandas:DataFrame
    :param topology:          main system topology
    :type topology:           :interact:TopologyDataFrame
    :param hydroph_dist_max:  maximum distance for hydrophobic interactions (nm)
    :type hydroph_dist_max:   :py:float
    :param min_dist:          minimum interaction distance (nm)
    :type min_dist:           :py:float
    :param carbon_types:      SYBYL carbon atom types
    :type carbon_types:       :py:tuple

    :return:                  contacts DataFrame
    :rtype:                   :pandas:DataFrame
    """

    # Select all `carbon_types` atom pairs below hydroph_dist_max contact distance
    hfobdist = contact_frame[(contact_frame['source', 'attype'].isin(carbon_types)) &
                             (contact_frame['target', 'attype'].isin(carbon_types)) &
                             (contact_frame['target', 'distance'] > min_dist) &
                             (contact_frame['target', 'distance'] < hydroph_dist_max)]

    logger.info("Evaluate {0} possible hydrophobic contacts using: hydroph_dist_max={1} nm".format(
        len(hfobdist), hydroph_dist_max))

    # Only select atoms that have hfob_atom_list atoms as neighbours
    hfob_idx = []
    hfob_index = hfobdist.index
    hfob_atom_list = set(list(carbon_types) + ['H'])
    for n, pair in enumerate(hfobdist[[('source', 'serial'), ('target', 'serial')]].values):
        source_neighbours = topology[topology['serial'] == pair[0]].neighbours(cutoff=constants['max_covalent_bond_dist'])
        target_neighbours = topology[topology['serial'] == pair[1]].neighbours(cutoff=constants['max_covalent_bond_dist'])
        if len(set(source_neighbours['attype']).difference(hfob_atom_list)) == 0 and \
           len(set(target_neighbours['attype']).difference(hfob_atom_list)) == 0:
            hfob_idx.append(hfob_index[n])

    # For pairs involving same ligand atom and multiple atoms of same target
    # protein residue and vice versa, retain only the pair with the smallest
    # distance.
    hf = contact_frame.loc[hfob_idx]
    hf_filtered = filter_contacts(hf, source_sel='serial', target_sel='resSeq')
    logger.info('Filtered {0} redundant hydrophobic contacts'.format(len(hf) - len(hf_filtered)))

    for index in hf_filtered.index:
        contact_frame.loc[index, 'contact'] = set_contact_type(contact_frame.loc[index, 'contact'], 'hf')

    return contact_frame


def eval_pistacking(contact_frame, topology, pistack_dist_max=0.55, pistack_ang_dev=30.0, min_dist=0.05,
                    pistack_offset_max=0.20):
    """
    Evaluate pi- and T-stacking between aromatic rings

    A pi- or T-stacking interaction is identified if the geometric center of
    the two aromatic rings are within pistack_dist_max from each other; the
    offset distance between the geometric center of one ring projected onto
    the other is no more than pistack_offset_max and if the angle between the
    normals of both rings does not deviate more than pistack_ang_dev from 180
    deg for pi-stacking or 90+/-offset for T-stacking

    Identified pi- or T-stacking interactions are added to the contact_frame as
    contact between the geometric centers of the two rings represented by two
    dummy atoms Du, (element D) of residue X1. The target distance is the 3D
    euclidean distance between the two dummy atoms and the target angle the
    smallest angle between the two ring normals. Pi-stacking contacts are
    marked 'ps' and T-stacking 'ts'

    :param contact_frame:      contacts DataFrame
    :type contact_frame:       :pandas:DataFrame
    :param topology:           main system topology
    :type topology:            :interact:TopologyDataFrame
    :param pistack_dist_max:   Cutoff distance between aromatic
                               center-of-masses(nm)
    :type pistack_dist_max:    :py:float
    :param pistack_ang_dev:    Maximum angle variation between ring normals
                               (deg). Deviation from 180 for pi-stacking and 90
                               for T-stacking
    :type pistack_ang_dev:     :py:float
    :param pistack_offset_max: Cutoff distance between geometric centers
                               projected on top of each other.
    :type pistack_offset_max:  :py:float
    :param min_dist:           minimum interaction distance (nm)
    :type min_dist:            :py:float

    :return:                   contacts DataFrame
    :rtype:                    :pandas:DataFrame
    """

    # Search for rings in source and target based on all residues
    # containing C.ar and N.ar atom types
    ringsel = contact_frame[(contact_frame['source', 'attype'].isin(('C.ar', 'N.ar', 'C.2'))) &
                            (contact_frame['target', 'attype'].isin(('C.ar', 'N.ar'))) &
                            (contact_frame['target', 'distance'] > min_dist) &
                            (contact_frame['target', 'distance'] < pistack_dist_max)]

    # If 'ringsel' is empty, then no contacting aromatic ring pairs
    if ringsel.empty:
        return contact_frame

    source_rings = []
    for source_residue in topology[topology['serial'].isin(ringsel[('source', 'serial')])].residues(extend=True):
        source_rings.extend(source_residue.find_rings(aromatic=True))

    target_rings = []
    for target_residue in topology[topology['serial'].isin(ringsel[('target', 'serial')])].residues(extend=True):
        target_rings.extend(target_residue.find_rings(aromatic=True))

    logger.info("Evaluate pi/T-stacking detection on {0} ring pairs using: pistack_dist_max={1}, pistack_ang_dev={2},"
                 "pistack_offset_max={3}".format(len(source_rings) * len(target_rings), pistack_dist_max,
                                                 pistack_ang_dev, pistack_offset_max))

    # Loop over ring pairs
    normalize_hf_labels = []
    for sring in source_rings:

        for tring in target_rings:

            stack, data = is_pistack(sring, tring, min_dist=min_dist, pistack_dist_max=pistack_dist_max,
                                     pistack_ang_dev=pistack_ang_dev, pistack_offset_max=pistack_offset_max)

            if stack:

                if data['type'] == 'ps':
                    normalize_hf_labels.extend(data['sring'])
                    normalize_hf_labels.extend(data['tring'])

                newindex = max(contact_frame.index) + 1
                contact_frame.loc[newindex, 'contact'] = data['type']
                contact_frame.loc[newindex, ('target', 'distance')] = data['distance']
                contact_frame.loc[newindex, ('target', 'angle')] = data['angle']
                for label in ['segmentID', 'chainID', 'resName', 'resSeq']:
                    contact_frame.loc[newindex, ('source', label)] = sring[label].unique()[0]
                contact_frame.loc[newindex, ('source', 'name')] = 'X1'
                contact_frame.loc[newindex, ('source', 'attype')] = 'Du'
                contact_frame.loc[newindex, ('source', 'serial')] = -1
                contact_frame.loc[newindex, ('source', 'element')] = 'D'
                for label in ['segmentID', 'chainID', 'resName', 'resSeq']:
                    contact_frame.loc[newindex, ('target', label)] = tring[label].unique()[0]
                contact_frame.loc[newindex, ('target', 'name')] = 'X1'
                contact_frame.loc[newindex, ('target', 'attype')] = 'Du'
                contact_frame.loc[newindex, ('target', 'serial')] = -1
                contact_frame.loc[newindex, ('target', 'element')] = 'D'

    # If pi-stacking, remove 'hf' label for these atoms
    if normalize_hf_labels:
        logger.info("Remove 'hf' label for {0} atoms in pi-stacked rings".format(len(normalize_hf_labels)))
        for idx, atom in contact_frame[(contact_frame[('source', 'serial')].isin(normalize_hf_labels)) &
                                       (contact_frame['contact'] == 'hf')].iterrows():
            contact_frame.loc[idx, 'contact'] = remove_contact_type(contact_frame.loc[idx, 'contact'], 'hf')

    # TODO: Reset dtype on atnum and resnum to int64 again. They get changed to float64 somehow.
    for n in ('source', 'target'):
        contact_frame[(n, 'serial')] = contact_frame[(n, 'serial')].astype('int64')
        contact_frame[(n, 'resSeq')] = contact_frame[(n, 'resSeq')].astype('int64')

    return contact_frame
