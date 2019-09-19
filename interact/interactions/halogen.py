# -*- coding: utf-8 -*-

import logging

from interact import __module__
from interact.core.geometry import angle
from interact.core.helpers import set_contact_type

logger = logging.getLogger(__module__)


def eval_halogen_bonds(contact_frame, topology, min_dist=0.05, halogen_max_dist=0.41, halogen_don_angle=165.0,
                       halogen_acc_angle=120.0, halogen_angle_dev=30.0, halogens=('I', 'Br', 'Cl', 'F')):
    """
    Evaluate halogen interactions between source and target selections.

    Evalute halogen bonds according to the following rule set extracted
    from Auffinger (2004):
    1) Halogen bond donor should be an element of type `halogens` having
       only a single carbon neighbour.
    2) Halogen bond acceptor should be an carbonyl, hydroxyl, charged
       carboxylate, or phosphate group evaluated as an O,N,S element 
       with one single atom neighbour of element type P,C or S.
    3) Distance between donor and acceptor should be within `min_dist`
       and `halogen_max_dist`.
    4) The donor angle (C-X -- O) should not deviate more then
       `halogen_angle_dev` from `halogen_don_angle`
    5) The acceptor angle (Y-O -- X) should not deviate more then
       `halogen_angle_dev` from `halogen_acc_angle`

    Detected halogen bond are assigned contact type 'xb' in the returned
    contacts DataFrame.

    Reference: P. Auffinger, Halogen bonds in biological molecules (2004),
               PNAS: vol. 101 no. 38. vol. 16789-16794

    :param contact_frame:     contacts DataFrame
    :type contact_frame:      :pandas:DataFrame
    :param topology:          main system topology
    :type topology:           :interact:TopologyDataFrame
    :param halogen_max_dist:  maximum distance for halogen bonds (nm)
    :type halogen_max_dist:   :py:float
    :param halogen_don_angle: average donor angle (C-X -- O)
    :type halogen_don_angle:  :py:float
    :param halogen_acc_angle: average acceptor angle (Y-O -- X)
    :type halogen_acc_angle:  :py:float
    :param halogen_angle_dev: maximum angle deviation for donor and
                              acceptor angles
    :type halogen_angle_dev:  :py:float
    :param min_dist:          minimum interaction distance (nm)
    :type min_dist:           :py:float
    :param halogens:          Halogen element types
    :type halogens:           :py:tuple

    :return:                  contacts DataFrame
    :rtype:                   :pandas:DataFrame
    """

    # Evaluate source - target and target - source halogen bonds
    for pair in (('source', 'target'), ('target', 'source')):

        # Preselect all contacts between source halogen and target oxygen,nitrogen
        # or sulfur below halogen_max_dist not involving waters
        hadist = contact_frame[(contact_frame[pair[0], 'element'].isin(halogens)) &
                               (contact_frame[pair[1], 'element'].isin(('O', 'N', 'S'))) &
                               (contact_frame['target', 'distance'] > min_dist) &
                               (contact_frame['target', 'distance'] <= halogen_max_dist) &
                               (contact_frame['target', 'resName'] != 'HOH')]

        logging.info(
            "Evaluate {0} potential {1} halogen bonds using: halogen_max_dist={2:.3f}, halogen_don_angle={3:.2f},"
            "halogen_acc_angle={4:.2f}, halogen_angle_dev={5:.2f}".format(
                len(hadist), pair, halogen_max_dist, halogen_don_angle, halogen_acc_angle, halogen_angle_dev))

        if hadist.empty:
            return contact_frame

        for idx, n in hadist.iterrows():

            source = topology[topology['serial'] == n[pair[0], 'serial']]
            source_neighbours = source.neighbours(covalent=True, cutoff=0.185)
            target = topology[topology['serial'] == n[pair[1], 'serial']]
            target_neighbours = target.neighbours(covalent=True)

            # Ensure source (halogen) only has a single carbon neighbour
            c = source_neighbours[source_neighbours['element'] == 'C']
            if len(c) > 1:
                continue

            # Ensure the target (oxygen, nitrogen, sulfur) only has C,P,S as single neighbour
            y = target_neighbours[target_neighbours['element'].isin(('P', 'C', 'S'))]
            if len(y) > 1:
                continue

            # Calculate donor (C-X -- O) and acceptor (Y-O -- X) angles
            donor_angle = angle(c.coord, source.coord, target.coord)
            acceptor_angle = angle(y.coord, target.coord, source.coord)
            if ((halogen_don_angle - halogen_angle_dev) < abs(donor_angle) < (halogen_don_angle + halogen_angle_dev)) and (
                    (halogen_acc_angle - halogen_angle_dev) < abs(acceptor_angle) < (halogen_acc_angle + halogen_angle_dev)):
                contact_frame.loc[idx, 'contact'] = set_contact_type(contact_frame.loc[idx, 'contact'], 'xb')
                contact_frame.loc[idx, ('target', 'angle')] = donor_angle
                logging.info(
                    "Halogen bond between {0}-{1} {2}-{3} and {4}-{5} {6}-{7}. Distance D-A: {8:.3f}A, donor angle: {9:.2f}"
                    "deg. acceptor angle: {10:.2f}".format(source.resSeq.values[0], source.resName.values[0],
                        source.serial.values[0], source.name.values[0], target.resSeq.values[0], target.resName.values[0],
                        target.serial.values[0], target.name.values[0], contact_frame.loc[idx, 'target'].distance,
                        donor_angle, acceptor_angle))

    return contact_frame
