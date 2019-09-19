# -*- coding: utf-8 -*-

import copy
import logging

from interact import __module__
from interact.core.geometry import angle, distance
from interact.core.helpers import set_contact_type

logger = logging.getLogger(__module__)


def eval_hbonds(contacts, topology, max_hbond_dist=0.41, hbond_don_anglediv=50.0, hbond_acc_anglediv=90.0,
                optimize=True):
    """
    Evaluate the presence of hydrogen bonded contacts in the provided
    contact DataFrame. This function does not evaluate water bridges.

    Prerequisites:
    This function uses the SYBYL atom types to identify possible
    hydrogen bond donor and acceptor atoms. At least one covalently
    bonded hydrogen is expected for donors and subsequently the
    possible bonding geometry in terms of distances and angles is
    evaluated.

    As a result the function requires the input structure to be fully
    protonated or having at least polar hydrogens attached.
    The geometry of the attached hydrogens influences the contacts
    identified. If the structure is not, or partially protonated, the
    method used to add hydrogens will influence the identified contacts.
    Structure derived from moleculare dynamics will likely having their
    (polar) hydrogens oriented as such to reflect a hydrogen bond if
    present. If hydrogens are added with another program this may not
    be the case. OpenBabel for instance will add hydrogens in standard
    conformation not taking into account the environment of the atom
    to wich hydrogens are attached. The HBplus program (McDonald I K &
    Thornton J M (1994). J. Mol. Biol., 238, 777-793.) also part of
    the LIGPLOT program, will optimize local hydrogen geometry first.

    Differences in the geometry of added hydrogens will mostly affect
    the angle criteria rather than the distance. To correct for
    non-optimized H-atom geometry without the need for optimization,
    the function allows the angle criteria to be a function of the
    number of attached atoms using the 'optimize' option. Angles are
    then defined as:

    180 / number of non-isolated covalent neighbours - 1

    for all donors that are not trigonal planar (N.pl3, N.plc, N.ar,
    N.2, O.2, O.co2, S.a)

    Algorithm:
    1) Select all heavy atom contacts within max_hbond_dist.
    2) Identify donor-acceptor pairs for source to target and target to
       source based on SYBYL atom types (see below).
    3) Check if donor has at least one covalently bonded H-atom
    4) Check if angle between donor - H - acceptor does not deviate more
       than hbond_don_anglediv from it's ideal in-plane (180) degree
       orientation (cone fit), (Hubbard & Haider, 2001). The value for
       hbond_don_anglediv is either fixed or a function of the number
       of covalently attached atoms (see above) when 'optimize' is True.
    5) Check if the angle between the heavy atom acceptor neighbour -
       acceptor - H does not deviate more than hbond_acc_anglediv.
    6) Check that distance donor-acceptor heavy atom is larger than
       donor-H-acceptor.

    atom descriptor       base type   donor1  acceptor    directionality
    --------------------------------------------------------------------
    sp3 N                  N.3        y       y           along lone pair
    sp2 N                  N.2        y       y           along lone pair
    sp  N                  N.1        n       y           along lone pair
    Acidic N               N.acid     y       y           along lone pair 2
    Aromatic N             N.ar       y       y           along lone pair
    Amide N                N.am       y       n
    Quaternary N           N.4        y       n
    Uncharged trigonal N   N.pl3      y       n           3
    Charged trigonal N     N.plc      y       n           4
    Hydroxyl O             O.3        y       y           in plane of lone pair
    Ether O                O.3        n       y           in plane of lone pair
    Carboxylate O          O.co2      n       y           along lone pair
    Carbonyl O             O.2        n       y           in plane of lone pair
    Nitro O                O.2        n       y           along lone pair
    N-oxide O              O.2        n       y
    Amide O                O.2        n       y           in plane of lone pair
    Neutral sulfur-bound O O.2        n       y           5
    Charged sulfur-bound O O.co2      n       y           cone 6
    Phosphate O            O.co2      n       y           cone
    Borate O               O.co2      n       y           cone
    Other neg-charged O    O.co2      n       y
    Negative charged S     S.m        n       y           along lone pair
    sp2 S                  S.a        n       y           along lone pair

    1: Provided at least one H-atom covalently bound
    2: An acidic nitrogen is a nitrogen bound by at least two single bonds
    3: As in uncharged histidine residue
    4: As in a guanidino residue
    5: As in sulfonamides, sulfoxides, sulfones
    6: As in sulphate groups

    :param contacts:           contact DataFrame
    :type contacts:            :py:DataFrame
    :param topology:           Pandas DataFrame representing the structure
    :type topology:            :interact:ToplogyDataFrame
    :param max_hbond_dist:     Maximum hydrogen bond distance cutoff
    :type max_hbond_dist:      :py:float
    :param hbond_don_anglediv: Maximum hydrogen bond donor-H-acceptor
                               angle deviation.
    :type hbond_don_anglediv:  :py:float
    :param hbond_acc_anglediv: Maximum hydrogen bond acceptor'-acceptor-H
                               angle deviation.
    :type hbond_acc_anglediv:  :py:float
    :param optimize:           Rather to optimize angle cutoff based on
                               donor atom geometry.
    :type optimize:            :py:bool

    :return:                   Changes the 'contact' label in the contacts to
                               hb-ad (hydrogen bond acceptor-donor) or hb-da
                               (hydrogen bond donor-acceptor) for identified
                               hydrogen-bonded contacts. Also add the value of
                               the donor-H-acceptor angle.
    :rtype:                    :pandas:DataFrame
    """

    # Preselect all contacts below max_hbond_dist
    hbdist = contacts[(contacts['target', 'distance'] <= max_hbond_dist)]

    logger.info("Init eval_hbonds with {0} contacts using: max_hbond_dist={1}".format(len(hbdist), max_hbond_dist))

    # Query for potential hbond donor-acceptor pairs
    accpt_attypes = ('N.3', 'N.2', 'N.1', 'N.acid', 'N.ar', 'O.3', 'O.co2', 'O.2', 'S.m', 'S.a')
    donor_attypes = ('N.3', 'N.2', 'N.acid', 'N.am', 'N.ar', 'N.4', 'N.pl3', 'N.plc', 'O.3')
    donor_avoid = ('N.pl3', 'N.plc', 'N.ar', 'N.2', 'O.2', 'O.co2', 'S.a')

    # Define donor_acceptor pairs. Source donor - target acceptor and vice versa
    donor_acceptor_dict = dict()
    donor_acceptor_dict['source'] = hbdist[(hbdist['source', 'attype'].isin(donor_attypes)) &
                                           (hbdist['target', 'attype'].isin(accpt_attypes))]
    donor_acceptor_dict['target'] = hbdist[(hbdist['source', 'attype'].isin(accpt_attypes)) &
                                           (hbdist['target', 'attype'].isin(donor_attypes))]

    logger.info("{0} contacts after selecting for donor-acceptor pairs".format(len(donor_acceptor_dict['source'])))
    logger.info("{0} contacts after selecting for acceptor-donor pairs".format(len(donor_acceptor_dict['target'])))

    # Search for hbonds
    label_dict = {'source': 'hb-da', 'target': 'hb-ad'}
    anglediv = copy.copy(hbond_don_anglediv)
    for direction, selection in donor_acceptor_dict.items():

        target = 'source'
        if direction == 'source':
            target = 'target'

        # Ensure donor and acceptors have neighbours (e.a. not ions etc.)
        for idx, n in selection.iterrows():
            donor = topology[topology['serial'] == n[direction, 'serial']]
            acceptor = topology[topology['serial'] == n[target, 'serial']]

            donor_bonded = donor.neighbours(covalent=True)
            acceptor_bonded = acceptor.neighbours(covalent=True)
            acceptor_bonded = acceptor_bonded[(acceptor_bonded['attype'] != 'H') &
                                              (acceptor_bonded['resSeq'] == n[target, 'resSeq'])]

            # Get donor and acceptor heavy-atom coordinates
            donor_heavy_coor = donor.coord
            acceptor_heavy_coor = acceptor.coord

            # There should at least be covalent neighbours (e.a. not ions etc.)
            if donor_bonded.empty or acceptor_bonded.empty:
                logger.debug('No neighbours in contact pair {0}-{1}, skipping'.format(n[direction, 'serial'],
                                                                                      n[target, 'serial']))
                continue

            # Check if there are H-atoms attached and asses H-bond geometry criteria
            for idy, h in donor_bonded[(donor_bonded['attype'] == 'H') &
                                       (donor_bonded['resSeq'] == n[direction, 'resSeq'])].iterrows():

                # Get donor hydrogen coordinates
                donor_h_coor = h.coord

                # Angle donor - H - acceptor
                angle1 = angle(donor_heavy_coor, donor_h_coor, acceptor_heavy_coor)
                dist1 = distance(donor_h_coor, acceptor_heavy_coor)

                # Angle acceptor_neigh - acceptor - H, keep largest.
                angle2 = []
                for serial in acceptor_bonded['serial']:
                    acceptor_neigh = topology[topology['serial'] == serial].coord
                    angle2.append(angle(acceptor_neigh, acceptor_heavy_coor, donor_h_coor))

                # If optimize equals True, determine donor-H-acceptor angle deviation based on covalent bonding
                # geometry for all non trigonal planar donor atoms
                hbond_don_anglediv = anglediv
                if optimize and not donor['attype'].values[0] in donor_avoid:
                    substitutions = 0
                    for idz, i in donor_bonded.iterrows():
                        r = topology[topology['serial'] == i['serial']].neighbours(covalent=True)
                        if len(r) > 1:
                            substitutions += 1
                    try:
                        hbond_don_anglediv = (180 / float(substitutions))
                    except ZeroDivisionError:
                        hbond_don_anglediv = 0.0

                if (180 - hbond_don_anglediv < abs(angle1) < 180 + hbond_don_anglediv) and \
                   (180 - hbond_acc_anglediv < abs(max(angle2)) < 180 + hbond_acc_anglediv) and \
                   (contacts.loc[idx, 'target'].distance > dist1):

                    contacts.loc[idx, 'contact'] = set_contact_type(contacts.loc[idx, 'contact'], label_dict[direction])
                    contacts.loc[idx, ('target', 'angle')] = angle1
                    logger.info(
                        "H-bond between {0}-{1} {2}-{3} and {4}-{5} {6}-{7}. Distance D-A: {8:.3f}, "
                        "Distance DH-A: {9:.3f}, angle: {10:.2f} deg. hbond_don_anglediv: {11:.2f}".format(
                            contacts.loc[idx, 'source'].resSeq, contacts.loc[idx, 'source'].resName,
                            contacts.loc[idx, 'source'].serial,
                            contacts.loc[idx, 'source'].name, contacts.loc[idx, 'target'].resSeq,
                            contacts.loc[idx, 'target'].resName,
                            contacts.loc[idx, 'target'].serial, contacts.loc[idx, 'target'].name,
                            contacts.loc[idx, 'target'].distance, dist1, angle1, hbond_don_anglediv
                        ))

    return contacts


def eval_water_bridges(contacts, topology, min_wbridge_dist=0.25, max_wbridge_dist=0.40, min_omega_angle=75.0,
                       max_omega_angle=140.0, min_theta_angle=100.0, wbfilter=True):
    """
    Evaluate the presence of water mediated hydrogen bonded bridges
    in the provided contact DataFrame.

    Algorithm:
    1) Select all water oxygen atoms within the range defined by
       min_wbridge_dist (Jiang et al., 2005) - 0.01 nm and max_wbridge_dist
       (Jiang et al., 2005) + 0.04 nm
    2) For each water, get neighbouring atoms below max_wbridge_dist excluding
       other waters.
    3) In two loops look for ligand donor - water - other acceptor pairs and
       ligand acceptor - water - other donor pairs.
    4) For each pair check if there is at least one covalently bound hydrogen
       attached to the donor.
    5) Check the theta angle (water O - donor H - donor), should be larger than
       min_theta_angle (Jiang et al., 2005).
    6) Check the omega angle (acceptor - water O - donor H), should be in the
       range defined by min_omega_angle, max_omega_angle (Jiang et al., 2005).
    7) If wbfilter option is True: a water molecule is only allowed to
       participate as donor in two hydrogen bonds (two hydrogen atoms as
       donors). In the case of more than two possible hydrogen bonds for a
       water molecule as donor, only the two contacts with a water angle
       closest to 110 deg. and/or smaller H-bond distances are kept.

    :param contacts:           contact DataFrame
    :type contacts:            :py:DataFrame
    :param topology:           Pandas DataFrame representing the structure
    :type topology:            :interact:ToplogyDataFrame
    :param min_wbridge_dist:   Minimal distance for water bridged hydrogen bonds
    :type min_wbridge_dist:    :py:float
    :param min_wbridge_dist:   Maximum distance for water bridged hydrogen bonds
    :type min_wbridge_dist:    :py:float
    :param min_omega_angle:    Minimum omega angle, acceptor-waterO-donorH
    :type min_omega_angle:     :py:float
    :param max_omega_angle:    Maximum omega angle, acceptor-waterO-donorH
    :type max_omega_angle:     :py:float
    :param min_theta_angle:    Minimum theta angle, waterO-donorH-donor

    :return:                   Changes the 'contact' label in the contacts to
                               hb-ad (hydrogen bond acceptor-donor) or hb-da
                               (hydrogen bond donor-acceptor) for identified
                               hydrogen-bonded contacts. Also add the value of
                               the donor-H-acceptor angle.
    :rtype:                    :pandas:DataFrame

    """

    # Preselect all water oxygen's close to ligand
    wbdist = contacts[(contacts['target', 'distance'] > min_wbridge_dist) &
                      (contacts['target', 'distance'] <= max_wbridge_dist) &
                      (contacts['target', 'resName'] == 'HOH') & (contacts['target', 'attype'] == 'O.3')]

    if wbdist.empty:
        logger.debug('No water oxygen atoms detected close to the ligand')
        return contacts

    logger.info("Run water bridge detection on {0} possible contacts using: min_wbridge_dist={1},"
                "max_wbridge_dist={2}, min_omega_angle={3}, max_omega_angle={4}, min_theta_angle={5},"
                "wbfilter={6}".format(len(wbdist), min_wbridge_dist, max_wbridge_dist, min_omega_angle,
                                      max_omega_angle, min_theta_angle, wbfilter))

    # Query for potential hbond donor-acceptor pairs
    accpt_attypes = ('N.3', 'N.2', 'N.1', 'N.acid', 'N.ar', 'O.3', 'O.co2', 'O.2', 'S.m', 'S.a', 'F', 'Br', 'Cl')
    donor_attypes = ('N.3', 'N.2', 'N.acid', 'N.am', 'N.4', 'N.pl3', 'N.plc', 'O.3')

    # Loop over waters looking for water bridges
    ligresnum = wbdist['source', 'resSeq'].unique()
    for water in sorted(wbdist['target', 'serial'].unique()):

        # Get neighbouring atoms
        water = topology.loc[(topology['serial'] == water)]
        water_neigh = water.neighbours(cutoff=max_wbridge_dist)
        w = water.coord

        # Remove other waters
        if not water_neigh.empty:
            water_neigh = water_neigh[(water_neigh['resName'] != 'HOH')]

        # Query possible ligand donor - water - acceptor contacts
        dwa_pairs = []
        for idd, d in water_neigh[(water_neigh['resSeq'].isin(ligresnum)) &
                                  (water_neigh['attype'].isin(donor_attypes))].iterrows():

            donor = topology[topology['serial'] == d['serial']]
            covalent_neighbours = donor.neighbours(covalent=True)
            x = donor.coord

            # Check if there are H-atoms attached and asses H-bond geometry criteria
            for idy, h in covalent_neighbours[covalent_neighbours['attype'] == 'H'].iterrows():
                y = topology[topology['serial'] == h['serial']].coord

                # Check theta angle: donor - hydrogen - water oxygen
                theta = angle(x, y, w)
                if abs(theta) > min_theta_angle:

                    # Loop over possible acceptors
                    for ida, a in water_neigh[~(water_neigh['resSeq'].isin(ligresnum)) &
                                              (water_neigh['attype'].isin(accpt_attypes))].iterrows():

                        acceptor = topology[topology['serial'] == a['serial']].coord

                        # Check omega angle: acceptor - water oxygen - donor h
                        omega = 180 - angle(acceptor, w, y)
                        if min_omega_angle < abs(omega) < max_omega_angle:
                            dist_aw = distance(w, acceptor)
                            dist_wd = distance(w, x)
                            dwa_pairs.append(
                                (omega, theta, dist_aw, dist_wd, d['serial'], a['serial'], water.serial.values[0]))

                            logger.info(
                                "Water bridge: donor {0}-{1}-{2}, acceptor {3}-{4}-{5}, water {6}. "
                                "Dist d-w {7:.3f} a-w {8:.3f}. Omega: {8:.2f} Theta {9:.2f}".format(
                                    donor.resName.values[0], donor.resSeq.values[0], donor.name.values[0],
                                    a.resName, a.resSeq, a.name, water.resSeq.values[0], dist_wd, dist_aw,
                                    omega, theta))

        if wbfilter and len(dwa_pairs) > 1:
            dwa_pairs.sort(key=lambda v: (110 - v[0]) + (v[2] + v[3]))
            dwa_pairs = [dwa_pairs[0]]

        for bridge in dwa_pairs:
            cid = contacts[(contacts['source', 'serial'] == bridge[4]) & (contacts['target', 'serial'] == bridge[6])]
            tid = topology[topology['serial'] == bridge[5]]
            contacts.loc[cid.index, 'contact'] = set_contact_type(contacts.loc[cid.index, 'contact'], 'wb-da')
            contacts.loc[cid.index, ('target', 'angle')] = bridge[1]

            newindex = contacts.index.max() + 1
            for mdx in [col for col in contacts['target'].columns if not col == 'index']:
                contacts.loc[newindex, ('target', mdx)] = cid['target', mdx].values[0]
            for mdx in [col for col in contacts['source'].columns if not col == 'index']:
                contacts.loc[newindex, ('source', mdx)] = tid[mdx].values[0]

            contacts.loc[newindex, 'contact'] = 'wb-da'
            contacts.loc[newindex, ('target', 'distance')] = bridge[2]
            contacts.loc[newindex, ('target', 'angle')] = bridge[0]

        # Query possible ligand acceptor - water - donor contacts
        awd_pairs = []
        for ida, a in water_neigh[(water_neigh['resSeq'].isin(ligresnum)) &
                                  (water_neigh['attype'].isin(accpt_attypes))].iterrows():

            acceptor = topology[topology['serial'] == a['serial']].coord
            # Loop over possible donors
            for idd, d in water_neigh[~(water_neigh['resSeq'].isin(ligresnum)) &
                                      (water_neigh['attype'].isin(donor_attypes))].iterrows():

                donor = topology[topology['serial'] == d['serial']]
                covalent_neighbours = donor.neighbours(covalent=True)
                x = donor.coord

                # Check if there are H-atoms attached and asses H-bond geometry criteria
                for idy, h in covalent_neighbours[covalent_neighbours['attype'] == 'H'].iterrows():
                    y = topology[topology['serial'] == h['serial']].coord

                    # Check theta angle: donor - hydrogen - water oxygen
                    theta = angle(x, y, w)
                    if abs(theta) > min_theta_angle:

                        # Check omega angle: acceptor - water oxygen - donor h
                        omega = 180 - angle(acceptor, w, y)
                        if min_omega_angle < abs(omega) < max_omega_angle:
                            dist_aw = distance(w, acceptor)
                            dist_wd = distance(w, x)
                            awd_pairs.append(
                                (omega, theta, dist_aw, dist_wd, d['serial'], a['serial'], water.serial.values[0]))

                            logger.info(
                                "Water bridge: donor {0}-{1}-{2}, acceptor {3}-{4}-{5}, water {6}. "
                                "Dist d-w {7:.3f} a-w {8:.3f}. Omega: {8:.2f} Theta {9:.2f}".format(
                                    donor.resName.values[0], donor.resSeq.values[0], donor.name.values[0],
                                    a.resName, a.resSeq, a.name, water.resSeq.values[0], dist_wd, dist_aw,
                                    omega, theta))

        if wbfilter and len(awd_pairs) > 2:
            awd_pairs.sort(key=lambda x: (110 - x[0]) + (x[2] + x[3]))
            awd_pairs = awd_pairs[:2]

        for bridge in awd_pairs:
            cid = contacts[(contacts['source', 'serial'] == bridge[5]) & (contacts['target', 'serial'] == bridge[6])]
            tid = topology[topology['serial'] == bridge[4]]
            contacts.loc[cid.index, 'contact'] = set_contact_type(contacts.loc[cid.index, 'contact'], 'wb-ad')
            contacts.loc[cid.index, ('target', 'angle')] = bridge[1]

            newindex = contacts.index.max() + 1
            for mdx in [col for col in contacts['target'].columns if not col == 'index']:
                contacts.loc[newindex, ('target', mdx)] = cid['target', mdx].values[0]
            for mdx in [col for col in contacts['source'].columns if not col == 'index']:
                contacts.loc[newindex, ('source', mdx)] = tid[mdx].values[0]

            contacts.loc[newindex, 'contact'] = 'wb-ad'
            contacts.loc[newindex, ('target', 'distance')] = bridge[2]
            contacts.loc[newindex, ('target', 'angle')] = bridge[0]

    # TODO: Reset dtype on atnum and resnum to int64 again. They get changed to float64 somehow.
    for n in ('source', 'target'):
        contacts[(n, 'serial')] = contacts[(n, 'serial')].astype('int64')
        contacts[(n, 'resSeq')] = contacts[(n, 'resSeq')].astype('int64')

    return contacts
