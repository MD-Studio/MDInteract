# -*- coding: utf-8 -*-

import logging
import numpy
import itertools

from interact import __module__
from interact.core.helpers import set_contact_type
from interact.core.geometry import plane_fit, distance, angle, projection, vector_angle
from interact.interactions.utils import is_pication

logger = logging.getLogger(__module__)


def eval_saltbridge(contact_frame, topology, max_charge_dist=0.55):
    """
    Evaluate contacts between centers of positive and negative charge.
    Physiological relevant pH is assumed.

    Algorithm:
    1) Primary selection is between all source and target atoms that are
       max_charge_dist apart according to (Barlow and Thornton, 1983) +
       0.15 nm
    2) Select all residues in previous selection that have a formal
       positive or negative charge according to the sum of partial charges
       in the 'charge' column. The latter charges are Gasteiger partial
       charges by default.
    3) Select all atoms that are likely a part of the charged group in
       in the residues from step 2 based on SYBYL atom types following:

    amino-acid                  type   atom         charge
    --------------------------------------------------------------------
    Arginine - Arg - R          N.pl3  RNHC(NH2)2+  +
    Lysine - Lys - K            N.4    RNH3         +
    Histidine - His - H         N.ar   ND1, NE2     +
    Aspartic acid - Asp - D     O.co2  RCOO-        -
    Glutamic acid - Glu - E     O.co2  RCOO-        -

    Ligands                     type   atom         charge
    --------------------------------------------------------------------
    quaterny ammonium           N.4                 +
    tertiary amines             N.am                +
    sulfonium groups            S.3                 +
    guanidine groups            C.cat               +
    phosphate                   O.co2  PO4          -
    sulfonate                   S.3    RSO2O-       -
    sulfonic acid               S.O2                -
    carboxylate                 O.co2               -

    4) Select neighbouring atoms not of element type 'C' or 'H' in the
       selection from step 3 to define all atoms part of the charged group
    5) Evaluate salt-bridges by the distance between the geometrical centers
       of two charged groups of opposite sign being smaller or equal to
       max_charge_dist.

    Although multiple atoms of both charged groups take part in the salt-bridge
    only the pair with the shortest atom-atom distance is reported using the
    labels:

    - 'sb-pn': for a positive charged source to negative target contact.
    - 'sb-np': for a negative charged source to positive target contact.

    Because salt-bridges are composed of hydrogen-bonded and charged
    interactions, the reported atom pair often was reported before as taking
    part in a hydrogen-bonded interactions when the `eval_hbonds` function was
    used. The salt-bridge label will be added to the contact column maintaining
    the hydrogen bond label.

    :param contact_frame:      contact DataFrame
    :type contact_frame:       :py:DataFrame
    :param topology:           Pandas DataFrame representing the structure
    :type topology:            :interact:ToplogyDataFrame
    :param max_charge_dist:    maximum distance cutoff between charge centers
    :type max_charge_dist:     :py:float

    :return:                   Adds the labels 'sb-np' or 'sb-pn' to the
                               'contact' column of the input contact frame.
    :rtype:                    :pandas:DataFrame
    """

    # Preselect all contacts below max_charge_dist
    chdist = contact_frame[contact_frame['target', 'distance'] <= max_charge_dist]

    # Select all charged source and target residues
    charged_groups = {'source': [], 'target': []}
    for group in charged_groups.keys():
        for charge_group in topology[topology['serial'].isin(chdist[group, 'serial'])].find_charged_centers():
            if charge_group[1] <= -1:
                charged_groups[group].append(('n', charge_group[0]))
            else:
                charged_groups[group].append(('p', charge_group[0]))

    if not len(charged_groups['source']) or not len(charged_groups['target']):
        logger.info('Not running salt-bridge detection. Charged groups in source: {0}, target: {1}'.format(
            len(charged_groups['source']), len(charged_groups['target'])))
        return contact_frame

    logger.info(
        "Run salt-bridge detection on {0} source and {1} target charged groups using: max_charge_dist={2}".format(
            len(charged_groups['source']), len(charged_groups['target']), max_charge_dist))

    # Loop over combinations of source and target charged groups
    for s, l in itertools.product(charged_groups['source'], charged_groups['target']):
        if s[0] != l[0]:
            center_distance = distance(s[1].center(), l[1].center())
            sb_type = 'sb-{0}{1}'.format(s[0], l[0])
            source_center = repr(list(s[1]['serial'])).strip('[]')
            target_center = repr(list(l[1]['serial'])).strip('[]')
            if center_distance <= max_charge_dist:
                logger.info('{0} between {1}-{2} and {3}-{4}. D: {5:.3f} nm between groups {6} and {7}'.format(
                    sb_type, s[1]['resSeq'].unique()[0], s[1]['resName'].unique()[0], l[1]['resSeq'].unique()[0],
                    l[1]['resName'].unique()[0], center_distance, source_center, target_center))

            # Report salt-bridges
            sb_selection = contact_frame[(contact_frame['source', 'serial'].isin(s[1]['serial'])) &
                                          contact_frame['target', 'serial'].isin(l[1]['serial'])]
            report_to = sb_selection.sort_values(by=('target', 'distance')).head(n=1)
            contact_frame.loc[report_to.index, 'contact'] = set_contact_type(
                contact_frame.loc[report_to.index, 'contact'], sb_type)

    return contact_frame


def eval_heme_coordination(contact_frame, topology, rings=None, heme_dist_prefilter=0.55, heme_dist_max=0.35,
                           heme_dist_min=0, min_heme_coor_angle=105, max_heme_coor_angle=160, fe_ox_dist=0.16,
                           exclude=('H', 'O.3', 'O.2', 'O.co2', 'O.spc', 'O.t3p', 'C.cat', 'S.o2')):
    """
    Evaluate heme coordination of ligand atoms
    """

    rings = rings or []

    # Select all atoms within heme_dist_prefilter distance from Fe excluding atoms in exclude list
    fedist = contact_frame[(contact_frame['target', 'name'] == 'FE') &
                           (~contact_frame['source', 'attype'].isin(exclude)) &
                           (contact_frame['target', 'distance'] < heme_dist_prefilter)]
    if fedist.empty:
        return contact_frame

    # Get Fe atom
    fe = topology[(topology['resName'] == 'HEM') & (topology['name'] == 'FE')]
    if fe.empty:
        logger.warn("Unable to asses heme coordination. Fe atom not found")
        return contact_frame

    # Get four nitrogen atoms coordinating the Fe atom
    fe_neigh = fe.neighbours(cutoff=0.3)
    fe_coordinating = fe_neigh[(fe_neigh['resName'] == 'HEM') & (fe_neigh['element'] == 'N')].sort_values(by='name')
    if len(fe_coordinating) != 4:
        logger.warn("Unable to asses heme coordination. Found {0} nitrogen atoms coordinating Fe. Expected 4".format(
            len(fe_coordinating)))
        return contact_frame

    logger.debug("Run heme coordination detection on {0} possible contacts using: heme_dist_prefilter={1:.2f}, "
                 "heme_dist_min={2:.2f}, heme_dist_max={3:.2f}, min_heme_coor_angle={4:.2f}, max_heme_coor_angle={5:.2f}, "
                 "fe_ox_dist={6:.2f}".format(fedist.shape[0], heme_dist_prefilter, heme_dist_min, heme_dist_max,
                                             min_heme_coor_angle, max_heme_coor_angle, fe_ox_dist))

    # Calculate normals between Nitrogens -> Fe vectors
    fe_coor = fe.coord
    n_coor = fe_coordinating.coord - fe_coor

    m1 = numpy.cross(n_coor[0], n_coor[1])
    m2 = numpy.cross(n_coor[1], n_coor[2])
    m3 = numpy.cross(n_coor[2], n_coor[3])
    m4 = numpy.cross(n_coor[3], n_coor[0])

    # Is there an Oxygen above the heme (complex I) or do we need to place a dummy
    close_fe_neigh = fe.neighbours(cutoff=0.2)
    dummyox = close_fe_neigh[(close_fe_neigh['resName'] == 'HEM') & (close_fe_neigh['element'] == 'O')]
    mv = numpy.mean(numpy.vstack((m1, m2, m3, m4)), axis=0)
    if len(dummyox) == 1:
        dummyox = dummyox.coord
        logger.info('Oxygen atom bonded to Fe (complex I)')

    else:
        # Calculate dummy O atom from the average of the four normals
        # Normalize normal mean, change vector size to 1.6 A and set point
        dummyox = ((mv / numpy.linalg.norm(mv)) * fe_ox_dist) + fe_coor
        logger.info("Reconstructed oxygen atom placed {0}nm above Heme Fe at position {1}".format(fe_ox_dist, ' '.join(
            ['{0:.3f}'.format(c) for c in dummyox])))

    # Check the coordination of the Fe atom by the SG atom of the Cys below Heme
    sg = fe_neigh[(fe_neigh['resName'] == 'CYS') & (fe_neigh['name'] == 'SG')]
    if not sg.empty:
        sg_angle = angle(dummyox, fe_coor, sg.coord)
        if not 160 < sg_angle < 200:
            logger.warn("Angle between reconstructed oxygen -> Fe -> Cys SG has unusual value {0:.3f}".format(sg_angle))
    else:
        logger.warn("No CYS SG atom in a distance of 0.3nm of the Heme Fe atom")

    # Check if there are rings with there center of mass below heme_dist_prefilter from heme FE.
    # Calculate ring normals
    ring_normals = []
    for aromatic in rings:
        aromatic_center = aromatic.center()
        aromatic_fe_dist = distance(fe_coor, aromatic_center)
        if aromatic_fe_dist < heme_dist_prefilter:
            aromatic_norm = plane_fit(aromatic.coord, center=aromatic_center)
            aromatic_norm_angle = vector_angle(aromatic_norm, mv, deg=True)
            aromatic_norm_angle = min(aromatic_norm_angle,
                                      180 - aromatic_norm_angle if not 180 - aromatic_norm_angle < 0 else
                                      aromatic_norm_angle)

            ring = aromatic.index.tolist()
            ring_normals.append((aromatic_center, aromatic_norm, aromatic_norm_angle, ring))

            logger.info("Ring {0} close to heme Fe: distance center-Fe {1:.2f}nm, normal angle heme plane-ring:"
                        "{2:.2f} deg.".format(ring, aromatic_fe_dist, aromatic_norm_angle))

    # Get ligand atoms coordinated
    for idx, n in fedist.iterrows():

        source = topology[topology.index == n['source', 'index']]
        source_atom_type = n['source', 'attype']
        z = source.coord

        # Check for heme coordination by aromatic nitrogens. label as 'hc'
        if source_atom_type in ('N.ar', 'N.2', 'N.3'):
            ar_norm_angle = 90
            for ring in ring_normals:
                if n['source', 'index'] in ring[-1]:
                    ar_norm_angle = ring[2]
                    break
            fe_dist = distance(z, fe_coor)
            fe_offset = distance(projection(mv, fe_coor, z), fe_coor)
            if 45 < ar_norm_angle < 95 and fe_dist < 0.35 and fe_offset < 0.1:
                contact_frame.loc[idx, 'contact'] = set_contact_type(contact_frame.loc[idx, 'contact'], 'hc')
                contact_frame.loc[idx, ('target', 'angle')] = ar_norm_angle
                logger.info(
                    "Heme Fe coordination with {0} {1}. Distance: {2:.2f} A. offset: {3:.2f} A plane normal angle: {4:.2f}".format(
                        n['source', 'serial'],
                        n['source', 'name'], fe_dist, fe_offset, ar_norm_angle))

        # Check for possible sites of metabolism and label as 'hm'.
        # Filter on covalent neighbours and apply knowledge based rules.
        if source_atom_type in ('C.2', 'C.3', 'C.ar', 'N.1', 'N.2', 'N.4', 'N.pl3', 'S.3'):
            cutoff = 0.16
            if source_atom_type == 'S.3': cutoff = 0.18
            neigh = source.neighbours(cutoff=cutoff)
            neigh_atom_types = set(neigh['attype'])

            # If ligand atom is of type C.3 or C.ar it should contain at least one covalently bonded atom
            # of type ['H','Cl','I','Br','F','Hal']
            if source_atom_type in ('C.3', 'C.ar') and len(
                    neigh_atom_types.intersection({'H', 'Cl', 'I', 'Br', 'F', 'Hal'})) == 0:
                logger.debug(
                    "Ligand target atom {0}-{1} excluded. Atom type {2} not covalently bonded to: H,Cl,I,Br,F or Hal".format(
                        n['source', 'serial'], n['source', 'name'], source_atom_type))
                continue

            # If ligand atom is of type N.4 it should contain at least one covalently bonded atom of type H
            if source_atom_type == 'N.4' and not 'H' in neigh_atom_types:
                logger.debug(
                    "Ligand target atom {0}-{1} excluded. Atom type N.4 not covalently bonded to hydrogen".format(
                        n['source', 'serial'], n['source', 'name']))
                continue

            # Exclude carbons that are a part of ketone or carboxylate
            if source_atom_type == 'C.2' and 'O.2' in neigh_atom_types or 'O.co2' in neigh_atom_types:
                logger.debug(
                    "Ligand target atom {0}-{1} excluded. Atom type C.2 part of ketone or carboxylate group.".format(
                        n['source', 'serial'], n['source', 'name']))
                continue

            # Additional check on S.O2 wrongly labeled as S.3 (PLANTS?)
            if source_atom_type == 'S.3' and neigh[neigh['attype'] == 'O.2'].shape[0] == 2:
                logger.debug(
                    "Ligand target atom {0}-{1} excluded. Atom labeled as S.3 but probably S.O2 as it is covalently bonded to two O.2".format(
                        n['source', 'serial'], n['source', 'name']))
                continue

            # If N.pl3 or N.2 check for nitro- or nitrate group.
            if source_atom_type in ('N.pl3', 'N.2') and (
                    'O.co2' in neigh_atom_types or len(neigh_atom_types.intersection(set(['O.2', 'O.3']))) == 2):
                logger.debug("Ligand target atom {0}-{1} excluded. Atom type {2}, exclude nitro- nitrate group".format(
                    n['source', 'serial'], n['source', 'name'], source_atom_type))
                continue

            # Exclude (iso)-nitrile group
            if source_atom_type in 'N.1' and 'C.1' in neigh_atom_types:
                logger.debug(
                    "Ligand target atom {0}-{1} excluded. Carbon with Sp hybridized N".format(n['source', 'serial'],
                                                                                              n['source', 'name']))
                continue

            # Check Heme-Nitrogen coordination (Type II binding).
            if source_atom_type in ('C.ar', 'N.ar'):
                ar_norm_angle = None
                for ring in ring_normals:
                    if n['source', 'serial'] in ring[-1]:
                        ar_norm_angle = ring[2]
                        break
                if ar_norm_angle and not (45 < ar_norm_angle < 85 or 95 < ar_norm_angle < 135):
                    logger.debug(
                        "Ligand target atom {0}-{1} excluded. Aromatic C or N part of ring with angle of {2:.2f} with respect to Heme plane".format(
                            n['source', 'serial'], n['source', 'name'], ar_norm_angle))
                    continue

        fe_ox_angle = angle(fe_coor, dummyox, z[0])
        dist = distance(dummyox, z)
        if min_heme_coor_angle < abs(fe_ox_angle) < max_heme_coor_angle and heme_dist_min < dist < heme_dist_max:
            contact_frame.loc[idx, 'contact'] = set_contact_type(contact_frame.loc[idx, 'contact'], 'hm')
            contact_frame.loc[idx, ('target', 'angle')] = fe_ox_angle
            logger.info(
                "Heme Fe possible som with {0} {1}. Distance: {2:.3f} A. FE-O-X angle: {3:.3f}".format(
                    n['source', 'serial'], n['source', 'name'], dist, fe_ox_angle))
        else:
            logger.debug(
                "Ligand target atom {0}-{1} excluded. Angle ({2:.3f}) or distance ({3:.3f}) criteria violated".format(
                    n['source', 'serial'], n['source', 'name'], fe_ox_angle, dist))

    return contact_frame


def eval_pication(contact_frame, topology, min_dist=0.05, pication_dist_max=0.6, pication_offset_max=0.2,
                  pication_amine_angle_dev=30.0, cation_attypes=('N.3', 'N.4', 'C.cat')):
    """
    Evaluate pi-cation interaction between aromatic rings and positively
    charged groups.

    Detect cations and aromatic rings in both source and target selections.
    Compute te pi-cation interaction using the `is_pication` function.

    :param contact_frame:      contacts DataFrame
    :type contact_frame:       :pandas:DataFrame
    :param topology:           main system topology
    :type topology:            :interact:TopologyDataFrame
    :param pication_dist_max:  maximum distance for pi-cation interactions (nm)
    :type pication_dist_max:   :py:float
    :param pication_offset_max Cutoff distance between geometric centers
                               projected on top of each other.
    :type pication_offset_max: :py:float
    :param pication_amine_angle_dev:  Maximum angle deviation between amine and
                                      ring normals.
    :type pication_amine_angle_dev:   :py:float
    :param min_dist:           minimum interaction distance (nm)
    :type min_dist:            :py:float
    :param cation_attypes:     SYBYL cation atom types
    :type cation_attypes:      :py:tuple

    :return:                   contacts DataFrame
    :rtype:                    :pandas:DataFrame
    """

    # Preselect all distances within pication_dist_max that involve cation and
    # aromatic ring atom types
    attypes = list(cation_attypes) + ['C.ar', 'N.ar']
    pcsel = contact_frame[(contact_frame['target', 'distance'] < pication_dist_max) &
                          (contact_frame['target', 'attype'].isin(attypes)) &
                          (contact_frame['source', 'attype'].isin(attypes))]

    if pcsel.empty:
        return contact_frame

    # Get source and/or target aromatic rings
    rings = {'source': [], 'target': []}
    for source_residue in topology[topology['serial'].isin(pcsel[('source', 'serial')])].residues(extend=True):
        rings['source'].extend(source_residue.find_rings(aromatic=True))

    for target_residue in topology[topology['serial'].isin(pcsel[('target', 'serial')])].residues(extend=True):
        rings['target'].extend(target_residue.find_rings(aromatic=True))

    if not rings['source'] and not rings['target']:
        return contact_frame

    logger.info('Evaluate {0} potential pi-cation interactions: pication_dist_max={1:.2f}, pication_offset_max={2:.1f},'
                ' pication_amine_angle_min={3:.1f}'.format(
        len(pcsel), pication_dist_max, pication_offset_max, pication_amine_angle_dev))

    # Loop over source and target aromatic ring sets and evaluate pi-cation interactions
    for group in (('source', 'target'), ('target', 'source')):

        # Get possible cations or continue
        atom_index = set(pcsel[group[0], 'serial'])
        cation = topology[(topology['serial'].isin(atom_index)) &
                          (topology['attype'].isin(cation_attypes))]
        if cation.empty:
            logger.debug('No cations found in {0}'.format(group[0]))
            continue

        # Do we have aromatic rings as target? else continue
        if not rings[group[1]]:
            logger.debug('No aromatic rings found in {0}'.format(group[1]))
            continue

        # Determine pi-cation
        for idc, cat in cation.iterrows():

            for ring in rings[group[1]]:

                iscation, data = is_pication(cat, ring,
                                             min_dist=min_dist,
                                             pication_dist_max=pication_dist_max,
                                             pication_offset_max=pication_offset_max,
                                             pication_amine_angle_dev=pication_amine_angle_dev)

                if iscation:
                    newindex = max(contact_frame.index) + 1
                    contact_frame.loc[newindex, 'contact'] = 'pc'
                    contact_frame.loc[newindex, ('target', 'distance')] = data['distance']
                    for label in ['segmentID', 'chainID', 'resName', 'resSeq']:
                        contact_frame.loc[newindex, (group[1], label)] = ring[label].unique()[0]
                    contact_frame.loc[newindex, ('target', 'serial')] = -1
                    contact_frame.loc[newindex, ('target', 'name')] = 'X1'
                    contact_frame.loc[newindex, ('target', 'attype')] = 'Du'
                    contact_frame.loc[newindex, ('target', 'element')] = 'D'

                    for label in ['serial', 'name', 'element', 'resSeq', 'resName', 'chainID', 'segmentID',
                                  'attype', 'charge']:
                        contact_frame.loc[newindex, (group[0], label)] = cat[label]

    return contact_frame
