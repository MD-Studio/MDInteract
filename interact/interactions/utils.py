# -*- coding: utf-8 -*-

import logging

from interact import __module__
from interact.core.geometry import distance, plane_fit, projection, vector_angle

logger = logging.getLogger(__module__)


def is_pication(cation, ring, min_dist=0.05, pication_dist_max=0.6, pication_offset_max=0.2,
                pication_amine_angle_dev=30.0):
    """
    Compute if a pi-cation interaction could exist between a cation and an
    aromatic ring

    Rules for pi-Cation interactions:
    1) Distance between cation and aromatic ring center should be between
       `min_dist` and `pication_dist_max`
    2) The distance offset between the ring center and the cation after
       projection onto the ring plane should not be more then
       `pication_offset_max`.
    3) If the cation is an amine that is linked with the remainder of the
       residue through more than one covalent bond (less degrees of freedom),
       the ring center should be positioned above the amine valence electrons
       computed as the angle between the normal of the ring and the normal
       defined by the heavy atom neighbours of the amine nitrogen that should
       be no more then `pication_amine_angle_dev`

    The function returns a boolean indicating if there is a pi-cation
    interaction and a dictionary with the computed metrics.
    notably: ring_center, ring_normal, distance, offset and cat_angle and
    cat_normal if it concerns an amine.

    :param cation:                    Cation atom selection
    :type cation:                     :interact:TopologyDataFrame,
                                      :interact:TopologySeries
    :param ring:                      Atom selection of the aromatic ring
    :type ring:                       :interact:TopologyDataFrame
    :param min_dist:                  minimum interaction distance (nm)
    :type min_dist:                   :py:float
    :param pication_dist_max:         maximum distance for pi-cation
                                      interactions (nm)
    :type pication_dist_max:          :py:float
    :param pication_offset_max:       Cutoff distance between geometric centers
    :type pication_offset_max:        :py:float
    :param pication_amine_angle_dev:  Maximum angle deviation between amine and
                                      ring normals.
    :type pication_amine_angle_dev:   :py:float

    :return:                          pi-cation interaction or not + data
    :rtype:                           :py:bool, py:dict
    """

    ispicat = False
    data = {}

    # Calculate distance between cation and ring center.
    # Calculate offset between ring center and cation projected onto ring plane
    ring_center = ring.center()
    ring_normal = plane_fit(ring.coord, center=ring_center)

    pcdist = distance(cation.coord, ring_center)
    pcoffset = distance(projection(ring_normal, ring_center, cation.coord), ring_center)

    data.update({'ring_center': ring_center, 'ring_normal': ring_normal, 'distance': pcdist, 'offset': pcoffset})
    if min_dist < pcdist < pication_dist_max and pcoffset < pication_offset_max:
        ispicat = True

        # If it concerns an tertiary or quarternary amine. Check angles.
        # Otherwise, we might have have a pi-cation interaction 'through' the ligand
        if cation.attype in ('N.3', 'N.4'):
            neigh = cation.neighbours(covalent=True)
            nonh = neigh[neigh['attype'] != 'H']

            # Count number of heavy atom neighbours that themselves are linked
            links = 0
            for i, n in nonh.iterrows():
                links += (len(n.neighbours(covalent=True)) -1 )

            if len(nonh) > 2 and links >= 2:

                # Calculate normal to plane defined by covalent neighbours of cation
                # Calculate angle between ring and cation normal
                cation_normal = plane_fit(nonh.coord, center=cation.coord)
                cation_angle = vector_angle(ring_normal, cation_normal)
                cation_angle = min(cation_angle, 180 - cation_angle if not 180 - cation_angle < 0 else cation_angle)

                data['cat_angle'] = cation_angle
                data['cat_normal'] = cation_normal

                # Vector angle should not deviate more then pication_amine_angle_dev
                if not cation_angle < pication_amine_angle_dev:
                    ispicat = False

                logging.debug('Cation likely an amine. Angle to ring normal: {0:.2f} deg.'.format(cation_angle))

    if ispicat:
        logger.info('Cation-pi interaction between {0}-{1} and ring {2}-{3}. Distance: {4:.3f} nm'
                    ' Offset:{5:.2f} nm'.format(cation.resName, cation.name, ring['resName'].values[0],
                                               ring['resSeq'].values[0], pcdist, pcoffset))

    return ispicat, data


def is_pistack(sring, tring, pistack_dist_max=0.55, pistack_ang_dev=30.0, min_dist=0.05, pistack_offset_max=0.20):

    pistack = False
    pistack_data = {'type': None}

    # Calculate source ring geometric center and normal to it
    sring_center = sring.center()
    sring_normal = plane_fit(sring.coord, center=sring_center)

    # Calculate target ring geometric center and normal to it
    tring_center = tring.center()
    tring_normal = plane_fit(tring.coord, center=tring_center) * -1

    # Calculate distance between ring centers
    dist = distance(sring_center, tring_center)

    # Rule 1: distance between ring centers
    if min_dist < dist < pistack_dist_max:

        # Calculate ring offset, (project each ring center into the other ring)
        proj1 = projection(sring_normal, sring_center, tring_center)
        proj2 = projection(tring_normal, tring_center, sring_center)
        offset = min(distance(proj1, sring_center), distance(proj2, tring_center))

        # Calculate angles between normals
        # Select smallest  of the two depending on direction
        a = vector_angle(sring_normal, tring_normal, deg=True)
        a = min(a, 180 - a if not 180 - a < 0 else a)

        pistack_data.update({'distance': dist, 'angle': a, 'offset': offset, 'sring_center': sring_center,
                             'tring_center': tring_center, 'sring_normal': sring_normal,
                             'tring_normal': tring_normal, 'sring': tuple(sring['serial']),
                             'tring': tuple(tring['serial'])})

        # Rule 2: pi-stacking
        if 0 < a < pistack_ang_dev and offset < pistack_offset_max:
            pistack_data['type'] = 'ps'
            pistack = True

        # Rule 3: T-stacking
        elif 90 - pistack_ang_dev < a < 90 + pistack_ang_dev and offset < pistack_offset_max:
            pistack_data['type'] = 'ts'
            pistack = True

        if pistack:
            logger.info('Pi-stacking type: {0} between {1}-{2} {3} and {4}-{5} {6}'.format(pistack_data['type'],
                                                                                         set(sring['resSeq']),
                                                                                         set(sring['resName']),
                                                                                         tuple(sring['serial']),
                                                                                         set(tring['resSeq']),
                                                                                         set(tring['resName']),
                                                                                         tuple(tring['serial'])))
            logger.info('Distance: {distance:.3f} nm, angle: {angle:.2f} deg, offset: {offset:.3f} nm'.format(
                **pistack_data))

    return pistack, pistack_data


def filter_contacts(cf, source_sel, target_sel, filter_sel=('target', 'distance')):
    """
    Filter contact dataframe for multiple occurrences of source-target
    selections based on a filter criterion.

    :param cf:          Contact DataFrame
    :type cf:           :pandas:DataFrame
    :param source_sel:  Column label for source selection
    :type source_sel:   :py:str
    :param target_sel:  Column label for target selection
    :type target_sel:   :py:sel
    :param filter_sel:  Column label to sort selection on
                        and filter all but the first row

    :return:            Filtered contact DataFrame
    :rtype:             :pandas:DataFrame
    """

    exclude_index = []
    for pair in (('source', 'target'), ('target', 'source')):

        # Get source selection
        source_mult = cf[(pair[0], source_sel)].value_counts()
        source_mult = set(source_mult[source_mult > 1].index)

        # Get target selection
        target_mult = cf[(pair[1], target_sel)].value_counts()
        target_mult = set(target_mult[target_mult > 1].index)

        # Filter
        for a in source_mult:
            for b in target_mult:
                sel = cf[(cf[(pair[0], source_sel)] == a) & (cf[(pair[1], target_sel)] == b)]
                if len(sel) > 1:
                    exclude_index.extend(list(sel.sort_values(by=filter_sel).iloc[1:].index))

    return cf.loc[~cf.index.isin(exclude_index)]
