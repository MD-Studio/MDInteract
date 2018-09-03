# -*- coding: utf-8 -*-

from interact.core.geometry import is_planar

def set_contact_type(current, add):

    current = current.values[0].split()
    if 'nd' in current:
        current.remove('nd')
    add = add.split()

    return ' '.join(set(current + add))


def remove_contact_type(current, remove):

    current = current.values[0].split()
    if remove in current:
        current.remove(remove)

    if current:
        return ' '.join(set(current))
    return 'nd'

def is_aromatic(ring, max_div=7.5, aromatic_attypes={'C.ar', 'N.ar', 'C.3', 'C.2', 'N.3', 'N.2', 'N.pl3', 'O.2'}):
    """
    Evaluate ring aromaticity

    Assumes that provided structure is a closed ring

    :param ring:  ring atom selection
    :type ring:   :interact:TopologyDataFrmae

    :return:      is aromatic or not
    :rtype:       :py:bool
    """

    # Rule 1: ring should be planar
    if not is_planar(ring.coord, max_div=max_div):
        return False

    # Rule 2: fully conjugated
    if not set(ring['attype']).issubset(aromatic_attypes):
        return False

    # Rule 3: the molecule must have (4n+2) Pi electrons
    # TODO: how to check?
    return True