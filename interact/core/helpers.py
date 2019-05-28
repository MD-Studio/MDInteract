# -*- coding: utf-8 -*-

from interact.core.geometry import is_planar


def renumber(topology, start=None, column='resSeq', retain_gap=True):
    """
    Renumber a numeric column in the TopologyDataFrame

    Select the column to renumber by column header name using the `column`
    argument. Renumber starting from a new `start` number or the first number
    in the column.
    If there is a gab in the numeric sequence the `retain_gab` argument will
    correct for this in the renumbered sequence.

    :param topology:
    :type topology:     :interact:TopologyDataFrame
    :param start:       start number
    :type:              :py:int
    :param column:      column name of numeric column to renumber
    :type column:       :py:str
    :param retain_gap:  preserve gaps in the numbering
    :type retain_gap:   :py:bool

    :return:
    """

    if column not in topology.columns:
        raise TypeError('No such column: {0}'.format(column))

    numbers = topology[column].values
    if start is None:
        start = int(numbers[0])

    renumbered = []
    current = None
    for nr in numbers:
        if current is None:
            current = nr

        if current != nr:

            if retain_gap:
                diff = nr - current
                if diff > 1:
                    start += diff
                else:
                    start += 1
            else:
                start += 1

            current = nr

        renumbered.append(start)

    if len(renumbered) != len(numbers):
        raise ArithmeticError('Renumbering failed, length mismatch')

    topology[column] = renumbered
    return topology


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


def is_aromatic(ring, max_div=7.5, aromatic_attypes=('C.ar', 'N.ar', 'C.3', 'C.2', 'N.3', 'N.2', 'N.pl3', 'O.2')):
    """
    Evaluate ring aromaticity

    Assumes that provided structure is a closed ring

    :param ring:             ring atom selection
    :type ring:              :interact:TopologyDataFrame
    :param max_div:          maximum ring planarity deviation
    :type max_div:           :py:float
    :param aromatic_attypes: sybyl atom types in aromatic ring
    :type aromatic_attypes:  :py:tuple

    :return:                 is aromatic or not
    :rtype:                  :py:bool
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
