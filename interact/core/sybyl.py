# -*- coding: utf-8 -*-

"""
file: sybyl.py

Functions for working with Tripos SYBYL atom types.
More info at:

    http://www.sdsc.edu/CCMS/Packages/cambridge/pluto/atom_types.html
"""

from numpy import nan

from interact.core.data import AA_SYBYL_TYPES


def parse_tripos_sybyl(mol2file):
    """
    Parse SYBYL atom types from a Tripos mol2 format file.

    :param mol2file:    Tripos mol2 file
    :type mol2file:     :py:str

    :return:            Dictionary with '<residue ID>-<atom name>' as
                        key and the SYBYL atom type as value.
    :rtype:             :py:dict
    """

    sybyldict = {}
    with open(mol2file) as molfile:
        read = False
        for line in molfile.readlines():
            if line.startswith('@<TRIPOS>ATOM'):
                read = True
                continue
            if line.startswith('@<TRIPOS>BOND'):
                break
            if read:
                line = line.split()
                atom = line[1]
                sybyl = line[5]
                resid = line[6]
                sybyldict['{0}-{1}'.format(resid, atom)] = sybyl

    return sybyldict


def assign_standard_sybyl_types(dataframe, residue_name='resName', atom_name='name'):
    """
    Assign default Tripos SYBYL atom types for amino acids and nucleic acids.

    Assigns based on residue and atom name identified in the DataFrame by the
    `residue_name` and `atom_name` arguments.
    The function adds a new 'attype' column to the dataframe with resolved
    SYBYL types or 'nan' if not.

    :param dataframe:       System DataFrame
    :type dataframe:        :pandas:DataFrame
    :param residue_name:    Residue name column header name in DataFrame
    :type residue_name:     :py:str
    :param atom_name:       Atom name column header name in DataFrame
    :type atom_name:        :py:str

    :return:                System DataFrame with `attype` column added
    :rtype:                 :pandas:DataFrame
    """

    resolved_sybyl_types = []
    for data in dataframe[[residue_name, atom_name]].values:
        resolved_sybyl_types.append(AA_SYBYL_TYPES.get('{0}-{1}'.format(*data), nan))

    dataframe['attype'] = resolved_sybyl_types

    return dataframe
