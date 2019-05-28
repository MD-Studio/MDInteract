# -*- coding: utf-8 -*-

import os
import re
import logging

from pandas import DataFrame


def mol2_to_dataframe(mol2_file, parse_multi_model=False, parse_coord=False,
                      columns=('serial', 'name', 'x', 'y', 'z', 'resSeq', 'resName', 'attype', 'charge', 'model')):
    """
    Parse a Tripos MOL2 file format to a Pandas DataFrame

    Uses the same column headers as the main TopologyDataFrame.
    Both atom coordinates and multiple models of a multi model MOL2 file are
    not parsed by default. To include them enable the `parse_coord` and
    'parse_multi_model` arguments respectively.

    MOL2 is a free format (no fixed column width). Their should be 9 columns
    and at least one empty space between each subsequent value on a line.
    The parser will raise an exception if this is not the case.

    :param mol2_file:           path to mol2 file
    :type mol2_file:            :py:str
    :param parse_multi_model:   parse all models from a multi-model MOL2 file
    :type parse_multi_model:    :py:bool
    :param parse_coord:         parse structure coordinates
    :type parse_coord:          :py:bool
    :param columns:             default DataFrame column headers
    :type columns:              :py:str

    :return:                    MOL2 as Pandas DataFrame
    :rtype:                     :pandas:DataFrame
    :raises:                    IOError on failed line parsing
    """

    # Check for file
    if not os.path.isfile(mol2_file):
        raise IOError('Tripos MOL2 file does not exist: {0}'.format(mol2_file))

    # Prepare DataFrame columns
    mol2_dict = dict([(n, []) for n in columns])

    if not parse_multi_model:
        del mol2_dict['model']
    if not parse_coord:
        for col in ('x', 'y', 'z'):
            del mol2_dict[col]

    read = False
    model = 0
    with open(mol2_file) as mf:
        for line in mf.readlines():

            # Start parsing after '@<TRIPOS>ATOM'
            if line.startswith('@<TRIPOS>ATOM'):

                # If more then one model and parse_multi_model False, stop
                if model == 1 and not parse_multi_model:
                    logging.debug('MOL2 likely contains multiple models. Only parsing the first')
                    break

                read = True
                model += 1
                continue

            # Stop parsing at '@<TRIPOS>BOND'
            elif line.startswith('@<TRIPOS>BOND'):
                read = False
                break

            # Parse MOL2 atom lines
            if read:

                # Proper MOL2 file should have 9 columns
                split_line = line.split()
                if len(split_line) < 9:
                    raise IOError('FormatError in mol2. Line: {0}'.format(line))

                try:
                    mol2_dict['serial'].append(int(split_line[0]))
                    mol2_dict['name'].append(split_line[1].upper())

                    # Parse atom coordinates or not
                    if parse_coord:
                        mol2_dict['x'].append(float(split_line[2]))
                        mol2_dict['y'].append(float(split_line[3]))
                        mol2_dict['z'].append(float(split_line[4]))

                    mol2_dict['attype'].append(split_line[5])
                    mol2_dict['resSeq'].append(int(split_line[6]))
                    mol2_dict['resName'].append(re.sub('{0}$'.format(split_line[6]), '', split_line[7]))
                    mol2_dict['charge'].append(float(split_line[8]))

                except ValueError as e:
                    raise IOError('FormatError in mol2. Line: {0}, error {1}'.format(line, e))

                # Add model identifier
                if parse_multi_model:
                    mol2_dict['model'].append(model)

    # Prepare the Pandas DataFrame
    return DataFrame(mol2_dict)
