# -*- coding: utf-8 -*-

import os

__module__ = 'interact'
__docformat__ = 'restructuredtext'
__version__ = '{major:d}.{minor:d}'.format(major=0, minor=1)
__author__ = 'Marc van Dijk'
__status__ = 'pre-release beta1'
__date__ = '13 August 2018'
__licence__ = 'Apache Software License 2.0'
__url__ = 'https://github.com/MD-Studio/MDInteract'
__copyright__ = "Copyright (c) VU University, Amsterdam"
__rootpath__ = os.path.dirname(os.path.abspath(__file__))
__all__ = ['System', 'reference_data', 'constants', '__module__']

from glob import glob
from pandas import read_csv

# Import constants
from interact.constants import constants

# Import reference data sets
# The datasets are loaded once into Pandas DataFrame's and may be changed
# by the user.
reference_data = {}
for ref_data_file in glob(os.path.join(__rootpath__, 'data/*.csv')):
    dataset_name = os.path.basename(ref_data_file).split('.')[0]
    reference_data[dataset_name] = read_csv(ref_data_file, na_filter=False)

from interact.md_system import System