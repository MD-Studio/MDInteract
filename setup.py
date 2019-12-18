# -*- coding: utf-8 -*-

# package: interact
# file: setup.py
#
# Python library for biomolecular interaction profiling with support for
# Molecular Dynamics trajectories and WAMP API for use as microservice
# in the MDStudio environment.
#
# Copyright © 2016 Marc van Dijk, VU University Amsterdam, the Netherlands
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from setuptools import setup, find_packages

distribution_name = 'interact'

setup(
    name='mdinteract',
    version=0.2,
    description='A Pandas interface to the analysis of biomolecular structures with MD trajectory support',
    author='Marc van Dijk - VU University - Amsterdam',
    author_email='m4.van.dijk@vu.nl',
    url='https://github.com/MD-Studio/MDInteract',
    download_url='https://github.com/MD-Studio/MDInteract',
    license='Apache Software License 2.0',
    keywords='MDStudio biomolecules interactions trajectories pandas',
    platforms=['Any'],
    packages=find_packages(),
    py_modules=[distribution_name],
    test_suite="tests",
    install_requires=['mdtraj', 'pandas', 'scipy', 'psutil'],
    tests_require=['numpy'],
    package_data={distribution_name: ['data/*.csv']},
    include_package_data=True,
    zip_safe=True,
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: Apache Software License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 3',
        'Topic :: Software Development :: Libraries',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry'
    ]
)
