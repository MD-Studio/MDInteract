# -*- coding: utf-8 -*-

"""
file: md_system.py

Primary functions to parse Molecular Dynamics trajectory files or structure
files as molecular system and return TopologyDataFrame objects.
"""

import logging

from mdtraj import iterload
from pandas import concat
from psutil import virtual_memory

from interact import __module__
from interact.core.fileio import mol2_to_dataframe
from interact.core.topology_dataframe import TopologyDataFrame
from interact.core.sybyl import assign_standard_sybyl_types

logger = logging.getLogger(__module__)


class System(object):

    def __init__(self, topfile, mol2file=None, **kwargs):

        self.topfile = topfile
        self.mol2file = mol2file
        self.config = kwargs

        self.framecount = -1
        self.timestep = None
        self.starttime = 0
        self.nframes = None
        self.topology = None

        list(self.iter_frames(start=0, stop=1, step=1, chunk=1, auto_chunk=False))

    def __getitem__(self, frame):
        """
        Implement class __getitem__

        Load a single or selection of frames from the trajectory

        A single frame ID returns a single TopologyDataFrame object.
        A Python slice object in the form of [start:stop:step] will return a
        `iter_frames` generator.

        :param frame: single frame id or slice to load
        :type frame:  :py:int, :py:slice

        :return:      frame topology
        :rtype:       :interact:TopologyDataFrame
        :raises:      IndexError, frame index does not exist
        """

        if isinstance(frame, int):
            frames = list(self.iter_frames(start=frame, stop=frame+1))

            if not len(frames):
                raise IndexError('Frame {0} not in trajectory of {1} frames'.format(frame, self.framecount))

            return frames[0]

        return self.iter_frames(start=frame.start or 0, stop=frame.stop, step=frame.step or 1)

    def __iter__(self):
        """
        Implement class __iter__

        Returns the default `iter_frames` generator starting at frame 0 until
        the last frame of the trajectory with steps of 1 frame.

        :return:  trajectory frame iterator
        :rtype:   :interact:md_system:System:iter_frames
        """

        return self.iter_frames()

    def __len__(self):
        """
        Implement class __len__

        TODO: we now need to parse the full trajectory to count frames.
              Faster solutions?

        :return: number of frames in the trajectory.
        :rtype:  :py:int
        """

        if self.nframes is not None:
            return self.nframes

        frame = 0
        for top, frame in self.iter_frames(auto_chunk=True):
            pass

        return frame

    def _init_topology_dataframe(self, trajectory):
        """
        Initialize TopologyDataFrame from MDTraj trajectory

        A topology is initialized once based on a MDTraj trajectory chunk created
        by `mdtraj.iterload`

        :param trajectory:  MDTraj 'chunk'
        :type trajectory:   :mdtraj:core:trajectory:Trajectory

        :return:            TopologyDataFrame object
        :rtype:             :interact:TopologyDataFrame
        """

        # Get MDTraj trajectory info if possible
        try:
            self.timestep = trajectory.timestep
            self.starttime = trajectory.time[0]
        except (AttributeError, ValueError):
            self.nframes = 1

        # Create bonds if needed
        if not trajectory.topology.n_bonds:
            trajectory.topology.create_standard_bonds()

        # Export MDTraj topology to regular pandas DataFrame
        mdtraj_dataframe, standard_bonds = trajectory.topology.to_dataframe()

        # Add SYBYL atom type information to dataframe
        if self.mol2file:

            mol2_df = mol2_to_dataframe(self.mol2file)
            if len(mol2_df) == len(mdtraj_dataframe):
                mdtraj_dataframe = concat([mdtraj_dataframe, mol2_df[['attype', 'charge']]], axis=1)
            else:
                # First assign standard SYBYL types for AA and NA
                mdtraj_dataframe = assign_standard_sybyl_types(mdtraj_dataframe)

                # Then assign SYBYL types from the mol2 file
                df_slice = mdtraj_dataframe[(mdtraj_dataframe['resSeq'] == mol2_df['resSeq']) &
                                            (mdtraj_dataframe['name'] == mol2_df['name'])]
                df_slice['attype'] = mol2_df['attype']

        # Init TopologyDataFrame based on mdtraj_dataframe
        self.topology = TopologyDataFrame(mdtraj_dataframe.copy())

    def iter_frames(self, start=0, stop=None, step=1, chunk=100, auto_chunk=True, ave_byte_size=1000000):
        """
        Iteratively load frames from a trajectory

        Enables efficient out-of-core computation on frames by loading a set of
        `chunk` frames at a time into memory. The method uses `mdtraj.iterload`
        to do the heavy lifting of iteratively parsing the trajectory files.
        The chunk size can be set automatically using `auto_chunk`.
        This conservative approach takes half of the available free memory in
        bytes divided by an average frame size of 1000 bytes as chunk size.

        The frames to return can be controlled using the `start`, `stop` and
        `step` arguments. Frame count starts from frame 0. When `stop` is
        defined the iterator will read until the stop frame number if it
        exists but NOT include it. The `step` argument allows to skip a certain
        number of frames.

        :param start:           start reading at frame number
        :type start:            :py:int
        :param stop:            stop reading when reaching frame number
        :type stop:             :py:int
        :param step:            skip n number of frames
        :type step:             :py:int
        :param chunk:           number of frames read from disk into memory at
                                each load iterations.
        :type chunk:            :py:int
        :param auto_chunk:      automatically determine chunk size as function
                                of available memory
        :type auto_chunk:       :py:bool
        :param ave_byte_size:   average byte size of topology frame coordinate
                                set used for automatic definition of chunk size
        :type ave_byte_size:    :py:int

        :returns:               TopologyDataFrame and frame count
        :rtype:                 :interact:TopologyDataFrame, :py:int
        """

        # Determine chunk size as function of available memory
        if auto_chunk:
            free_memory = virtual_memory().free
            chunk = (free_memory / 2) / ave_byte_size
            chunk = int((chunk / step) * step)
            logger.info('set chunk size to: {0}'.format(chunk))

        # Iterload trajectory
        self.framecount = -1
        next_frame = start
        parse_chunks = True
        chunk_count = 0
        for traj_chunk in iterload(self.topfile, chunk=chunk, **self.config):

            chunk_count += chunk
            if not parse_chunks:
                break

            # continue if start frame not in current chunk range
            if next_frame > chunk_count + len(traj_chunk):
                self.framecount += chunk
                continue

            # Parse coordinates from every frame in chunk
            for counter, coords in enumerate(traj_chunk.xyz):

                self.framecount += 1
                if self.framecount != next_frame:
                    continue

                # 'stop' frame reached
                if stop:
                    if self.framecount >= stop:
                        parse_chunks = False
                        break

                # Set frame coordinates
                if self.topology is None:
                    self._init_topology_dataframe(traj_chunk)

                # Set frame coordinates in TopologyDataFrame
                self.topology.set_coord(coords)

                # Copy MDtraj unitcell information to TopologyDataFrame
                for arg in ('unitcell_vectors', 'unitcell_lengths', 'unitcell_angles', 'time'):
                    value = getattr(traj_chunk, arg, None)
                    if value is not None:
                        setattr(self.topology, arg, value[counter])

                # Define next frame to parse
                next_frame = self.framecount + step

                yield self.topology, self.framecount
