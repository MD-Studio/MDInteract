# -*- coding: utf-8 -*-

"""
file: geometry.py

Basic functions for calculating geometrical parameters in a atomic system.
All functions operate on a single or set of atomic coordinates defined as
`numpy.ndarray`.
"""

import numpy
import math
import itertools

from pandas import Series

from interact import reference_data


def angle(coor1, coor2, coor3, deg=True):
    """
    Compute angle between three coordinates.

    :param coor1: first coordinate
    :type coor1:  :numpy:ndarray
    :param coor2: second coordinate
    :type coor2:  :numpy:ndarray
    :param coor3: third coordinate
    :type coor3:  :numpy:ndarray
    :param deg:   return angle in degrees or radians
    :type deg:    :py:bool

    :return:      angle
    :rtype:       :py:float
    """

    # Calculate vectors
    v1 = coor1 - coor2
    v2 = coor3 - coor2

    return vector_angle(v1, v2, deg=deg)


def dihedral(coor1, coor2, coor3, coor4, deg=True):
    """
    Compute dihedral angle between four coordinates.

    Compute dihedral angle between four coordinates defining the bond for
    which the torsion is calculated (~) as:

        V1 - V2 ~ V3 - V4

    :param coor1: first coordinate
    :type coor1:  :numpy:ndarray
    :param coor2: second coordinate
    :type coor2:  :numpy:ndarray
    :param coor3: third coordinate
    :type coor3:  :numpy:ndarray
    :param coor4: fourth coordinate
    :type coor4:  :numpy:ndarray
    :param deg:   return angle in degrees or radians
    :type deg:    :py:bool

    :return:      dihedral angle
    :rtype:       :py:float
    :raises:      ArithmeticError if the dihedral angle cant be calculated
                  (because vectors are collinear)
    """

    b0 = -1.0 * (coor2 - coor1)
    b1 = coor3 - coor2
    b2 = coor4 - coor3

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= numpy.linalg.norm(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - numpy.dot(b0, b1) * b1
    w = b2 - numpy.dot(b2, b1) * b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = numpy.dot(v, w)
    y = numpy.dot(numpy.cross(b1, v), w)
    rad = numpy.arctan2(y, x)

    if deg:
        return numpy.degrees(rad)
    return rad


def distance(coor1, coor2):
    """
    Calculate Euclidean distance between two coordinates.

    :param coor1: first coordinate
    :type coor1:  :numpy:ndarray
    :param coor2: second coordinate
    :type coor2:  :numpy:ndarray

    :return:      distance
    :rtype:       :py:float
    """

    return numpy.linalg.norm(coor1 - coor2)


def plane_fit(coordinates, center=None):
    """
    Compute n-dimensional plane to a series or coordinates and return the
    normal to that plane

    :param coordinates: coordinates to fit a plane to
    :type coordinates:  :numpy:ndarray
    :param center:      point on the plane used as origin for the normal
                        geometric center by default.
    :type center:       :numpy:ndarray

    :return:            normal to the fitted plane
    :rtype:             :numpy:ndarray
    """

    if center is None:
        center = coordinates.mean(axis=0)

    x = coordinates.T - center[:, None]

    return numpy.linalg.svd(x)[0][:, -1]


def projection(pnormal1, ppoint, tpoint):
    """
    Compute the centroid from a 3D point cloud and returns the coordinates

    :param pnormal1: normal of plane
    :type pnormal1:  :numpy:ndarray
    :param ppoint:   coordinates of point in the plane
    :type ppoint:    :numpy:ndarray
    :param tpoint:   coordinates of point to be projected
    :type tpoint:    :numpy:ndarray

    :returns:        coordinates of point orthogonally projected on the plane
    :rtype:          :numpy:ndarray
    """

    # Choose the plane normal pointing to the point to be projected
    pnormal2 = pnormal1 * -1
    d1 = distance(tpoint, pnormal1 + ppoint)
    d2 = distance(tpoint, pnormal2 + ppoint)
    pnormal = pnormal1 if d1 < d2 else pnormal2

    # Calculate the projection of tpoint to the plane
    sn = -numpy.dot(pnormal, tpoint - ppoint)
    sd = numpy.dot(pnormal, pnormal)
    sb = sn / sd

    return numpy.array([c1 + c2 for c1, c2 in zip(tpoint, [sb * pn for pn in pnormal])])


def vector_angle(v1, v2, deg=True):
    """
    Compute angle between two vectors

    :param v1:    first vector
    :type v1:     :numpy:ndarray
    :param v2:    second vector
    :type v2:     :numpy:ndarray
    :param deg:   return angle in degrees or radians
    :type deg:    :py:bool

    :return:      angle
    :rtype:       :py:float
    :raises:      ArithmeticError, if two vectors are identical
    """

    # Calculate angle in radians
    cosine_angle = numpy.dot(v1, v2) / (numpy.linalg.norm(v1) * numpy.linalg.norm(v2))
    if numpy.round(cosine_angle, 5) == 1.0:
        raise ArithmeticError('Unable to calculate angle between two identical vectors')

    rad = numpy.arccos(cosine_angle)

    if deg:
        return numpy.degrees(rad)
    return rad


def is_planar(coordinates, max_div=7.5):
    """
    Asses planarity of a ring structure based on the dihedral angles spanning
    the ring. None of the dihedrals should deviate more that `max_div` from
    an ideal 0/180 degree angle

    :param coordinates: ring coordinates
    :type coordinates:  :numpy:ndarray
    :param max_div:     maximum angle deviation in degrees
    :type max_div:      :py:float

    :return:            is planar or not
    :rtype:             :py:bool
    """

    for n in itertools.combinations(coordinates, 4):
        dangle = abs(dihedral(n[0], n[1], n[2], n[3]))
        if 0 <= dangle <= max_div or 180 - max_div <= dangle <= 180:
            continue
        else:
            return False

    return True


def radius_gyration(topology, mass=True):
    """
    Compute the radius of gyration of the atom selection

    :param topology: atom selection to compute radius of gyration for
    :type topology:  :interact:TopologyDataFrame
    :param mass:     use atomic masses otherwise equal masses of 1
    :type mass:      :py:bool

    :return:         radius of gyration
    :rtype:          :py:float
    """

    coords = topology.coord

    if mass:
        elements = reference_data['element_data']
        atom_mass = Series(elements.atomicMass.values, index=elements.symbol).to_dict()
        scale = numpy.array([atom_mass.get(element, 12.0) for element in topology['element']])
    else:
        scale = numpy.ones((len(coords), 1))

    weights = scale / scale.sum()

    mu = coords.mean(1)
    centered = (coords.transpose((1, 0, 2)) - mu).transpose((1, 0, 2))
    squared_dists = (centered ** 2).sum(2)

    return (squared_dists * weights).sum(1) ** 0.5


def rotation_matrix(rot_angle, direction):
    """
    Define rotation matrix over a given angle and direction vector.

    TODO: Beter make this one use quaternion
    """

    sina = math.sin(rot_angle)
    cosa = math.cos(rot_angle)
    direction = unit_vector(direction[:3])

    # rotation matrix around unit vector
    r = numpy.diag([cosa, cosa, cosa])
    r += numpy.outer(direction, direction) * (1.0 - cosa)
    direction *= sina
    r += numpy.array([[0.0, -direction[2], direction[1]],
                      [direction[2], 0.0, -direction[0]],
                      [-direction[1], direction[0], 0.0]])

    return r


def unit_vector(data, axis=None, out=None):
    """
    Return the unit vector
    """

    if out is None:
        data = numpy.array(data, dtype=numpy.float64, copy=True)
        if data.ndim == 1:
            data /= math.sqrt(numpy.dot(data, data))
            return data
    else:
        if out is not data:
            out[:] = numpy.array(data, copy=False)
        data = out
    length = numpy.atleast_1d(numpy.sum(data * data, axis))
    numpy.sqrt(length, length)
    if axis is not None:
        length = numpy.expand_dims(length, axis)
    data /= length
    if out is None:
        return data
