"""For reading, writing and manipulating MULgraph geometry grids.

Copyright 2011 University of Auckland.

This file is part of PyTOUGH.

PyTOUGH is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

PyTOUGH is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with PyTOUGH.  If not, see <http://www.gnu.org/licenses/>."""

from __future__ import print_function

from string import ascii_lowercase, ascii_uppercase
from geometry import *
from fixed_format_file import *

def padstring(s, length = 80): return s.ljust(length)

def int_to_chars(i, st = '', chars = ascii_lowercase,
                 spaces = True, length = 0):
    """Converts a number into a string of characters, using the specified
    characters. If spaces is False, no spaces are used in the name,
    and it is padded out with the first character in chars to the
    specified length."""
    def pad_to_length(st):
        return ''.join([chars[0] * (length - len(st)), st])
    if i > 0:
        n = len(chars)
        char_index = i - 1 if spaces else i
        st = int_to_chars(char_index // n,
                          ''.join([chars[char_index % n], st]),
                          chars, spaces, length = 0)
    return pad_to_length(st) if length and not spaces else st

def new_dict_key(d, istart = 0, justfn = str.rjust, length = 5,
                 chars = ascii_lowercase, spaces = True):
    """Returns an unused key for dictionary d, using the specified
    characters, plus the corresponding next starting index."""
    i = istart
    used = True
    while used:
        i += 1
        name = justfn(int_to_chars(i, chars = chars,
                                   spaces = spaces,
                                   length = length),
                      length)
        used = name in d
    return name, i

def uniqstring(s): return ''.join(sorted(set(s), key = s.index))

def fix_blockname(name):
    """Fixes blanks in 4th column of block names, caused by TOUGH2
    treating names as (a3, i2)"""
    if name[2].isdigit() and name[4].isdigit() and name[3] == ' ':
        return '0'.join((name[0:3], name[4:5]))
    else: return name

def unfix_blockname(name):
    """The inverse of fix_blockname()."""
    return "%3s%2d" % (name[0:3], int(name[3:5]))

def fix_block_mapping(blockmap):
    """Fixes block names in specified block mapping."""
    keys_to_fix = {}
    for k, v in blockmap.items():
        fixedk = fix_blockname(k)
        if k != fixedk: keys_to_fix[k] = fixedk
        blockmap[k] = fix_blockname(v)
    for k,v in keys_to_fix.items():
        item = blockmap[k]
        del blockmap[k]
        blockmap[v] = item

def valid_blockname(name):
    """Tests if a 5-character string is a valid blockname.  Allows names
    with the first three characters either letters, numbers, spaces or
    punctuation, the fourth character a digit or a space and the last
    character a digit.
    """
    from string import ascii_letters, digits, punctuation
    digit_space = digits + ' '
    letter_digit_space_punct = ascii_letters + digit_space + punctuation
    return all([s in letter_digit_space_punct for s in name[0:3]]) and \
        (name[3] in digit_space) and (name[4] in digits)

class NamingConventionError(Exception):
    """Used to raise exceptions when grid naming convention is not
    respected- e.g. when column or layer names are too long.
    """
    pass

class quadtree(object):
    """Quadtree for spatial searching in 2D grids."""

    def __init__(self, bounds, elements, parent = None):
        self.parent = parent
        self.bounds = bounds
        self.elements = elements
        self.child = []
        if self.parent:
            self.generation = self.parent.generation + 1
            self.all_elements = self.parent.all_elements
        else:
            self.generation = 0
            self.all_elements = set(elements)
        if self.num_elements > 1:
            rects = sub_rectangles(self.bounds)
            rect_elements = [[], [], [], []]
            for elt in self.elements:
                for irect, rect in enumerate(rects):
                    if in_rectangle(elt.centre, rect):
                        rect_elements[irect].append(elt)
                        break
            for rect, elts in zip(rects, rect_elements):
                if len(elts) > 0: self.child.append(quadtree(rect, elts, self))

    def __repr__(self): return self.bounds.__repr__()

    def get_num_elements(self): return len(self.elements)
    num_elements = property(get_num_elements)

    def get_num_children(self): return len(self.child)
    num_children = property(get_num_children)

    def search_wave(self, pos):
        from copy import copy
        todo = copy(self.elements)
        done = []
        while len(todo) > 0:
            elt = todo.pop(0)
            if elt.contains_point(pos): return elt
            done.append(elt)
            for nbr in elt.neighbour & self.all_elements:
                if rectangles_intersect(nbr.bounding_box, self.bounds) and \
                   not ((nbr in done) or (nbr in todo)):
                    todo.append(nbr)
        return None

    def search(self, pos):
        leaf = self.leaf(pos)
        if leaf: return leaf.search_wave(pos)
        else: return None

    def leaf(self, pos):
        if in_rectangle(pos, self.bounds):
            for child in self.child:
                childleaf = child.leaf(pos)
                if childleaf: return childleaf
            return self
        else: return None

    def plot(self, plt = None):
        if plt is None: import matplotlib.pyplot as plt
        x = [self.bounds[0][0], self.bounds[1][0], self.bounds[1][0],
             self.bounds[0][0], self.bounds[0][0]]
        y = [self.bounds[0][1], self.bounds[0][1], self.bounds[1][1],
             self.bounds[1][1], self.bounds[0][1]]
        plt.plot(x, y, '.--')
        for child in self.child: child.plot(plt)

mulgrid_format_specification = {
    'header': [['type', '_convention', '_atmosphere_type',
                'atmosphere_volume', 'atmosphere_connection', 
                'unit_type', 'gdcx', 'gdcy', 'cntype', 'permeability_angle'], 
               ['5s', '1d', '1d',
                '10.2e', '10.2e',
                '5s', '10.2f', '10.2f', '1d', '10.2f']], 
    'node': [['name', 'x', 'y'],  ['3s'] + ['10.2f'] * 2], 
    'column': [['name', 'centre_specified', 'num_nodes', 'xcentre', 'ycentre'],
               ['3s', '1d', '2d'] + ['10.2f'] * 2],
    'column_node': [['name'], ['3s']],
    'connection': [['name1', 'name2'], ['3s', '3s']],
    'layer': [['name', 'bottom', 'centre'], ['3s'] + ['10.2f'] * 2],
    'surface': [['name', 'elevation'], ['3s', '10.2f']],
    'well': [['name', 'x', 'y', 'z'], ['5s'] + ['10.1f'] * 3]}

class node(object):
    """Grid node class"""
    def __init__(self, name = '   ', pos = None):
        if pos is None: pos = np.array([0.0, 0.0])
        self.name = name
        if isinstance(pos, (tuple, list)): pos = np.array(pos)
        self.pos = pos
        self.column = set([])
    def __repr__(self): return self.name

class column(object):
    """Grid column class"""
    def __init__(self, name = '   ', node = None, centre = None, surface = None):
        if node is None: node = []
        self.name = name
        self.node = node
        if centre is None:
            self.centre_specified = 0
            if self.num_nodes > 0: self.centre = self.centroid
            else: self.centre = None
        else:
            self.centre_specified = 1
            self.centre = centre
        self.surface = surface
        self.get_area()
        if self.area < 0.: # check node numbering orientation
            self.node.reverse()
            self.area = -self.area
        self.neighbour = set([])
        self.connection = set([])
        self.num_layers = 0

    def get_num_nodes(self): return len(self.node)
    num_nodes = property(get_num_nodes)
    def get_num_neighbours(self): return len(self.neighbour)
    num_neighbours = property(get_num_neighbours)

    def get_surface(self): return self._surface
    def set_surface(self, val):
        self._surface = val
        if val is None: self.default_surface = True
        else: self.default_surface = False
    surface = property(get_surface, set_surface)

    def is_against(self, col):
        """Returns True if the column is against the specified other column-
        that is, if it shares more than one node with it."""
        return len(set(self.node).intersection(set(col.node))) > 1

    def get_polygon(self):
        """Returns polygon formed by node positions."""
        return [node.pos for node in self.node]
    polygon = property(get_polygon)

    def get_area(self):
        """Calculates column area"""
        self.area = polygon_area(self.polygon)

    def get_centroid(self):
        """Returns column centroid"""
        return polygon_centroid(self.polygon)
    centroid = property(get_centroid)

    def get_bounding_box(self):
        """Returns (horizontal) bounding box of the column."""
        return bounds_of_points([node.pos for node in self.node])
    bounding_box = property(get_bounding_box)

    def get_neighbourlist(self):
        """Returns a list of neighbouring columns corresponding to each column side (None if
        the column side is on a boundary)."""
        nbrlist = []
        for i, nodei in enumerate(self.node):
            i1 = (i + 1) % self.num_nodes
            nodes = set([nodei, self.node[i1]])
            con = [cn for cn in self.connection if set(cn.node) == nodes]
            if con: col = [c for c in con[0].column if c != self][0]
            else: col = None
            nbrlist.append(col)
        return nbrlist
    neighbourlist = property(get_neighbourlist)

    def near_point(self, pos):
        """Returns True if pos is within the bounding box of the column."""
        return in_rectangle(pos, self.bounding_box)

    def contains_point(self, pos):
        """Determines if specified point is inside the column."""
        return in_polygon(pos, self.polygon)

    def in_polygon(self, polygon):
        """Returns true if the centre of the column is inside the specified polygon."""
        if len(polygon) == 2: return in_rectangle(self.centre, polygon) # for rectangles
        else: return in_polygon(self.centre, polygon)

    def get_exterior_angles(self):
        """Returns list of exterior angle for each node in the column."""
        side = [self.node[i].pos - self.node[i - 1].pos for i in range(self.num_nodes)]
        h = [vector_heading(s) for s in side]
        angles = [np.pi - (h[(i + 1) % self.num_nodes] - h[i]) for i in range(self.num_nodes)]
        angles = [a % (2 * np.pi) for a in angles]
        return angles
    exterior_angles = property(get_exterior_angles)

    def get_interior_angles(self):
        return [2. * np.pi - a for a in self.exterior_angles]
    interior_angles = property(get_interior_angles)

    def get_angle_ratio(self):
        """Returns the angle ratio for the column, defined as the ratio of the
        largest interior angle to the smallest interior angle.
        """
        angles = self.interior_angles
        return max(angles) / min(angles)
    angle_ratio = property(get_angle_ratio)

    def get_side_lengths(self):
        "Returns list of side lengths for the column"
        return np.array([norm(self.node[(i + 1) % self.num_nodes].pos -
                              self.node[i].pos) for i in range(self.num_nodes)])
    side_lengths = property(get_side_lengths)

    def get_side_ratio(self):
        """Returns the side ratio for the column, defined as the ratio of the
        largest side length to the smallest side length (a
        generalisation of the aspect ratio for quadrilateral
        columns).
        """
        l = self.side_lengths
        return np.max(l) / np.min(l)
    side_ratio = property(get_side_ratio)

    def bisection_sides(self, direction = None):
        """Returns indices of column sides which should be used to bisect the
        column.  If direction is specified as 'x' or 'y', the column
        in bisected across the sides most closely aligned with that
        direction; otherwise, bisection is done for triangles across
        the two longest sides of the column, and for quadrilaterals
        across the longest side and its opposite.
        """
        if direction is None:
            l = self.side_lengths
            isort = np.argsort(l)
            if self.num_nodes == 3: return (isort[-1], isort[-2])
            elif self.num_nodes == 4:
                imax = isort[-1]
                iopp = (imax + 2) % self.num_nodes
                return (imax, iopp)
            else: return None
        else:
            n = np.array([[1., 0.], [0., 1.]][direction == 'y'])
            d, iside = [], []
            nn = self.num_nodes
            if nn in [3, 4]:
                for i in range(nn):
                    x1 = 0.5 * (self.node[i].pos + self.node[(i + 1) % nn].pos)
                    if self.num_nodes == 3: i2 = (i + 1) % nn
                    else: i2 = (i + 2) % nn
                    x2 = 0.5 * (self.node[i2].pos + self.node[(i2 + 1) % nn].pos)
                    d.append(abs(np.dot(x2 - x1, n)))
                    iside.append((i, i2))
                imax = np.argsort(d)
                return iside[imax[-1]]
            else: return None

    def basis(self, xi):
        """Returns bilinear 2D finite element basis functions for the column
        at the specified local coordinate."""
        if self.num_nodes == 3: return np.array([xi[0], xi[1], 1. - xi[0] - xi[1]])
        elif self.num_nodes == 4: # over [-1, 1]
            a0, a1, b0, b1 = 1. - xi[0], 1. + xi[0], 1. - xi[1], 1. + xi[1]
            return 0.25 * np.array([a0 * b0, a1 * b0, a1 * b1, a0 * b1])
        else: return None

    def basis_derivatives(self, xi):
        """Returns bilinear 2D finite element basis function derivatives for
        the column at the specified local coordinate."""
        if self.num_nodes == 3: return np.array([[1., 0.], [0., 1.], [-1., -1.]])
        elif self.num_nodes == 4:
            a0, a1, b0, b1 = 1. - xi[0], 1. + xi[0], 1. - xi[1], 1. + xi[1]
            return 0.25 * np.array([[-b0, -a0], [b0, -a1], [b1, a1], [-b1, a0]])
        else: return None

    def Jacobian(self, xi):
        """Returns bilinear 2D finite element Jacobian matrix for the column
        at the specified local coordinate."""
        dpsi = self.basis_derivatives(xi)
        J = np.zeros((2, 2))
        for i in range(2):
            for j in range(2):
                for k, nodek in enumerate(self.node): J[i, j] += dpsi[k, j] * nodek.pos[i]
        return J

    def global_pos(self, xi):
        """Returns global coordinates of the local point xi in the column."""
        psi = self.basis(xi)
        return sum([psi[i] * nodei.pos for i, nodei in enumerate(self.node)])

    def local_inside(self, xi):
        """Returns true if a local point is inside the column."""
        if self.num_nodes == 3: return all([x >= 0. for x in xi]) and (np.sum(xi) <= 1.)
        elif self.num_nodes == 4: return all([abs(x) <= 1. for x in xi])
        else: return None

    def local_pos(self, x):
        """Finds local coordinates of global point x in the column."""
        if self.num_nodes in [3, 4]:
            tolerance, max_iterations = 1.e-8, 15
            if self.num_nodes == 3: xi = np.array([1 / 3., 1 / 3.])
            else: xi = np.zeros(2)
            found = False
            for n in range(max_iterations): # Newton iteration
                dx = self.global_pos(xi) - x
                if np.linalg.norm(dx) <= tolerance:
                    found = True
                    break
                else:
                    J = self.Jacobian(xi)
                    try:
                        xi -= np.linalg.solve(J, dx)
                    except np.linalg.LinAlgError: break
            if not found: return None
            else:
                if self.local_inside(xi): return xi
                else: return None
        else: return None

    def index_plus(self, i, d):
        """Adds d to index i around column."""
        return (i + d) % self.num_nodes

    def index_minus(self, i, d):
        """Subtracts d from index i around column."""
        result = i - d
        if result < 0: result += self.num_nodes
        return result

    def index_dist(self, i1, i2):
        """Returns distance between two integer indices around the column."""
        d = abs(i1 - i2)
        if 2 * d > self.num_nodes: d = self.num_nodes - d
        return d

    def __repr__(self): return self.name

class connection(object):
    """Column connection class"""

    def __init__(self, col = None, nod = None):
        if col is None: col = [column(), column()]
        if nod is None: nod = [node(), node()]
        self.column = col
        self.node = nod

    def __repr__(self): return self.column[0].name + ':' + self.column[1].name

    def get_angle_cosine(self):
        """Returns cosine of angle between the connection face and the line
        joining the two columns in the connection.  Ideally want this
        to be zero, i.e. connecting line perpendicular to the face.
        """
        n = self.node[1].pos - self.node[0].pos
        n = n / norm(n)
        dcol = self.column[1].centre - self.column[0].centre
        d = dcol / norm(dcol)
        return np.dot(n, d)
    angle_cosine = property(get_angle_cosine)

class layer(object):
    """Grid layer class"""

    def __init__(self, name = '  ', bottom = 0.0, centre = 0.0, top = 0.0):
        self.name = name
        self.bottom = bottom
        self.centre = centre
        self.top = top

    def __repr__(self):
        return self.name + '(' + str(self.bottom) + ':' + str(self.top) + ')'

    def contains_elevation(self, z):
        return self.bottom <= z <= self.top

    def translate(self, shift):
        """Translates a layer up or down by specified distance"""
        self.top += shift
        self.bottom += shift
        self.centre += shift

    def get_thickness(self): return self.top - self.bottom
    thickness = property(get_thickness)

class well(object):
    """Well class"""

    def __init__(self, name = '  ', pos = None):
        if pos is None: pos = []
        self.name = name
        for i, p in enumerate(pos):
            if isinstance(p, (list, tuple)): pos[i] = np.array(p)
        self.pos = pos

    def __repr__(self): return self.name

    def get_num_pos(self): return len(self.pos)
    num_pos = property(get_num_pos)

    def get_num_deviations(self): return self.num_pos - 1
    num_deviations = property(get_num_deviations)

    def get_deviated(self): return self.num_deviations > 1
    deviated = property(get_deviated)

    def get_head(self): return self.pos[0]
    head = property(get_head)

    def get_bottom(self): return self.pos[-1]
    bottom = property(get_bottom)

    def pos_coordinate(self, index):
        """Returns array of specified coordinate in pos array."""
        return np.array([pos[index] for pos in self.pos])

    def get_pos_depth(self):
        """Returns array of downhole depths corresponding to pos array."""
        return np.cumsum([0.] + [np.linalg.norm(pos - self.pos[i]) for
                                 i, pos in enumerate(self.pos[1:])])
    pos_depth = property(get_pos_depth)

    def elevation_depth(self, elevation):
        """Returns downhole depth corresponding to a given elevation (or None
        if the specified elevation is outside the well).
        """
        epos = self.pos_coordinate(2)
        # NB: np.interp() needs abcissa to be increasing, so have to
        # reverse the arrays here:
        if epos[-1] <= elevation <= epos[0]:
            return np.interp(elevation, epos[::-1], self.pos_depth[::-1])
        else: return None

    def depth_elevation(self, depth):
        """Returns elevation corresponding to a given downhole depth (or None if the specified
        depth is outside the well)."""
        dpos = self.pos_depth
        if dpos[0] <= depth <= dpos[-1]:
            return np.interp(depth, dpos, self.pos_coordinate(2))
        else: return None

    def elevation_pos(self, elevation, extend = False):
        """Returns 3D position in well, given an elevation.  If extend is
        True, return extrapolated positions for elevations below the
        bottom of the well.
        """
        poscoord = [self.pos_coordinate(i) for i in range(3)]
        epos = poscoord[2]
        if epos[-1] <= elevation <= epos[0]:
            return np.array([np.interp(elevation, epos[::-1],
                                       poscoord[i][::-1]) for i in range(3)])
        elif elevation < epos[-1] and extend:
            # extrapolate last deviation:
            pbot = self.pos[-1]
            if self.num_pos > 1: ptop = self.pos[-2]
            else: ptop = np.array(list(pbot[0:2]) + [pbot[2] + 1.])
            ebot, etop = pbot[2], ptop[2]
            alpha = (elevation - ebot) / (etop - ebot)
            return (1. - alpha) * pbot + alpha * ptop
        else: return None

    def depth_pos(self, depth):
        """Returns 3D position in well, given a depth."""
        elevation = self.depth_elevation(depth)
        if elevation: return self.elevation_pos(elevation)
        else: return None

class mulgrid(object):
    """MULgraph grid class"""

    def __init__(self, filename = '', type = 'GENER', convention = 0,
                 atmos_type = 0, atmos_volume = 1.e25,
                 atmos_connection = 1.e-6, unit_type = '', permeability_angle = 0.0,
                 read_function = default_read_function):
        self.filename = filename
        self.type = type  # geometry type- only GENER supported
        self._convention = convention  # naming convention:
        # 0: 3-char column + 2-digit layer
        # 1: 3-char layer + 2-digit column
        # 2: 2-char layer + 3-digit column
        self._atmosphere_type = atmos_type  # atmosphere type:
        # 0: single atmosphere block
        # 1: one atmosphere block per column
        # else: no atmosphere blocks
        self.set_secondary_variables()
        self.atmosphere_volume = atmos_volume
        self.atmosphere_connection = atmos_connection
        self.unit_type = unit_type
        self.gdcx, self.gdcy = None, None
        self.cntype = None # not supported
        self.permeability_angle = permeability_angle
        self.read_function = read_function
        self.empty()
        if self.filename: self.read(filename)

    def set_secondary_variables(self):
        """Sets variables dependent on naming convention and atmosphere type"""
        if self.atmosphere_type == 0:
            self.atmosphere_column_name = ['ATM', ' 0', '  0'][self.convention]
        self.colname_length = [3, 2, 3][self.convention]
        self.layername_length = [2, 3, 2][self.convention]

    def get_convention(self):
        """Get naming convention"""
        return self._convention
    def set_convention(self, convention):
        """Set naming convention"""
        self._convention = convention
        self.set_secondary_variables()
        self.setup_block_name_index()
        self.setup_block_connection_name_index()
    convention = property(get_convention, set_convention)

    def get_atmosphere_type(self):
        """Get atmosphere type"""
        return self._atmosphere_type
    def set_atmosphere_type(self, atmos_type):
        """Set atmosphere type"""
        self._atmosphere_type = atmos_type
        self.set_secondary_variables()
        self.setup_block_name_index()
        self.setup_block_connection_name_index()
    atmosphere_type = property(get_atmosphere_type, set_atmosphere_type)

    def get_unit_type(self):
        """Get unit type"""
        return self._unit_type
    def set_unit_type(self, unit_type):
        """Set unit type"""
        self._unit_type = unit_type
        self.unit_scale = {'': 1.0, 'FEET ': 0.3048}[unit_type]
    unit_type = property(get_unit_type, set_unit_type)

    def get_area(self):
        """Grid area- sum of column areas"""
        return sum([col.area for col in self.columnlist])
    area = property(get_area)

    def get_centre(self):
        """Returns grid centre- approximated as area-weighted average of column centres"""
        if self.num_columns > 0:
            return sum([col.area * col.centre for col in self.columnlist]) / self.area
        else: return None
    centre = property(get_centre)

    def get_num_blocks(self):
        """Returns number of blocks in the tough2 grid represented by the geometry."""
        return len(self.block_name_list)
    num_blocks = property(get_num_blocks)

    def get_num_atmosphere_blocks(self):
        """Returns number of atmosphere blocks in the tough2 grid represented by the geometry."""
        return [1, self.num_columns, 0][self.atmosphere_type]
    num_atmosphere_blocks = property(get_num_atmosphere_blocks)

    def get_num_underground_blocks(self):
        """Returns number of blocks under the ground surface
        (i.e. non-atmosphere blocks) in the tough2 grid represented by
        the geometry.
        """
        return self.num_blocks - self.num_atmosphere_blocks
    num_underground_blocks = property(get_num_underground_blocks)

    def get_num_block_connections(self):
        """Returns number of connections between blocks in the TOUGH2 grid
        represented by the geometry."""
        return len(self.block_connection_name_list)
    num_block_connections  =  property(get_num_block_connections)

    def get_column_angle_ratio(self):
        """Returns an array of angle ratios for each column."""
        return np.array([col.angle_ratio for col in self.columnlist])
    column_angle_ratio = property(get_column_angle_ratio)

    def get_column_side_ratio(self):
        """Returns an array of side ratios for each column."""
        return np.array([col.side_ratio for col in self.columnlist])
    column_side_ratio = property(get_column_side_ratio)

    def get_connection_angle_cosine(self):
        """Returns an array of connection angle cosines, for each connection."""
        return np.array([con.angle_cosine for con in self.connectionlist])
    connection_angle_cosine = property(get_connection_angle_cosine)

    def get_tilt_vector(self):
        """Returns a tilt vector, used to calculate gravity cosines of TOUGH2
        grid connections when the GDCX or GDCY grid tilting options
        are used.
        """
        from math import sqrt
        gdcx = 0. if self.gdcx is None else self.gdcx
        gdcy = 0. if self.gdcy is None else self.gdcy
        def cosfromsin(sinangle): return sqrt(1. - min(sinangle * sinangle, 1.))
        sintheta = -gdcy
        costheta = cosfromsin(sintheta)
        try:
            sinphi = gdcx / costheta
            cosphi = cosfromsin(sinphi)
            return np.array([costheta * sinphi, -sintheta, -costheta * cosphi])
        except ZeroDivisionError: return np.array([0., -sintheta, 0.]) # theta = pi/2
    tilt_vector = property(get_tilt_vector)

    def empty(self):
        """Empties grid contents."""
        self.nodelist = []
        self.columnlist = []
        self.layerlist = []
        self.connectionlist = []
        self.welllist = []
        self.node = {}
        self.column = {}
        self.layer = {}
        self.connection = {}
        self.well = {}
        self.block_name_list = []
        self.block_name_index = {}
        self.block_connection_name_list = []
        self.block_connection_name_index = {}

    def __repr__(self):
        conventionstr = [
            '3 characters for column, 2 digits for layer',
            '3 characters for layer, 2 digits for column',
            '2 characters for layer, 3 digits for column'][self.convention]
        atmstr = [
            'single atmosphere block',
            'one atmosphere block over each column',
            'no atmosphere blocks'][self.atmosphere_type]
        return str(self.num_nodes) + ' nodes; ' + \
            str(self.num_columns) + ' columns; ' + \
            str(self.num_layers) + ' layers; ' + \
            str(self.num_blocks) + ' blocks; ' + \
            str(self.num_wells) + ' wells' + '\n' + \
            'Naming convention: ' + str(self.convention) + \
            ' (' + conventionstr + ')\n' + \
            'Atmosphere type  : ' + str(self.atmosphere_type) + \
            ' (' + atmstr + ')'

    def get_default_surface(self):
        return all([col.default_surface for col in self.columnlist])
    default_surface = property(get_default_surface)

    def get_num_nodes(self):
        return len(self.node)
    num_nodes = property(get_num_nodes)

    def get_num_columns(self):
        return len(self.column)
    num_columns = property(get_num_columns)

    def get_num_layers(self):
        return len(self.layer)
    num_layers = property(get_num_layers)

    def get_num_connections(self):
        return len(self.connectionlist)
    num_connections = property(get_num_connections)

    def get_num_wells(self):
        return len(self.well)
    num_wells = property(get_num_wells)

    def get_layer_index(self):
        return dict([(lay.name, i) for i, lay in enumerate(self.layerlist)])
    layer_index = property(get_layer_index)

    def get_column_index(self):
        return dict([(col.name, i) for i, col in enumerate(self.columnlist)])
    column_index = property(get_column_index)

    def connects(self, col1, col2):
        """Returns True if the geometry contains a connection connecting the
        two specified columns."""
        return any([(col1 in con.column) and (col2 in con.column) for
                    con in self.connectionlist])

    def setup_block_name_index(self):
        """Sets up list and dictionary of block names and indices for the
        tough2 grid represented by the geometry."""
        self.block_name_list = []
        if self.num_layers > 0:
            if self.atmosphere_type  ==  0: # one atmosphere block
                self.block_name_list.append(
                    self.block_name(self.layerlist[0].name, self.atmosphere_column_name))
            elif self.atmosphere_type == 1: # one atmosphere block per column
                for col in self.columnlist:
                    self.block_name_list.append(
                        self.block_name(self.layerlist[0].name, col.name))
            for lay in self.layerlist[1:]:
                for col in [col for col in self.columnlist if col.surface > lay.bottom]:
                    self.block_name_list.append(self.block_name(lay.name, col.name))
        self.block_name_index = dict([(blk, i) for i, blk in enumerate(self.block_name_list)])

    def setup_block_connection_name_index(self):
        """Sets up list and dictionary of connection names and indices for
        blocks in the TOUGH2 grid represented by the geometry."""
        self.block_connection_name_list = []
        for ilay, lay in enumerate(self.layerlist[1:]):
            layercols = [col for col in self.columnlist if col.surface > lay.bottom]
            for col in layercols: # vertical connections
                thisblkname = self.block_name(lay.name,  col.name)
                if (ilay == 0) or (col.surface <= lay.top): # connection to atmosphere
                    abovelayer = self.layerlist[0]
                    if self.atmosphere_type == 0:
                        aboveblkname = self.block_name_list[0]
                    elif self.atmosphere_type == 1:
                        aboveblkname = self.block_name(abovelayer.name, col.name)
                    else: continue
                else:
                    abovelayer = self.layerlist[ilay]
                    aboveblkname = self.block_name(abovelayer.name, col.name)
                self.block_connection_name_list.append((thisblkname,  aboveblkname))
            layercolset = set(layercols) # horizontal connections:
            cons = [con for con in self.connectionlist if
                    set(con.column).issubset(layercolset)]
            for con in cons:
                conblocknames = tuple([self.block_name(lay.name, concol.name) for
                                       concol in con.column])
                self.block_connection_name_list.append(conblocknames)
        self.block_connection_name_index = dict(
            [(con, i) for i, con in enumerate(self.block_connection_name_list)])

    def column_name(self, blockname):
        """Returns column name of block name."""
        if self.convention == 0: return blockname[0: 3]
        elif self.convention == 1: return blockname[3: 5]
        elif self.convention == 2: return blockname[2: 5]
        else: return None

    def layer_name(self, blockname):
        """Returns layer name of block name."""
        if self.convention == 0: return blockname[3: 5]
        elif self.convention == 1: return blockname[0: 3]
        elif self.convention == 2: return blockname[0: 2]
        else: return None

    def node_col_name_from_number(self, num, justfn = str.rjust,
                                  chars = ascii_lowercase, spaces = True):
        """Returns node or column name from number."""
        if self.convention == 0:
            name = justfn(int_to_chars(num, chars = chars, spaces = spaces,
                                       length = self.colname_length), self.colname_length)
        else: name = str.rjust(str(num), self.colname_length)
        return name

    def column_name_from_number(self, num, justfn = str.rjust,
                                chars = ascii_lowercase, spaces = True):
        """Returns column name from column number."""
        name = self.node_col_name_from_number(num, justfn, chars, spaces)
        if len(name) > self.colname_length:
            raise NamingConventionError(
                "Column name is too long for the grid naming convention.")
        return name

    def node_name_from_number(self, num, justfn = str.rjust,
                              chars = ascii_lowercase, spaces = True):
        """Returns node name from node number."""
        name = self.node_col_name_from_number(num, justfn, chars, spaces)
        if len(name) > self.colname_length:
            raise NamingConventionError(
                "Node name is too long for the grid naming convention.")
        return name

    def layer_name_from_number(self, num, justfn = str.rjust, chars = ascii_lowercase,
                               spaces = True):
        """Returns layer name from layer number."""
        if self.convention == 0:
            name = justfn(str(num), self.layername_length)
        else:
            name = justfn(int_to_chars(num, chars = chars, spaces = spaces,
                                       length = self.layername_length),
                          self.layername_length)
        if len(name) > self.layername_length:
            raise NamingConventionError(
                "Layer name is too long for the grid naming convention.")
        return name

    def get_uppercase_names(self):
        """Returns True if character part of block names are uppercase."""
        return all([(blkname[0:3] == blkname[0:3].upper()) for
                    blkname in self.block_name_list])
    uppercase_names = property(get_uppercase_names)

    def get_right_justified_names(self):
        """Returns True if character part of block names are right-justified."""
        return all([(blkname[0:3] == blkname[0:3].rjust(3)) for
                    blkname in self.block_name_list])
    right_justified_names = property(get_right_justified_names)

    def new_node_name(self, istart = 0, justfn = str.rjust, chars = ascii_lowercase,
                      spaces = True):
        name, i = new_dict_key(self.node, istart, justfn, self.colname_length,
                               chars, spaces)
        if len(name) > self.colname_length:
            raise NamingConventionError(
                "Node name is too long for the grid naming convention.")
        else: return name, i

    def new_column_name(self, istart = 0, justfn = str.rjust,
                        chars = ascii_lowercase, spaces = True):
        name, i = new_dict_key(self.column, istart, justfn, self.colname_length,
                               chars, spaces)
        if len(name) > self.colname_length:
            raise NamingConventionError(
                "Column name is too long for the grid naming convention.")
        else: return name, i

    def column_bounds(self, columns):
        """Returns horizontal bounding box for a list of columns."""
        nodes = self.nodes_in_columns(columns)
        return bounds_of_points([node.pos for node in nodes])

    def get_bounds(self):
        """Returns horizontal bounding box for grid."""
        return bounds_of_points([node.pos for node in self.nodelist])
    bounds = property(get_bounds)

    def add_node(self, nod = None):
        """Adds node to the geometry. If a node with the specified name
        already exists in the geometry, no new node is added."""
        if nod is None: nod = node()
        if nod.name not in self.node:
            self.nodelist.append(nod)
            self.node[nod.name] = self.nodelist[-1]

    def delete_node(self, nodename):
        """Deletes node from the geometry."""
        node = self.node[nodename]
        del self.node[nodename]
        self.nodelist.remove(node)

    def add_column(self, col = None):
        """Adds column to the geometry. If a column with the specified
        name already exists in the geometry, no new column is added."""
        if col is None: col = column()
        if col.name not in self.column:
            self.columnlist.append(col)
            self.column[col.name] = self.columnlist[-1]
            for node in col.node: node.column.add(col)

    def delete_column(self, colname):
        """Deletes a column from the geometry."""
        col = self.column[colname]
        cons = [con for con in self.connectionlist if col in con.column]
        for con in cons:
            self.delete_connection(tuple([c.name for c in con.column]))
        for nbr in col.neighbour: nbr.neighbour.remove(col)
        for node in col.node: node.column.remove(col)
        del self.column[colname]
        self.columnlist.remove(col)

    def split_column(self, colname, nodename, chars = ascii_lowercase):
        """Splits the specified quadrilateral column into two triangles,
        splitting at the specified node.  Returns True if the
        operation was successful.
        """
        chars = uniqstring(chars)
        justfn = [str.ljust, str.rjust][self.right_justified_names]
        if colname in self.column:
            col = self.column[colname]
            nn = col.num_nodes
            if nn == 4:
                nodenames = [node.name for node in col.node]
                try:
                    i0 = nodenames.index(nodename)
                    i = [(i0 + j) % nn for j in range(nn)]
                    colname2, iname = self.new_column_name(justfn = justfn, chars = chars)
                    col2 = column(colname2,
                                  node = [col.node[i[2]], col.node[i[3]], col.node[i[0]]],
                                  surface = col.surface)
                    # switch connections and neighbours from col to col2 as needed:
                    n3 = col.node[i[3]]
                    n3cols = [c for c in list(col.neighbour) if n3 in c.node]
                    swapcons, swapnbrs = [], []
                    for con in list(col.connection):
                        if con.column[0] in n3cols:
                            con.column[1] = col2
                            swapcons.append(con)
                            swapnbrs.append(con.column[0])
                        elif con.column[1] in n3cols:
                            con.column[0] = col2
                            swapcons.append(con)
                            swapnbrs.append(con.column[1])
                    for con in swapcons:
                        col.connection.remove(con)
                        col2.connection.add(con)
                    for c in swapnbrs:
                        col.neighbour.remove(c)
                        c.neighbour.remove(col)
                        col2.neighbour.add(c)
                        c.neighbour.add(col2)
                    del col.node[i[3]]
                    col.centre = col.centroid
                    self.add_column(col2)
                    self.add_connection(connection([col, col2]))
                    self.setup_block_name_index()
                    self.setup_block_connection_name_index()
                    return True
                except ValueError: return False # node not in column
        return False

    def rename_column(self, oldcolname, newcolname):
        """Renames a column or list of columns."""
        if isinstance(oldcolname, str) and isinstance(newcolname, str):
            oldcolname, newcolname = [oldcolname], [newcolname]
        try:
            for olditem, newitem in zip(oldcolname, newcolname):
                i = self.columnlist.index(self.column[olditem])
                self.columnlist[i].name = newitem
                self.column[newitem] = self.column.pop(olditem)
            self.setup_block_name_index()
            self.setup_block_connection_name_index()
            return True
        except ValueError: return False

    def clear_layers(self):
        """Deletes all layers from the grid."""
        self.layer = {}
        self.layerlist = []

    def add_layer(self, lay = None):
        """Adds layer to the grid. If a layer with the same name
        already exists in the geometry, no new layer is added."""
        if lay is None: lay = layer()
        if lay.name not in self.layer:
            self.layerlist.append(lay)
            self.layer[lay.name] = self.layerlist[-1]

    def delete_layer(self, layername):
        """Deletes a layer from the geometry."""
        layer = self.layer[layername]
        del self.layer[layername]
        self.layerlist.remove(layer)

    def rename_layer(self, oldlayername, newlayername):
        """Renames a layer or list of layers."""
        if isinstance(oldlayername, str) and isinstance(newlayername, str):
            oldlayername, newlayername = [oldlayername], [newlayername]
        try:
            for olditem, newitem in zip(oldlayername, newlayername):
                i = self.layerlist.index(self.layer[olditem])
                self.layerlist[i].name = newitem
                self.layer[newitem] = self.layer.pop(olditem)
            self.setup_block_name_index()
            self.setup_block_connection_name_index()
            return True
        except ValueError: return False

    def add_connection(self, con = None):
        """Adds connection to the grid. If a connection with the same
        names already exists in the geometry, no new connection is added."""
        if con is None: con = connection()
        names = (con.column[0].name, con.column[1].name)
        if names not in self.connection:
            self.connectionlist.append(con)
            self.connection[names] = self.connectionlist[-1]
            self.connectionlist[-1].node = self.connection_nodes(con.column)
            for col in self.connectionlist[-1].column:
                col.connection.add(self.connectionlist[-1])

    def connection_nodes(self, cols):
        """Identifies nodes on the connection between a pair of two columns.
        The node ordering in the first column is preserved.
        """
        connodes = None
        if cols[0].num_nodes > 2: a, b = 0, 1
        elif cols[1].num_nodes > 2: a, b = 1, 0
        else: return None
        for i, n in enumerate(cols[a].node):
            nextn = cols[a].node[(i + 1) % cols[a].num_nodes]
            if n in cols[b].node and nextn in cols[b].node:
                connodes = [n, nextn]
                break
        if a == 0: return connodes
        else: return connodes[::-1] if connodes else connodes

    def delete_connection(self, colnames):
        """Deletes a connection from the geometry."""
        con = self.connection[colnames]
        for col in con.column: col.connection.remove(con)
        del self.connection[colnames]
        self.connectionlist.remove(con)

    def add_well(self, wl = None):
        """Adds well to the geometry. If a well with the specified name
        already exists in the geometry, no new well is added."""
        if wl is None: wl = well()
        if wl.name not in self.well:
            self.welllist.append(wl)
            self.well[wl.name] = self.welllist[-1]

    def delete_well(self, wellname):
        """Deletes a well from the geometry."""
        well = self.well[wellname]
        del self.well[wellname]
        self.welllist.remove(well)

    def delete_orphan_wells(self):
        """Deletes any wells with wellheads not inside the grid."""
        delwells = []
        for well in self.welllist:
            wh = well.pos[0][0:2]
            if self.column_containing_point(wh) is None: delwells.append(well.name)
        for wellname in delwells: self.delete_well(wellname)

    def identify_neighbours(self):
        """Identify neighbour columns"""
        for con in self.connectionlist:
            for i in range(2): con.column[i].neighbour.add(con.column[not i])

    def identify_layer_tops(self):
        """Identifies top elevations of grid layers"""
        self.layerlist[0].top = self.layerlist[0].bottom
        for i, this in enumerate(self.layerlist[1:]):
            above = self.layerlist[i]
            this.top = above.bottom

    def set_default_surface(self):
        """Sets default column surface elevations"""
        ground = self.layerlist[0].bottom
        for col in self.columnlist:
            col.surface = ground
            col.default_surface = True
            col.num_layers = self.num_layers - 1

    def set_column_num_layers(self, col):
        """Sets col.num_layers property according to the layers in the grid."""
        col.num_layers = len([layer for layer in self.layerlist[1:] if
                              layer.bottom < col.surface])

    def column_surface_layer_index(self, col):
        """Returns the index in the layerlist of the surface layer for the
        given column."""
        return self.num_layers - col.num_layers

    def column_surface_layer(self, col):
        """Returns the surface layer for the given column."""
        return self.layerlist[self.column_surface_layer_index(col)]

    def copy_layers_from(self, geo):
        """Copies layer structure from another geometry."""
        self.clear_layers()
        from copy import deepcopy
        for lay in geo.layerlist: self.add_layer(deepcopy(lay))
        for col in self.columnlist: self.set_column_num_layers(col)
        self.setup_block_name_index()
        self.setup_block_connection_name_index()

    def copy_wells_from(self, geo):
        """Copies wells from another geometry."""
        self.well, self.welllist = {}, []
        from copy import deepcopy
        for w in geo.welllist: self.add_well(deepcopy(w))

    def get_min_surface_block_thickness(self):
        """Returns the minimum surface block thickness, and the column name it occurs in."""
        surfcols = [col for col in self.columnlist if col.surface is not None]
        thick = np.array([col.surface - self.column_surface_layer(col).bottom for
                          col in surfcols])
        imin = np.argmin(thick)
        return thick[imin], surfcols[imin].name
    min_surface_block_thickness = property(get_min_surface_block_thickness)

    def columns_in_polygon(self, polygon):
        """Returns a list of all columns with centres inside the specified polygon."""
        return [col for col in self.columnlist if col.in_polygon(polygon)]

    def nodes_in_polygon(self, polygon):
        """Returns a list of all nodes inside the specified polygon."""
        if len(polygon) == 2:
            return [node for node in self.nodelist if in_rectangle(node.pos, polygon)]
        else: return [node for node in self.nodelist if in_polygon(node.pos, polygon)]

    def column_quadtree(self, columns = None):
        """Returns a quadtree structure for searching the grid for columns
        containing particular points.  If the columns parameter is
        specified, a quadtree is returned just for those columns,
        otherwise it is for all columns.
        """
        if columns is None:
            bounds = self.bounds
            columns = self.columnlist
        else: bounds = self.column_bounds(columns)
        return quadtree(bounds, columns)

    def get_node_kdtree(self):
        """Returns a kd-tree structure for searching the grid for particular nodes."""
        from scipy.spatial import cKDTree
        return cKDTree([node.pos for node in self.nodelist])
    node_kdtree = property(get_node_kdtree)

    def node_nearest_to(self, point, kdtree = None):
        """Returns the node nearest to the specified point.  A kd-tree can be
        specified to speed searching- useful if searching for a lot of
        points.
        """
        if isinstance(point, (list, tuple)): point = np.array(point)
        if kdtree:
            r, i = kdtree.query(point)
            return self.nodelist[i]
        else:
            d = np.array([np.linalg.norm(node.pos - point) for node in self.nodelist])
            isort = np.argsort(d)
            return self.nodelist[isort[0]]

    def read_header(self, geo):
        """Reads grid header info from file geo"""
        geo.read_value_line(self.__dict__, 'header')
        self.convention = self._convention
        self.atmosphere_type = self._atmosphere_type
        self.unit_type = self._unit_type
        if self.cntype is not None and self.cntype != 0:
            print('CNTYPE option = %d not supported.' % (self.cntype))

    def read_nodes(self, geo):
        """Reads grid nodes from file geo"""
        line = padstring(geo.readline())
        while line.strip():
            [name, x, y] = geo.parse_string(line, 'node')
            name = name.strip().rjust(self.colname_length)
            pos = np.array([x, y]) * self.unit_scale
            newnode = node(name, pos)
            self.add_node(newnode)
            line = geo.readline()

    def read_columns(self, geo):
        """Reads grid columns from file geo"""
        line = padstring(geo.readline())
        while line.strip():
            [colname, centre_specified,
             nnodes, centrex, centrey] = geo.parse_string(line, 'column')
            colname = colname.strip().rjust(self.colname_length)
            if centre_specified:
                centre = np.array([centrex, centrey]) * self.unit_scale
            else: centre = None
            nodes = []
            for each in range(nnodes):
                [nodename] = geo.read_values('column_node')
                nodename = nodename.strip().rjust(self.colname_length)
                colnode = self.node[nodename]
                nodes.append(colnode)
            self.add_column(column(colname, nodes, centre))
            line = geo.readline()

    def read_connections(self, geo):
        """Reads grid connections from file geo"""
        line = padstring(geo.readline())
        while line.strip():
            names = geo.parse_string(line, 'connection')
            names = [name.strip().rjust(self.colname_length) for name in names]
            cols = [self.column[name] for name in names]
            self.add_connection(connection(cols))
            line = geo.readline()
        self.identify_neighbours()

    def read_layers(self, geo):
        """Reads grid layers from file geo"""
        line = padstring(geo.readline())
        while line.strip():
            name, bottom, centre = geo.parse_string(line, 'layer')
            name = name.strip().rjust(self.layername_length)
            bottom *= self.unit_scale
            newlayer = layer(name, bottom)
            self.add_layer(newlayer)
            if centre: centre *= self.unit_scale
            else:
                nlayers = len(self.layer)
                if nlayers > 1:
                    centre = 0.5 * (newlayer.bottom +
                                    self.layerlist[nlayers - 2].bottom)
                else: centre = newlayer.bottom
            newlayer.centre = centre
            line = geo.readline()
        self.identify_layer_tops()
        self.set_default_surface()

    def read_surface(self, geo):
        """Reads grid surface from file geo"""
        line = padstring(geo.readline())
        while line.strip():
            [name, surface] = geo.parse_string(line, 'surface')
            name = name.strip().rjust(self.colname_length)
            surface *= self.unit_scale
            col = self.column[name]
            col.surface = surface
            self.set_column_num_layers(col)
            line = geo.readline()

    def read_wells(self, geo):
        """Reads grid wells from file geo"""
        line = padstring(geo.readline())
        while line.strip():
            [name, x, y, z] = geo.parse_string(line, 'well')
            p = np.array([x, y, z]) * self.unit_scale
            if name in self.well: self.well[name].pos.append(p)
            else: self.add_well(well(name, [p]))
            line = geo.readline()

    def read(self, filename):
        """Reads MULgraph grid from file"""
        self.empty()
        geo = fixed_format_file(filename, 'rU',
                                mulgrid_format_specification, self.read_function)
        self.read_header(geo)
        if self.type == 'GENER':
            read_fn = {
                'VERTI':self.read_nodes,
                'GRID': self.read_columns,
                'CONNE': self.read_connections,
                'LAYER': self.read_layers,
                'SURFA': self.read_surface,
                'SURF': self.read_surface,
                'WELLS': self.read_wells}
            more = True
            while more:
                line = geo.readline().strip()
                if line:
                    keyword = line[0:5].rstrip()
                    read_fn[keyword](geo)
                else: more = False
            self.setup_block_name_index()
            self.setup_block_connection_name_index()
        else: print('Grid type', self.type, 'not supported.')
        geo.close()
        return self

    def block_surface(self, lay, col):
        """Returns elevation of top of block for given layer and column"""
        if lay.name == self.layerlist[0].name:
            if self.atmosphere_type == 1: return lay.top # atmosphere block
            else: return None
        else:
            if col.surface is None: return lay.top
            else:
                if col.surface < lay.top:
                    if lay.bottom < col.surface:
                        return col.surface # surface layer with surface below layer top
                    else: return None # outside grid
                elif col.surface > self.layerlist[0].top:
                    if lay.name == self.layerlist[1].name:
                        return col.surface # surface layer with surface above layer top
                    else: return lay.top
                else: return lay.top # subsurface layer

    def block_volume(self, lay, col):
        """Returns volume of block at specified layer and column"""
        if lay.name == self.layerlist[0].name:
            if (self.atmosphere_type == 0) and \
               (col.name == self.atmosphere_column_name):
                return self.atmosphere_volume
            elif self.atmosphere_type == 1: return self.atmosphere_volume
            else: return None
        else:
            surf = self.block_surface(lay, col)
            if surf is not None: return (surf - lay.bottom) * col.area
            else: return None

    def block_centre(self, lay, col):
        """Returns centre of block at specified layer and column.  The
        vertical centre is always the layer centre, except for surface
        blocks with column surface lower than the layer top.  For
        surface blocks with column surface higher than the layer top,
        the vertical centre is still the layer centre, to give a
        uniform pressure reference, even though this is not
        geometrically correct.
        """
        if isinstance(lay, str): lay = self.layer[lay]
        if isinstance(col, str): col = self.column[col]
        if lay.name == self.layerlist[0].name:
            if self.atmosphere_type == 1: midelev = lay.centre
            else: return None
        else:
            if (lay.bottom < col.surface <= lay.top):
                midelev = 0.5 * (lay.bottom + col.surface)
            else: 
                if col.surface <= lay.bottom: return None # outside grid
                else: midelev = lay.centre
        return np.array([col.centre[0], col.centre[1], midelev])

    def connection_params(self, con, lay):
        """Returns connection parameters (distances, interface area) between
        specified columns, for the given layer"""
        if con in self.connectionlist:
            sidelength = norm(con.node[0].pos - con.node[1].pos)
            height = min([self.block_surface(lay, c) - lay.bottom for c in con.column])
            area = sidelength * height
            nodeline = [node.pos for node in con.node]
            dist = [norm(line_projection(c.centre, nodeline) - c.centre) for c in con.column]
            return [dist, area]
        else: return [[0.0, 0.0], 0.0] # no connection

    def block_name(self, layername, colname, blockmap = {}):
        """Returns block name from layer and column names, depending on the
        naming convention.  An optional block mapping can be applied.
        """
        if self.convention == 0: name = colname[0:3] + layername[0:2]
        elif self.convention == 1: name = layername[0:3] + colname[0:2]
        else: name = layername[0:2] + colname[0:3]
        blkname = fix_blockname(name)
        if blkname in blockmap: blkname = blockmap[blkname]
        return blkname

    def write(self, filename = ''):
        """Writes a MULgraph grid to file"""
        if filename: self.filename = filename
        if self.filename == '': self.filename = 'geometry.dat'
        geo = fixed_format_file(self.filename, 'w', mulgrid_format_specification)
        self.write_header(geo)
        self.write_nodes(geo)
        self.write_columns(geo)
        self.write_connections(geo)
        self.write_layers(geo)
        if not self.default_surface: self.write_surface(geo)
        if self.num_wells > 0: self.write_wells(geo)
        geo.write('\n')
        geo.close()

    def write_header(self, geo):
        """Writes MULgraph grid header to file"""
        geo.write_value_line(self.__dict__, 'header')

    def write_nodes(self, geo):
        """Writes MULgraph grid nodes to file"""
        geo.write('VERTICES\n')
        for node in self.nodelist:
            name, pos = node.name.ljust(3), node.pos / self.unit_scale
            geo.write_values([name, pos[0], pos[1]], 'node')
        geo.write('\n')

    def write_columns(self, geo):
        """Writes MULgraph grid columns to file"""
        geo.write('GRID\n')
        for col in self.columnlist:
            name = col.name.ljust(3)
            vals = [name, col.centre_specified, col.num_nodes]
            if col.centre_specified: centre = list(col.centre / self.unit_scale)
            else: centre = [None] * 2
            vals += centre
            geo.write_values(vals, 'column')
            for node in col.node: geo.write_values([node.name.ljust(3)], 'column_node')
        geo.write('\n')

    def write_connections(self, geo):
        """Writes MULgraph grid connections to file"""
        geo.write('CONNECTIONS\n')
        for con in self.connectionlist:
            names = [col.name.ljust(3) for col in con.column]
            geo.write_values(names, 'connection')
        geo.write('\n')

    def write_layers(self, geo):
        """Writes MULgraph grid layers to file"""
        geo.write('LAYERS\n')
        for lay in self.layerlist:
            vals = [lay.name.ljust(3), lay.bottom / self.unit_scale, lay.centre / self.unit_scale]
            geo.write_values(vals, 'layer')
        geo.write('\n')

    def write_surface(self, geo):
        """Writes MULgraph grid surface to file"""
        geo.write('SURFA\n')
        for col in [col for col in self.columnlist if not col.default_surface]:
            geo.write_values([col.name.ljust(3), col.surface / self.unit_scale], 'surface')
        geo.write('\n')

    def write_wells(self, geo):
        """Writes MULgraph wells to file"""
        geo.write('WELLS\n')
        for wl in self.welllist:
            for pos in wl.pos:
                vals = [wl.name] + list(pos / self.unit_scale)
                geo.write_values(vals, 'well')
        geo.write('\n')

    def rectangular(self, xblocks, yblocks, zblocks,
                    convention = 0, atmos_type = 2, origin = None,
                    justify = 'r', case = None,  chars = ascii_lowercase,
                    spaces = True):
        """Returns a rectangular MULgraph grid with specified block sizes.
        The arguments are arrays of the block sizes in each dimension
        (x,y,z).  Naming convention, atmosphere type and origin can
        optionally be specified.  The optional justify and case
        parameters specify the format of the character part of the
        block names (whether they are right or left justified, and
        lower or upper case).
        """
        if origin is None: origin = [0., 0., 0.]
        if isinstance(xblocks, (list, tuple)): xblocks = np.array(xblocks)
        if isinstance(yblocks, (list, tuple)): yblocks = np.array(yblocks)
        if isinstance(zblocks, (list, tuple)): zblocks = np.array(zblocks)
        if isinstance(origin, (list, tuple)): origin = np.array(origin)
        grid = mulgrid(type = 'GENER', convention = convention, atmos_type = atmos_type)
        grid.empty()
        xverts = np.array([0.] + np.cumsum(xblocks).tolist()) + origin[0]
        yverts = np.array([0.] + np.cumsum(yblocks).tolist()) + origin[1]
        nxv = len(xverts)
        nxb, nyb = len(xblocks), len(yblocks)
        justfn = [str.rjust, str.ljust][justify == 'l']
        if case is not None:
            casefn = [str.upper, str.lower][case == 'l']
            chars = casefn(chars)
        chars = uniqstring(chars)
        # create nodes:
        num = 1
        y = origin[1]
        for y in yverts:
            for x in xverts:
                name = grid.node_name_from_number(num, justfn, chars, spaces)
                grid.add_node(node(name, np.array([x, y])))
                num += 1
        # create columns:
        num = 1
        for j in range(nyb):
            for i in range(nxb):
                colname = grid.column_name_from_number(num, justfn, chars, spaces)
                colverts = [
                    j * nxv + i + 1,
                    (j + 1) * nxv + i + 1,
                    (j + 1) * nxv + i + 2,
                    j * nxv + i + 2]
                nodenames = [grid.node_name_from_number(v, justfn, chars, spaces)
                             for v in colverts]
                colnodes = [grid.node[name] for name in nodenames]
                grid.add_column(column(colname, colnodes))
                num += 1
        # x-connections:
        for j in range(nyb):
            for i in range(nxb - 1):
                num1, num2 = j * nxb + i + 1, j * nxb + i + 2
                name1 = grid.column_name_from_number(num1, justfn, chars, spaces)
                name2 = grid.column_name_from_number(num2, justfn, chars, spaces)
                grid.add_connection(connection([grid.column[name1], grid.column[name2]]))
        # y-connections:
        for i in range(nxb):
            for j in range(nyb - 1):
                num1, num2 = j * nxb + i + 1, (j + 1) * nxb + i + 1
                name1 = grid.column_name_from_number(num1, justfn, chars, spaces)
                name2 = grid.column_name_from_number(num2, justfn, chars, spaces)
                grid.add_connection(connection([grid.column[name1], grid.column[name2]]))
        # create layers:
        grid.add_layers(zblocks, origin[2], justify, chars, spaces)
        grid.set_default_surface()
        grid.identify_neighbours()
        grid.setup_block_name_index()
        grid.setup_block_connection_name_index()
        return grid

    def add_layers(self, thicknesses, top_elevation = 0, justify = 'r',
                   chars  =  ascii_lowercase, spaces = True):
        """Adds layers of specified thicknesses and top elevation."""
        justfn = [str.rjust, str.ljust][justify == 'l']
        chars = uniqstring(chars)
        num = 0
        self.clear_layers()
        z = top_elevation
        surfacelayername = [' 0', 'atm', 'at'][self.convention]
        self.add_layer(layer(surfacelayername, z, z))
        for thickness in thicknesses:
            z -= thickness
            centre = z + 0.5 * thickness
            name = surfacelayername
            while name == surfacelayername:
                # make sure layer name is different from surface layer name
                num += 1
                name = self.layer_name_from_number(num, justfn, chars, spaces)
            self.add_layer(layer(name, z, centre))
        self.identify_layer_tops()

    def from_gmsh(self, filename, layers, convention = 0, atmos_type = 2,
                  top_elevation = 0, chars = ascii_lowercase, spaces = True):
        """Returns a MULgraph grid constructed from a 2D gmsh grid and the
        specified layer structure."""
        grid = mulgrid(type = 'GENER', convention = convention, atmos_type = atmos_type)
        grid.empty()
        gmsh = open(filename, 'rU')
        line = ''
        chars = uniqstring(chars)
        while not '$Nodes' in line: line = gmsh.readline()
        num_nodes = int(gmsh.readline().strip())
        for i in range(num_nodes):
            items = gmsh.readline().strip().split(' ')
            name, x, y = items[0], float(items[1]), float(items[2])
            name = self.node_name_from_number(int(name), chars = chars, spaces = spaces)
            grid.add_node(node(name, np.array([x, y])))
        while not '$Elements' in line: line = gmsh.readline()
        num_elements = int(gmsh.readline().strip())
        for i in range(num_elements):
            items = gmsh.readline().strip().split(' ')
            element_type = int(items[1])
            if element_type in [2, 3]: # triangle or quadrilateral
                name = items[0]
                name = self.column_name_from_number(int(name), chars = chars,
                                                    spaces = True)
                ntags = int(items[2])
                colnodenumbers = [int(item) for item in items[3 + ntags:]]
                colnodenames = [[self.node_name_from_number(nodeno,
                                                            chars = chars,
                                                            spaces = spaces),
                                 nodeno][convention > 0] for nodeno in colnodenumbers]
                colnodes = [grid.node[v] for v in colnodenames]
                grid.add_column(column(name, colnodes))
        gmsh.close()
        for con in grid.missing_connections: grid.add_connection(con)
        grid.delete_orphans()
        grid.add_layers(layers, top_elevation, chars, spaces)
        grid.set_default_surface()
        grid.identify_neighbours()
        grid.setup_block_name_index()
        grid.setup_block_connection_name_index()
        return grid

    def translate(self, shift, wells = False):
        """Translates a grid by specified shift.  If wells is True, they
        will also be translated."""
        if isinstance(shift, (list, tuple)): shift = np.array(shift)
        for node in self.nodelist: node.pos += shift[0:2]
        for col in self.columnlist:
            col.centre += shift[0:2]
            if col.surface is not None: col.surface += shift[2]
        for layer in self.layerlist: layer.translate(shift[2])
        if wells:
            for well in self.welllist:
                for pos in well.pos: pos += shift

    def rotate(self, angle, centre = None, wells = False):
        """Rotates grid horizontally by specified angle (degrees clockwise).
        If centre is not specified, the centre of the grid is used.
        If wells is True, they will also be rotated."""
        if centre is not None:
            if isinstance(centre, (list, tuple)): centre = np.array(centre)
            c = centre
        else: c = self.centre
        R = linear_trans2().rotation(angle, c)
        for node in self.nodelist: node.pos = R(node.pos)
        for col in self.columnlist: col.centre = R(col.centre)
        if wells:
            for well in self.welllist:
                for pos in well.pos: pos[0:2] = R(pos[0:2])

    def get_missing_connections(self):
        """Returns a set of connections for columns that have shared faces but
        don't have a connection defined between them."""
        missing = set([])
        for node in self.nodelist:
            nodecols = list(node.column)
            for i, coli in enumerate(nodecols[0:-1]):
                for colj in nodecols[i + 1:]:
                    if coli.is_against(colj) and not self.connects(coli, colj):
                        if coli.name <= colj.name: mincol, maxcol = coli, colj
                        else: mincol, maxcol = colj, coli
                        missing.add((mincol.name, maxcol.name))
        return set([connection([self.column[colname] for
                                colname in m]) for m in missing])
    missing_connections = property(get_missing_connections)

    def get_extra_connections(self):
        """Returns a set of pairs of column names defined between columns that
        aren't against each other."""
        extra = set([])
        for con in self.connectionlist:
            if not con.column[0].is_against(con.column[1]):
                extra.add(tuple([col.name for col in con.column]))
        return extra
    extra_connections = property(get_extra_connections)

    def get_orphans(self):
        """Returns a set of 'orphaned' nodes, i.e. nodes that do not belong to
        any column."""
        return set([node for node in self.nodelist if len(node.column) == 0])
    orphans = property(get_orphans)

    def delete_orphans(self):
        """Deletes any orphaned nodes."""
        for node in self.orphans: self.delete_node(node.name)

    def get_bad_columns(self):
        """Returns a set of columns that do not contain their own centres."""
        return set([col for col in self.columnlist if not col.contains_point(col.centre)])
    bad_columns = property(get_bad_columns)

    def get_bad_layers(self):
        """Returns a set of layers that do not contain their own centres."""
        return set([layer for layer in self.layerlist[1:] if
                    not layer.bottom <= layer.centre <= layer.top])
    bad_layers = property(get_bad_layers)

    def check(self, fix = False, silent = False):
        """Checks a grid for errors, and optionally fixes them.  Errors checked for are:
        - missing connections
        - extra connections
        - orphaned nodes
        - columns and layers that do not contain their own centres.
        Returns True if no errors were found, and False otherwise.
        If silent is True, there is no printout."""
        ok = True
        mc = self.missing_connections
        if len(mc) > 0:
            ok = False
            if not silent: print('Missing connections:', list(mc))
            if fix:
                for c in mc: self.add_connection(c)
                if not silent: print('Missing connections fixed.')
        ec = self.extra_connections
        if len(ec) > 0:
            ok = False
            if not silent: print('Extra connections:', list(ec))
            if fix:
                for c in ec: self.delete_connection(c)
                if not silent: print('Extra connections fixed.')
        orphans = self.orphans
        if len(orphans) > 0:
            ok = False
            if not silent: print('Orphaned nodes:', list(orphans))
            if fix:
                self.delete_orphans()
                if not silent: print('Orphaned nodes deleted.')
        bc = self.bad_columns
        if len(bc) > 0:
            ok = False
            if not silent: print('Bad columns:', list(bc))
            if fix:
                for c in bc: c.centre = sum([n.pos for n in c.node]) / c.num_nodes
                if not silent: print('Columns fixed.')
        bl = self.bad_layers
        if len(bl) > 0:
            ok = False
            if not silent: print('Bad layers:', list(bl))
            if fix:
                for layer in bl:
                    layer.bottom = min(layer.bottom, layer.top)
                    layer.top = max(layer.bottom, layer.top)
                    layer.centre = 0.5 * (layer.bottom + layer.top)
                if not silent: print('Layers fixed.')
        if ok and not silent: print('No problems found.')
        return ok

    def column_values_to_block(self, x):
        """Takes an array of values for each column and extends it into an
        array of values for each block."""
        blkval = np.zeros(self.num_blocks, float64)
        colindex = self.column_index
        for i, blk in enumerate(self.block_name_list):
            colname = self.column_name(blk)
            if colname in self.column:
                ci = colindex[colname]
                blkval[i] = x[ci]
        return blkval

    def column_containing_point(self, pos, columns = None, guess = None,
                                bounds = None, qtree = None):
        """Returns column containing the specified horizontal position (or
        None if not found).  If the columns parameter is specified,
        search only within the given list of columns.  A starting
        guess of the column can also be optionally provided, in which
        case that column and (if necessary) its neighbours will be
        searched first.  A bounding polygon (or rectangle)for
        searching within can also optionally be supplied- this can,
        for example, be specified as the boundary polygon of the grid.
        A quadtree for searching the columns can also optionally be
        specified.
        """
        target = None
        if bounds is not None:
            if len(bounds) == 2: inbounds = in_rectangle(pos, bounds)
            else: inbounds = in_polygon(pos, bounds)
        else: inbounds = True
        if inbounds:
            if columns is None: searchcols = self.columnlist
            else: searchcols = columns
            donecols = set([])
            if guess is not None:
                if guess.contains_point(pos): return guess
                else: # search neighbours of guess, sorted by distance from pos:
                    donecols.add(guess)
                    from copy import copy
                    nbrcols = list(copy(guess.neighbour))
                    nearnbrcols = [col for col in nbrcols if
                                   col.near_point(pos) and col in searchcols]
                    sortindex = np.argsort([norm(col.centre - pos) for col in nearnbrcols])
                    for i in sortindex:
                        if nearnbrcols[i].contains_point(pos): return nearnbrcols[i]
                    donecols.update(set(nearnbrcols))
            # guess was no good- do full search on remaining columns:
            if qtree: return qtree.search(pos)
            else:
                nearcols = list(set([col for col in searchcols if
                                     col.near_point(pos)]) - donecols)
                sortindex = np.argsort([norm(col.centre - pos) for col in nearcols])
                for i in sortindex:
                    if nearcols[i].contains_point(pos):
                        target = nearcols[i]
                        break
        return target

    def layer_containing_elevation(self, z):
        """Returns layer containing the specified vertical elevation (or None
        if not found)."""
        target = None
        for layer in self.layerlist[1:]:
            if layer.contains_elevation(z):
                target = layer
                break
        return target

    def column_mapping(self, geo):
        """Returns a dictionary mapping each column name in a geometry object
         geo to the name of the nearest column in self.  If the SciPy
         library is available, a KDTree structure is used to speed
         searching.
        """
        if self.atmosphere_type == geo.atmosphere_type == 0:
            mapping = {geo.atmosphere_column_name:self.atmosphere_column_name}
        else: mapping = {}
        try:
            from scipy.spatial import cKDTree
            kdtree = cKDTree([col.centre for col in self.columnlist])
            def closest_col(col):
                r, i = kdtree.query(col.centre)
                return self.columnlist[i]
        except ImportError: # if don't have SciPy installed:
            def closest_col(col):
                coldist = np.array([norm(selfcol.centre - col.centre) for
                                    selfcol in self.columnlist])
                return self.columnlist[np.argmin(coldist)]
        for col in geo.columnlist: mapping[col.name] = closest_col(col).name
        return mapping

    def layer_mapping(self, geo):
        """Returns a dictionary mapping each layer name in a geometry object
        geo to the name of the nearest layer in self."""
        mapping = {geo.layerlist[0].name: self.layerlist[0].name}  # surface mapped to surface
        for layer in geo.layerlist[1:]:
            laydist = np.array([abs(selflay.centre - layer.centre) for
                                selflay in self.layerlist[1:]])
            # (1 added for surface layer, omitted from search):
            closest = self.layerlist[1 + np.argmin(laydist)]
            mapping[layer.name] = closest.name
        return mapping

    def block_mapping(self, geo, column_mapping = False):
        """Returns a dictionary mapping each block name in a geometry object
        geo to the name of the nearest block in self.  Columns are
        given priority over layers, i.e. first the nearest column is
        found, then the nearest layer for blocks in that column.  The
        associated column mapping can also optionally be returned.
        """
        mapping = {}
        col_mapping = self.column_mapping(geo)
        layer_mapping = self.layer_mapping(geo)
        for dest in geo.block_name_list:
            destcol, destlayer = geo.column_name(dest), geo.layer_name(dest)
            sourcecol, sourcelayer = col_mapping[destcol], layer_mapping[destlayer]
            if destlayer == geo.layerlist[0].name:
                sourcelayer = self.layerlist[0].name # atmosphere layer
                if self.atmosphere_type == 0:
                    sourcecol = self.atmosphere_column_name
            else:
                # if source block is above surface in column, use
                # first layer below surface instead:
                if self.column[sourcecol].surface <= self.layer[sourcelayer].bottom:
                    sourcelayer = self.column_surface_layer(self.column[sourcecol]).name
            mapping[dest] = self.block_name(sourcelayer, sourcecol)
        if column_mapping: return (mapping, col_mapping)
        else: return mapping

    def block_name_containing_point(self, pos, qtree = None, blockmap = {}):
        """Returns name of grid block containing 3D point (or None if the
        point is outside the grid)."""
        blkname = None
        col = self.column_containing_point(pos[0:2], qtree = qtree)
        if col:
            if self.layerlist[0].bottom < pos[2] <= col.surface:
                layer = self.layerlist[1] 
            else: layer = self.layer_containing_elevation(pos[2])
            if layer:
                if (col.surface > layer.bottom):
                    blkname = self.block_name(layer.name, col.name, blockmap)
        return blkname

    def block_contains_point(self, blockname, pos):
        """Returns True if the block with specified name contains the
        specified 3D point."""
        result = False
        colname = self.column_name(blockname)
        if colname in self.column:
            col = self.column[colname]
            layname = self.layer_name(blockname)
            if layname in self.layer:
                lay = self.layer[layname]
                if col.surface > lay.bottom:
                    if lay.contains_elevation(pos[2]):
                        result = col.contains_point(pos[0:2])
        return result

    def column_track(self, line):
        """Returns a list of tuples of (column,entrypoint,exitpoint)
        representing the horizontal track traversed by the specified
        line through the grid.  Line is a tuple of two 2D arrays.  The
        resulting list is ordered by distance from the start of the
        line.
        """

        def furthest_intersection(poly, line):
            """Returns furthest intersection point between line and poly."""
            pts, inds = line_polygon_intersections(poly, line,
                                                   bound_line = (True, False),
                                                   indices = True)
            if pts:
                d = np.array([np.linalg.norm(intpt - line[0]) for intpt in pts])
                i = np.argmax(d)
                return pts[i], inds[i]
            else: return None, None

        def find_track_start(line):
            """Finds starting point for track- an arbitrary point on the line that is inside
            the grid.  If the start point of the line is inside the grid, that is used;
            otherwise, a recursive bisection technique is used to find a point."""
            col, start_type = None, None
            for endpt, name in zip(line, ['start', 'end']):
                pos, col, start_type = endpt, self.column_containing_point(endpt), name
                if col: break
            if not col: # line ends are both outside the grid:
                start_type = 'mid'
                max_levels = 7

                def find_start(line, level = 0):
                    midpt = 0.5 * (line[0] + line[1])
                    col = self.column_containing_point(midpt)
                    if col: return midpt, col
                    else:
                        if level <= max_levels:
                            line0, line1 = [line[0], midpt], [midpt, line[1]]
                            pos, col = find_start(line0, level + 1)
                            if col: return pos, col
                            else:
                                pos, col = find_start(line1, level + 1)
                                if col: return pos, col
                                else: return None, None
                        else: return None, None

                pos, col = find_start(line)
            return pos, col, start_type

        def next_corner_column(col, pos, more, cols):
            """If the line has hit a node, determine a new column containing that node,
            not already visited."""
            node_tol = 1.e-12
            nextcol = None
            nearnodes = [n for n in col.node if np.linalg.norm(n.pos - pos) < node_tol]
            if nearnodes: # hit a node
                nearnode = nearnodes[0]
                nearcols = nearnode.column - cols
                if nearcols: nextcol = nearcols.pop()
                else: more = False
            return nextcol, more

        def next_neighbour_column(col, more, cols):
            """Determine a new neighbour column not already visited."""
            nbrs = col.neighbour - cols
            if nbrs: return nbrs.pop(), more
            else: return None, False

        def find_track_segment(linesegment, pos, col):
            """Finds track segment starting from the specified position and column."""
            track = []
            cols, more, inpos = set(), True, pos
            colnbr, nextcol = col.neighbourlist, None
            lined = np.linalg.norm(linesegment[1] - linesegment[0])
            while more:
                cols.add(col)
                outpos, ind = furthest_intersection(col.polygon, linesegment)
                if outpos is not None:
                    d = np.linalg.norm(outpos - linesegment[0])
                    if d >= lined: # gone past end of line
                        outpos = linesegment[1]
                        more = False
                    if np.linalg.norm(outpos - inpos) > 0.:
                        track.append(tuple([col, inpos, outpos]))
                    if more: # find next column
                        inpos = outpos
                        nextcol = colnbr[ind]
                        if nextcol:
                            if nextcol in cols:
                                nextcol, more = next_corner_column(col, outpos, more, cols)
                                if nextcol is None:
                                    nextcol, more = next_neighbour_column(col, more, cols)
                                    nbr_base_col = col
                        else: nextcol, more = next_corner_column(col, outpos, more, cols)
                else:
                    nextcol, more = next_neighbour_column(nbr_base_col, more, cols)
                col = nextcol
                if col: colnbr = col.neighbourlist
                else: more = False

            return track

        def reverse_track(track): return [tuple([tk[0], tk[2], tk[1]]) for tk in track][::-1]

        pos, col, start_type = find_track_start(line)
        if pos is not None and col:
            if start_type == 'start':
                track = find_track_segment(line, pos, col)
            elif start_type == 'end':
                track = find_track_segment(line[::-1], pos, col)
                track = reverse_track(track)
            else:
                track1 = find_track_segment([pos, line[0]], pos, col)
                track2 = find_track_segment([pos, line[1]], pos, col)
                # remove arbitrary starting point from middle of track, and join:
                midtk = tuple([track1[0][0], track1[0][2], track2[0][2]])
                track = reverse_track(track1)[:-1] + [midtk] + track2[1:]
            return track
        else: return []

    def layer_plot_wells(self, plt, ax, layer, wells, well_names,
                         hide_wells_outside, wellcolour, welllinewidth, wellname_bottom):
        """Draws wells on a layer plot, given the plot and axes.  For
        documentation of the parameters, see layer_plot()."""
        if wells is True: wells = self.welllist
        elif wells is False or wells is None: wells = []
        elif isinstance(wells, list):
            if isinstance(wells[0], str): wells = [self.well[name] for name in wells]
        if well_names is True: well_names = wells
        elif well_names is None or well_names is False: well_names = []
        elif isinstance(well_names, list):
            if isinstance(well_names[0], str): well_names = [self.well[name] for name in well_names]
        for well in wells:
            if layer is None: hitlayer, show_well = False, True # surface plot
            else:
                toppos = well.elevation_pos(layer.top)
                bottompos = well.elevation_pos(layer.bottom)
                hitlayer = toppos is not None
                show_well = hitlayer or not hide_wells_outside
            if show_well:
                plt.plot(well.head[0], well.head[1], 'o', color = wellcolour)
                wpos = [well.pos_coordinate(i) for i in range(3)]
                if hitlayer:
                    above = np.where(wpos[2] > toppos[2])
                    abovepos = [list(wpos[i][above]) + [toppos[i]] for i in range(2)]
                    if bottompos is not None: # well passes through layer
                        inside = np.where((toppos[2] >= wpos[2]) & (wpos[2] >= bottompos[2]))
                        insidepos = [[toppos[i]] + list(wpos[i][inside]) + \
                                     [bottompos[i]] for i in range(2)]
                        below = np.where(wpos[2] < bottompos[2])
                        belowpos = [[bottompos[i]] + list(wpos[i][below]) for i in range(2)]
                    else: # well stops in layer
                        inside = np.where(toppos[2] >= wpos[2])
                        insidepos = [[toppos[i]] + list(wpos[i][inside]) for i in range(2)]
                        belowpos =  np.array([])
                    if abovepos:
                        plt.plot(abovepos[0], abovepos[1], ':', color = wellcolour,
                                 linewidth = welllinewidth)
                    if insidepos:
                        plt.plot(insidepos[0], insidepos[1], '-', color = wellcolour,
                                 linewidth = welllinewidth)
                    if belowpos:
                        plt.plot(belowpos[0], belowpos[1], ':', color = wellcolour,
                                 linewidth = welllinewidth)
                else:
                    if layer is None: welllinestyle = '-'
                    else: welllinestyle =':'
                    plt.plot(wpos[0], wpos[1], welllinestyle, color = wellcolour,
                             linewidth = welllinewidth)
                if well in well_names:
                    if wellname_bottom: namepos = well.bottom
                    else: namepos = well.head
                    ax.text(namepos[0], namepos[1], well.name, clip_on = True,
                            horizontalalignment = 'center')

    def setup_rocktype_plot(self, grid, vals, colourmap, allrocks, rockgroup):
        """Sets up plotted values and colourbar parameters for rock types on
        layer or slice plots."""
        from matplotlib import cm
        if colourmap is None: colourmap = 'jet'
        if allrocks:
            num_shown_rocks = grid.num_rocktypes
            rocknames = [rt.name for rt in grid.rocktypelist]
        else:
            shown_rocks = list(set(vals))
            shown_rocks.sort()
            num_shown_rocks = len(shown_rocks)
            rockmap = dict(zip(shown_rocks, range(num_shown_rocks)))
            vals = [rockmap[val] for val in vals]
            rocknames = [grid.rocktypelist[irock].name for irock in shown_rocks]
        if rockgroup:
            if isinstance(rockgroup, str):
                rockgroup = [i for i, c in enumerate(rockgroup) if c != '*']
            def namegroup(name):
                return ''.join([c if i in rockgroup else '*' for
                                i, c in enumerate(name)])
            valgroup = dict([(i, namegroup(name)) for
                             i, name in enumerate(rocknames)])
            rocknames = list(set(valgroup.values()))
            rocknames.sort()
            num_shown_rocks = len(rocknames)
            rockmap = dict(zip(rocknames, range(num_shown_rocks)))
            vals = [rockmap[valgroup[val]] for val in vals]
        colourmap = cm.get_cmap(colourmap, num_shown_rocks)
        colourbar_limits = (0, num_shown_rocks)
        return vals, rocknames, colourmap, colourbar_limits

    def plot_colourbar(self, plt, col, scalelabel, rocktypes, rocknames):
        """Draws colour bar on a layer or slice plot."""
        cbar = plt.colorbar(col)
        cbar.set_label(scalelabel)
        if rocktypes:
            cbar.set_ticks([i + 0.5 for i in range(len(rocknames))])
            cbar.set_ticklabels(rocknames)
            cbar.ax.invert_yaxis() # to get in same top-down order as in the data file

    def plot_flows(self, plt, X, Y, U, V, flow_variable_name, flow_unit,
                   flow_scale, flow_scale_pos, arrow_width,
                   connection_flows = False):
        """Draws flows (and a key) on a layer or slice plot."""
        if len(X) > 0:
            maxq = max([np.linalg.norm(np.array([u, v])) for u, v in zip(U, V)])
            if flow_scale is None:
                if maxq > 0.:
                    from math import log10
                    flow_scale = 10 ** (int(round(log10(maxq)))) # order of magnitude
                else: flow_scale = 1.e-9
            flow_scale_factor = flow_scale * 10.
            if connection_flows:
                pivot = 'tail'
                key_str = flow_variable_name + ' = ' + str(flow_scale) + ' ' + flow_unit
            else:
                pivot = 'middle'
                key_str = flow_variable_name + '/area = ' + str(flow_scale) + \
                          ' ' + flow_unit + '/$m^2$'
            Q = plt.quiver(X, Y, U, V, units = 'width', pivot = pivot, scale = flow_scale_factor,
                           scale_units = 'width', width = arrow_width)
            qk = plt.quiverkey(Q, flow_scale_pos[0], flow_scale_pos[1], flow_scale, key_str)

    def layer_plot(self, layer = 0, variable = None, variable_name = None,
                   unit = None, column_names = None, node_names = None,
                   column_centres = None, nodes = None, colourmap = None,
                   linewidth = 0.2, linecolour = 'black', aspect = 'equal',
                   plt = None, subplot = 111, title = None,
                   xlabel = 'x (m)', ylabel = 'y (m)',
                   contours = False, contour_label_format = '%3.0f',
                   contour_grid_divisions = (100,100),
                   connections = None, colourbar_limits = None, plot_limits = None,
                   wells  =  None, well_names = True,
                   hide_wells_outside = False, wellcolour = 'blue',
                   welllinewidth = 1.0, wellname_bottom = True,
                   rocktypes = None, allrocks = False, rockgroup = None,
                   flow = None, grid = None, flux_matrix = None,
                   flow_variable_name = None, flow_unit = None, flow_scale = None,
                   flow_scale_pos = (0.5, 0.02),
                   flow_arrow_width = None, connection_flows = False,
                   blockmap = {}, block_names = None):
        """Produces a layer plot of a Mulgraph grid, shaded by the specified
        variable (an array of values for each block).  A unit string
        can be specified for annotation.  Column names, node names,
        column centres and nodes can be optionally superimposed, and
        the colour map, linewidth, aspect ratio, colour-bar limits and
        plot limits specified.  If no variable is specified, only the
        grid is drawn, without shading. If an elevation (float) is
        given instead of a layer name, the layer containing that
        elevation is plotted.  If layer is set to None, then the
        ground surface is plotted (i.e. the surface layer for each
        column).
        """
        import matplotlib
        if plt is None:
            import matplotlib.pyplot as plt
            loneplot = True
        else: loneplot = False
        matplotlib.rcParams.update({'mathtext.default': 'regular',
                                    'figure.figsize':(12, 9)})
        ax = plt.subplot(subplot, aspect = aspect)
        if isinstance(layer, (float, int)):
            layer_elev = layer
            layer = self.layer_containing_elevation(float(layer))
            if layer:
                default_title = 'layer ' + layer.name + ' (elevation ' + \
                                ("%4.0f" % float(layer_elev)).strip() + ' m)'
            else: raise Exception("Layer elevation out of range in layer_plot()")
        elif layer is None: default_title = 'surface layer'
        elif layer in self.layerlist: default_title = 'layer ' + layer.name
        elif layer in self.layer:
            layer = self.layer[layer]
            default_title = 'layer ' + layer.name
        else: raise Exception("Unknown layer in layer_plot()")
        if variable is not None:
            if len(variable) == self.num_columns:
                variable = self.column_values_to_block(variable)
        if variable_name: varname = variable_name
        else: varname = 'Value'
        if column_names:
            if not isinstance(column_names, list):
                column_names = self.column.keys()
        else: column_names = []
        if block_names:
            if block_names == True:
                block_names = [
                    blockmap[blk] if blk in blockmap else blk for
                    blk in self.block_name_list]
        else: block_names = []
        if node_names:
            if not isinstance(node_names, list):
                node_names = self.node.keys()
        else: node_names = []
        if column_centres:
            if not isinstance(column_centres, list):
                column_centres = self.column.keys()
        else: column_centres = []
        if nodes:
            if not isinstance(nodes, list): nodes = self.node.keys()
        else: nodes = []
        verts, vals = [], []
        if not isinstance(contours, bool): contours = list(contours)
        xc, yc = [], []
        if connections is not None:
            c = np.abs(self.connection_angle_cosine)
            ithreshold = np.where(c > connections)[0]
            from matplotlib.colors import colorConverter
            for i in ithreshold:
                colc = [col.centre for col in self.connectionlist[i].column]
                plt.plot([p[0] for p in colc],
                         [p[1] for p in colc],
                         color = colorConverter.to_rgb(str(1. - c[i])))
        if rocktypes:
            variable, varname = rocktypes.get_rocktype_indices(self, blockmap), 'Rock type'
        if flow is not None:
            if flow_variable_name is None: flow_variable_name = 'Flow'
            if flow_unit is None: flow_unit = 'units'
            if grid is None:
                from t2grids import t2grid
                grid = t2grid().fromgeo(self)
            if not grid.connection_centres_defined:
                grid.calculate_connection_centres(self)
            if not connection_flows:
                if flux_matrix is None: flux_matrix = grid.flux_matrix(self)
                blkflow = flux_matrix * flow
                blkflow = blkflow.reshape((self.num_underground_blocks, 3))
            natm = self.num_atmosphere_blocks
            U, V = [], []

        layercols = []
        for col in self.columnlist:
            if layer is None:
                layername = self.column_surface_layer(col).name
            else: layername = layer.name
            blkname = self.block_name(layername, col.name)
            if blkname in self.block_name_index:
                layercols.append(col)
                xc.append(col.centre[0])
                yc.append(col.centre[1])
                if variable is not None:
                    val = variable[self.block_name_index[blkname]]
                else: val = 0
                vals.append(val)
                verts.append(tuple([tuple([p for p in n.pos]) for n in col.node]))
                if col.name in column_names:
                    ax.text(col.centre[0], col.centre[1], col.name,
                            clip_on = True, horizontalalignment = 'center')
                if col.name in column_centres:
                    ax.text(col.centre[0], col.centre[1], ' + ', color = 'red', clip_on = True,
                            horizontalalignment = 'center', verticalalignment = 'center')
                mapped_blkname = blockmap[blkname] if blkname in blockmap else blkname
                if mapped_blkname in block_names:
                    ax.text(col.centre[0],  col.centre[1],  mapped_blkname,
                            color = 'red',  clip_on = True,
                            horizontalalignment = 'center',  verticalalignment = 'center')
                if flow is not None and not connection_flows:
                    blkindex = self.block_name_index[blkname] - natm
                    q = blkflow[blkindex]
                    U.append(q[0])
                    V.append(q[1])

        for node in [self.node[name] for name in node_names]:
                ax.text(node.pos[0], node.pos[1], node.name,
                        clip_on = True, horizontalalignment = 'center')
        for node in [self.node[name] for name in nodes]:
                ax.text(node.pos[0], node.pos[1], ' + ', color = 'red', clip_on = True,
                        horizontalalignment = 'center', verticalalignment = 'center')
        import matplotlib.collections as collections
        if variable is not None: facecolors = None
        else: facecolors = []
        if rocktypes: vals, rocknames, colourmap, colourbar_limits = \
                self.setup_rocktype_plot(rocktypes, vals, colourmap, allrocks, rockgroup)
        else: rocknames, rocktypes = None, None
        if block_names:
            if block_names == True: block_names = [blockmap[blk] if blk in blockmap
                                                   else blk for blk in self.block_name_list]
        else: block_names = []
        col = collections.PolyCollection(verts, cmap = colourmap,
                                         linewidth = linewidth,
                                         facecolors = facecolors, edgecolors = linecolour)
        if variable is not None: col.set_array(np.array(vals))
        if colourbar_limits is not None: col.norm.vmin, col.norm.vmax = tuple(colourbar_limits)
        ax.add_collection(col)
        if plot_limits is not None:
            plt.xlim(plot_limits[0])
            plt.ylim(plot_limits[1])
        else: ax.autoscale_view()
        if contours is not False:
            from matplotlib.mlab import griddata
            valc = np.array(vals)
            bds = self.bounds
            xgrid = np.linspace(bds[0][0], bds[1][0], contour_grid_divisions[0])
            ygrid = np.linspace(bds[0][1], bds[1][1], contour_grid_divisions[1])
            valgrid = griddata(xc, yc, valc, xgrid, ygrid, interp = 'linear')
            if isinstance(contours, list): cvals = contours
            else: cvals = False
            CS = plt.contour(xgrid, ygrid, valgrid, cvals, colors = 'k')
            if contour_label_format is not None:
                plt.clabel(CS, inline = 1, fmt = contour_label_format)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        scalelabel = varname
        if unit: scalelabel  += ' (' + unit + ')'
        if variable is not None:
            self.plot_colourbar(plt, col, scalelabel, rocktypes, rocknames)
            default_title = varname + ' in ' + default_title
        self.layer_plot_wells(plt, ax, layer, wells, well_names,
                              hide_wells_outside, wellcolour, welllinewidth,
                              wellname_bottom)
        if flow is not None:
            if connection_flows:
                xflow, yflow = [], []
                for geocon in self.connectionlist:
                    conblknames = tuple([self.block_name(layername, col.name) for
                                         col in geocon.column])
                    if conblknames in grid.connection:
                        con = grid.connection[conblknames]
                        xflow.append(con.midpoint[0])
                        yflow.append(con.midpoint[1])
                        con_index = self.block_connection_name_index[conblknames]
                        U.append(-flow[con_index] * con.normal[0])
                        V.append(-flow[con_index] * con.normal[1])
            else:
                ishow = [ind for ind, col in enumerate(layercols) if col.num_nodes >= 3]
                xflow, yflow = np.array(xc)[ishow], np.array(yc)[ishow]
                U, V = np.array(U)[ishow],  np.array(V)[ishow]
            self.plot_flows(plt, xflow, yflow, U, V, flow_variable_name,
                            flow_unit, flow_scale, flow_scale_pos,
                            flow_arrow_width, connection_flows)
        if title is None: title = default_title
        plt.title(title)
        if loneplot: plt.show()

    def slice_plot_wells(self, plt, ax, linespec, line, wells, well_names, hide_wells_outside,
                         wellcolour, welllinewidth, wellname_bottom):
        """Draws wells on a layer plot, given the plot and axes.  For
        documentation of the parameters, see slice_plot()."""
        if wells is True: wells = self.welllist
        elif wells is False or wells is None: wells = []
        elif isinstance(wells, list):
            if isinstance(wells[0], str):
                wells = [self.well[name] for name in wells]
        if well_names is True: well_names = wells
        elif well_names is None or well_names is False:
            well_names = []
        elif isinstance(well_names, list):
            if isinstance(well_names[0], str):
                well_names = [self.well[name] for name in well_names]
        if len(wells) > 0:
            if hide_wells_outside is False:
                def show_well(well): return True
            elif isinstance(hide_wells_outside, (float,int)):
                def show_well(well):
                    hpos = [pos[:2] for pos in well.pos]
                    return polyline_line_distance(hpos, line) <= hide_wells_outside
            else:
                raise Exception('slice_plot_wells() error: ' + \
                                'unrecognised value for parameter hide_wells_outside')

            def slice_project(pos):
                """Returns 2-D projection of a 3-D point onto the slice plane"""
                hppos = line_projection(pos[:2], line)
                return np.array([np.linalg.norm(hppos - line[0]), pos[2]])

            def slice_translate(pos):
                if linespec == 'x': pos  += np.array([line[0][0], 0.])
                elif linespec == 'y': pos += np.array([line[0][1], 0.])
                return pos

            for well in wells:
                if show_well(well):
                    pwellhead = slice_translate(slice_project(well.head))
                    plt.plot(pwellhead[0], pwellhead[1], 'o', color = wellcolour)
                    if hide_wells_outside is False:
                        wpos = slice_translate(np.array([slice_project(pos) for
                                                         pos in well.pos]))
                        plt.plot(wpos[:, 0], wpos[:, 1], '-',
                                 color = wellcolour, linewidth = welllinewidth)
                    else: # draw well sections outside as dotted lines
                        wellsections = {True: [], False: []}
                        top = well.head
                        dtop = point_line_distance(top[:2], line)
                        topinside = dtop <= hide_wells_outside
                        wellsections[topinside].append([top])
                        for bot in well.pos[1:]:
                            dbot = point_line_distance(bot[:2], line)
                            botinside = dbot <= hide_wells_outside
                            if botinside == topinside:
                                wellsections[topinside][-1].append(bot)
                            else:
                                try:
                                    xi = (hide_wells_outside - dtop) / (dbot - dtop)
                                except ZeroDivisionError: xi = 1.0
                                cross = (1. - xi) * top + xi * bot
                                wellsections[topinside][-1].append(cross)
                                wellsections[botinside].append([cross])
                            top, dtop, topinside = bot, dbot, botinside
                        linetype = {True:'-', False:':'}
                        for inside, sections in wellsections.items():
                            for section in sections:
                                wpos = slice_translate(np.array([slice_project(pos) for
                                                                 pos in section]))
                                plt.plot(wpos[:, 0], wpos[:, 1], linetype[inside],
                                         color = wellcolour, linewidth = welllinewidth)
                    if well in well_names:
                        if wellname_bottom:
                            namepos, namealign = well.bottom, 'top'
                        else:
                            namepos, namealign = well.head, 'bottom'
                        nameposp = slice_translate(slice_project(namepos))
                        ax.text(nameposp[0], nameposp[1], well.name, clip_on = True,
                                horizontalalignment = 'center', verticalalignment = namealign)

    def slice_plot(self, line = None, variable = None, variable_name = None,
                   unit = None, block_names = None, colourmap = None, linewidth = 0.2,
                   linecolour = 'black', aspect = 'auto', plt = None, subplot = 111,
                   title = None, xlabel = None, ylabel = 'elevation (m)',
                   contours = False, contour_label_format = '%3.0f',
                   contour_grid_divisions = (100,100), colourbar_limits = None,
                   plot_limits = None, column_axis = False, layer_axis = False,
                   wells = None, well_names = True, hide_wells_outside = False,
                   wellcolour = 'blue', welllinewidth = 1.0, wellname_bottom = False,
                   rocktypes = None, allrocks = False, rockgroup = None, flow = None,
                   grid = None, flux_matrix = None, flow_variable_name = None,
                   flow_unit = None, flow_scale = None, flow_scale_pos = (0.5, 0.02),
                   flow_arrow_width = None, connection_flows = False, blockmap = {}):
        """Produces a vertical slice plot of a Mulgraph grid, shaded by the
       specified variable (an array of values for each block).  A unit
       string can be specified for annotation.  Block names can be
       optionally superimposed, and the colour map, linewidth, aspect
       ratio, colour-bar limits and plot limits specified.  If no
       variable is specified, only the grid is drawn, without shading.
       If no line is specified, a slice through the grid bounds is
       made (bottom left to top right).  If a string 'x' or 'y' is
       passed in instead of a line, a plot is made through the centre
       of the grid along the x- or y-axes, and the coordinate along
       the slice represents the actual x- or y- coordinate.  If a
       northing (float, in degrees) is passed instead of a line, a
       plot is made through the centre along the specified northing
       direction.
        """
        if line is None:
            l = self.bounds
            default_title = 'vertical slice across grid bounds'
        elif isinstance(line, str):
            axislines = {'x':[np.array([self.bounds[0][0], self.centre[1]]),
                              np.array([self.bounds[1][0], self.centre[1]])],
                         'y':[np.array([self.centre[0], self.bounds[0][1]]),
                              np.array([self.centre[0], self.bounds[1][1]])]}
            if line in axislines:
                l = axislines[line]
                default_title = 'vertical slice along ' + line + ' axis'
            else:
                l = self.bounds
                default_title = 'vertical slice across grid bounds'
        elif isinstance(line, (float, int)):
            r = 0.5 * norm(self.bounds[1] - self.bounds[0])
            from math import radians, cos, sin
            theta = radians(line)
            d = r * np.array([sin(theta), cos(theta)])
            l = [self.centre - d, self.centre + d]
            default_title = 'vertical slice ' + ("%3.0f" % float(line)).strip() + \
                            '$^\circ$N'
        else:
            l = line
            default_title = 'vertical slice from (' + ("%7.0f" % l[0][0]).strip() + \
                            ', ' + ("%7.0f" % l[0][1]).strip() + ') to (' + \
                            ("%7.0f" % l[1][0]).strip() + ', ' + \
                            ("%7.0f" % l[1][1]).strip() + ')'
        if norm(l[1] - l[0]) > 0.0:
            import matplotlib
            if plt is None:
                import matplotlib.pyplot as plt
                loneplot = True
            else: loneplot = False
            matplotlib.rcParams.update({'mathtext.default': 'regular',
                                        'figure.figsize':(12, 9)})
            ax = plt.subplot(subplot, aspect = aspect)
            if variable is not None:
                if len(variable) == self.num_columns:
                    variable = self.column_values_to_block(variable)
            if variable_name: varname = variable_name
            else: varname = 'Value'
            if rocktypes:
                variable = rocktypes.get_rocktype_indices(self, blockmap)
                varname = 'Rock type'
            if block_names:
                if block_names == True:
                    block_names = [blockmap[blk] if blk in blockmap else blk for
                                   blk in self.block_name_list]
            else: block_names = []

            track = self.column_track(l)
            if track:
                if xlabel is None:
                    if line == 'x': xlabel = 'x (m)'
                    elif line == 'y': xlabel = 'y (m)'
                    else: xlabel = 'distance (m)'
                ax.set_xlabel(xlabel)
                if column_axis: colnames, colcentres = [], []
                verts, vals = [], []
                if not isinstance(contours, bool): contours = list(contours)
                xc, yc = [], []
                if flow is not None:
                    if flow_variable_name is None: flow_variable_name = 'Flow'
                    if flow_unit is None: flow_unit = 'units'
                    if grid is None:
                        from t2grids import t2grid
                        grid = t2grid().fromgeo(self)
                    if not grid.connection_centres_defined:
                        grid.calculate_connection_centres(self)
                    if not connection_flows:
                        if flux_matrix is None:
                            flux_matrix = grid.flux_matrix(self)
                        blkflow = flux_matrix * flow
                        blkflow = blkflow.reshape((self.num_underground_blocks, 3))
                    natm = self.num_atmosphere_blocks
                    U, V = [], []
                    slice_dirn = (l[1] - l[0]).T
                    # normal vector in slice direction:
                    slice_dirn /= np.linalg.norm(slice_dirn)

                ind, ishow, sliceblocks, coldist = 0, [], set(), {}
                for trackitem in track:
                    col, points = trackitem[0], trackitem[1:]
                    inpoint = points[0]
                    if len(points) > 1: outpoint = points[1]
                    else: outpoint = inpoint
                    if line == 'x':
                        din, dout = inpoint[0], outpoint[0]
                    elif line == 'y':
                        din, dout = inpoint[1], outpoint[1]
                    else:
                        din, dout = norm(inpoint - l[0]), norm(outpoint - l[0])
                    coldist[col.name] = (din, dout)
                    dcol = 0.5 * (din + dout)
                    if column_axis:
                        colnames.append(col.name); colcentres.append(dcol)
                    for lay in self.layerlist[1:]:
                        if col.surface > lay.bottom:
                            blkname = self.block_name(lay.name, col.name)
                            sliceblocks.add(blkname)
                            if variable is not None:
                                val = variable[self.block_name_index[blkname]]
                            else: val = 0
                            vals.append(val)
                            top = self.block_surface(lay, col)
                            centre = self.block_centre(lay, col)
                            verts.append(((din, lay.bottom),
                                          (din, top),
                                          (dout, top),
                                          (dout, lay.bottom)))
                            mapped_blkname = blockmap[blkname] if \
                                             blkname in blockmap else blkname
                            if mapped_blkname in block_names:
                                ax.text(dcol, centre[2], mapped_blkname,
                                        clip_on = True, horizontalalignment = 'center')
                            xc.append(dcol); yc.append(centre[2])
                            if flow is not None and not connection_flows:
                                blkindex = self.block_name_index[blkname] - natm
                                q = blkflow[blkindex]
                                qslice = np.dot(slice_dirn, q[:2])
                                U.append(qslice)
                                V.append(q[2])
                            if col.num_nodes >= 3: ishow.append(ind)
                            ind += 1

                import matplotlib.collections as collections
                if variable is not None: facecolors = None
                else: facecolors = []
                if rocktypes: vals, rocknames, colourmap, colourbar_limits = \
                        self.setup_rocktype_plot(rocktypes, vals, colourmap,
                                                 allrocks, rockgroup)
                else: rocknames, rocktypes = None, None
                col = collections.PolyCollection(verts, cmap = colourmap,
                                                 linewidth = linewidth,
                                                 facecolors = facecolors,
                                                 edgecolors = linecolour)
                if variable is not None:
                    col.set_array(np.array(vals))
                if colourbar_limits is not None:
                    col.norm.vmin, col.norm.vmax = tuple(colourbar_limits)
                ax.add_collection(col)
                if contours != False:
                    from matplotlib.mlab import griddata
                    valc = np.array(vals)
                    bds = ((np.min(xc), np.min(yc)), (np.max(xc), np.max(yc)))
                    xgrid = np.linspace(bds[0][0], bds[1][0], contour_grid_divisions[0])
                    ygrid = np.linspace(bds[0][1], bds[1][1], contour_grid_divisions[1])
                    valgrid = griddata(xc, yc, valc, xgrid, ygrid, interp = 'linear')
                    if isinstance(contours, list): cvals = contours
                    else: cvals = False
                    CS = plt.contour(xgrid, ygrid, valgrid, cvals, colors = 'k')
                    if contour_label_format is not None:
                        plt.clabel(CS,  inline = 1, fmt = contour_label_format)
                ax.set_ylabel(ylabel)
                scalelabel = varname
                if unit: scalelabel += ' (' + unit + ')'
                if variable is not None:
                    self.plot_colourbar(plt, col, scalelabel, rocktypes, rocknames)
                    default_title = varname + ' in ' + default_title
                if column_axis:
                    ax.set_xticks(colcentres)
                    ax.set_xticklabels(colnames)
                    ax.set_xlabel('column')
                if layer_axis:
                    ax.set_yticks([lay.centre for lay in self.layerlist])
                    ax.set_yticklabels([lay.name for lay in self.layerlist])
                    ax.set_ylabel('layer')
                self.slice_plot_wells(plt, ax, line, l, wells, well_names,
                                      hide_wells_outside, wellcolour,
                                      welllinewidth, wellname_bottom)
                if plot_limits is not None:
                    ax.set_xlim(plot_limits[0]); ax.set_ylim(plot_limits[1])
                else: ax.autoscale_view()
                if flow is not None:
                    if connection_flows:
                        xflow, yflow = [], []
                        sliceblocks = sliceblocks | set(self.block_name_list[:natm])
                        slicecons = [con for con in self.block_connection_name_list if
                                     all([blkname in sliceblocks for blkname in con])]
                        for conblknames in slicecons:
                            con = grid.connection[conblknames]
                            conlays = [self.layer_name(blkname) for blkname in conblknames]
                            concols = [self.column_name(blkname) for blkname in conblknames]
                            if conlays[0] == conlays[1]: # horizontal
                                segments = [coldist[concol] for concol in concols]
                                if segments[0][0] > segments[1][0]: d = segments[1][1]
                                else: d = segments[0][1]
                            else: # vertical
                                iblk = [i for i, blkname in enumerate(conblknames) if
                                        self.block_name_index[blkname] >= natm][0]
                                col = concols[iblk]
                                segment = coldist[col]
                                d = 0.5 * (segment[0] + segment[1])
                            xflow.append(d)
                            yflow.append(con.midpoint[2])
                            con_index = self.block_connection_name_index[conblknames]
                            nflow = -flow[con_index] * con.normal
                            qslice = np.dot(slice_dirn, nflow[:2])
                            U.append(qslice)
                            V.append(nflow[2])
                    else:
                        xflow, yflow = np.array(xc)[ishow], np.array(yc)[ishow]
                        U, V = np.array(U)[ishow], np.array(V)[ishow]
                    self.plot_flows(plt, xflow, yflow, U, V,
                                    flow_variable_name, flow_unit, flow_scale,
                                    flow_scale_pos, flow_arrow_width, connection_flows)
                if title is None: title = default_title
                plt.title(title)
                if loneplot: plt.show()
            else: print('Slice', str(line), 'does not intersect the grid.')

    def line_values(self, start, end, variable, divisions = 100,
                    coordinate = False, qtree = None):
        """Gets values of variable along specified line through geometry.
        Returns two arrays for distance along line (or specified
        coordinate) and value at each position.
        """
        if isinstance(start, (list, tuple)): start = np.array(start)
        if isinstance(end, (list, tuple)): end = np.array(end)
        x, y = [], []
        line_length = norm(end - start)
        if line_length > 0.0:
            for i in range(divisions + 1):
                xi = float(i) / divisions
                pos = (1. - xi) * start + xi * end
                dist = xi * line_length
                blkname = self.block_name_containing_point(pos, qtree = qtree)
                if blkname:
                    if coordinate != None: x.append(pos[coordinate])
                    else: x.append(dist)
                    y.append(variable[self.block_name_index[blkname]])
        return np.array(x), np.array(y)

    def polyline_values(self, polyline, variable, divisions = 100,
                        coordinate = False, qtree = None):
        """Gets values of a variable along a specified polyline, returning two
        arrays for distance along the polyline and value."""
        x, y = [], []
        for i in range(len(polyline) - 1):
            start, end = polyline[i], polyline[i + 1]
            xi, yi = self.line_values(start, end, variable, divisions, coordinate, qtree = qtree)
            if i > 0:
                xi = xi[1:]; yi = yi[1:]
            if not coordinate:
                if len(x) > 0: xi += x[-1]  # add end distance from last segment
            x += list(xi)
            y += list(yi)
        return np.array(x), np.array(y)

    def well_values(self, well_name, variable, divisions = 1,
                    elevation = False, deviations = False,
                    qtree = None, extend = False):
        """Gets values of a variable down a specified well, returning distance
        down the well (or elevation) and value.  Vertical coordinates
        can be taken from the nodes of the well deviations, or from
        the grid geometry layer centres (if deviations is False).  If
        extend is True, the well trace is extended to the bottom of
        the model.
        """
        if elevation: coordinate = 2  # return coordinate 2 (i.e. z)
        else: coordinate = False
        if well_name in self.well:
            well = self.well[well_name]
            if deviations:
                from copy import copy
                polyline = copy(well.pos)
                grid_bottom = self.layerlist[-1].bottom
                if extend and well.bottom[2] > grid_bottom:
                    polyline.append(well.elevation_pos(grid_bottom, extend = True))
            else:
                polyline = []
                for layer in self.layerlist:
                    p = well.elevation_pos(layer.centre, extend = extend)
                    if p is not None: polyline.append(p)
            return self.polyline_values(polyline, variable, divisions, coordinate, qtree = qtree)
        else: return None

    def column_values(self, col, variable, depth = False):
        """Gets values of a variable down a specified column in the grid,
        returning elevation (or depth) and value.
        """
        if isinstance(col, str):
            if col in self.column: col = self.column[col]
            else: return None
        itop = self.column_surface_layer_index(col)
        # atmosphere block:
        if self.atmosphere_type == 0:
            blks = [self.block_name(self.layerlist[0].name,
                                    self.atmosphere_column_name)]
        elif self.atmosphere_type ==1:
            blks = [self.block_name(self.layerlist[0].name, col.name)]
        else: blks = []
        # indices of subsurface blocks:
        lays = self.layerlist[itop:]
        blks += [self.block_name(lay.name, col.name) for lay in lays]
        blkindex = np.array([self.block_name_index[blk] for blk in blks])
        val = variable[blkindex]
        z = np.array([self.layer[self.layer_name(blk)].centre for blk in blks])
        if depth:
            return self.layer[self.layer_name(blks[0])].top - z, val
        else: return z, val

    def line_plot(self, start = None, end = None, variable = None,
                  variable_name = None, unit = None, divisions = 100,
                  plt = None, subplot = 111, title = '', xlabel = 'distance (m)'):
        """Produces a line plot of the specified variable through a Mulgraph grid."""
        if (start is None) or (end is None):
            [start, end] = self.bounds
            default_title = 'line plot across grid bounds'
        else:
            if isinstance(start, (list, tuple)): start = np.array(start)
            if isinstance(end, (list, tuple)): end = np.array(end)
            default_title = 'line plot from (' + ("%7.0f" % start[0]).strip() + \
                            ',' + ("%7.0f" % start[1]).strip() + ',' + \
                            ("%7.0f" % start[2]).strip() + \
                            ') to (' + ("%7.0f" % end[0]).strip() + \
                            ',' + ("%7.0f" % end[1]).strip() + \
                            ',' + ("%7.0f" % end[2]).strip() + ')'
        x, y = self.line_values(start, end, variable, divisions)
        import matplotlib
        if plt is None:
            import matplotlib.pyplot as plt
            loneplot = True
        else: loneplot = False
        matplotlib.rcParams.update({'mathtext.default': 'regular',
                                    'figure.figsize':(12, 9)})
        plt.subplot(subplot)
        if variable is not None:
            if len(variable) == self.num_columns:
                variable = self.column_values_to_block(variable)
        if variable_name: varname = variable_name
        else: varname = 'Value'
        plt.plot(x, y)
        plt.xlabel(xlabel)
        ylabel = varname
        if unit: ylabel +=  ' (' + unit + ')'
        plt.ylabel(ylabel)
        default_title += ' of ' + varname
        if title is None: title = default_title
        plt.title(title)
        if loneplot: plt.show()

    def optimize(self, nodenames = None, connection_angle_weight = 1.0,
                 column_aspect_weight = 0.0, column_skewness_weight = 0.0,
                 pest = False):
        """Adjusts positions of specified nodes to optimize grid.  If
        nodenames list is not specified, all node positions are
        optimized.  Grid quality can be defined as a combination of
        connection angle cosine, column aspect ratio and column
        skewness.  Increasing the weight for any of these increases
        its importance in the evaluation of grid quality.  Note that
        an error will result if the connection angle weight and either
        of the other weights is set to zero- in this case there are
        not enough constraints to fit the parameters.  If pest is set
        to True, the PEST parameter estimation software is used to
        perform the optimization.
        """
        if nodenames is None: nodenames = self.node.keys()
        # identify which columns are affected:
        colnames = [col.name for col in self.columnlist if
                    (set(nodenames) & set([node.name for node in col.node]))]
        for colname in colnames:
            self.column[colname].centre_specified = 0
        if connection_angle_weight > 0.0: # identify which connections are affected:
            cons = [con for con in self.connectionlist if
                    (set(col.name for col in con.column) & set(colnames))]
        if pest:
            original_filename = self.filename
            gridfilename = 'gpestmesh.dat'
            self.write(gridfilename)
            obsgroups = ['angle', 'aspect', 'skew']
            obsweight = {}
            nobs = 0
            if connection_angle_weight > 0.0:
                obsweight['angle'] = connection_angle_weight
                nobs += len(cons)
            if column_aspect_weight > 0.0:
                obsweight['aspect'] = column_aspect_weight
                nobs += len(colnames)
            if column_skewness_weight > 0.0:
                obsweight['skew'] = column_skewness_weight
                nobs += len(colnames)
            def write_pest_control_file():
                pst = open('pestmesh.pst', 'w')
                pst.write('\n'.join([
                            'pcf', '* control data', 'restart estimation',
                            str(2 * len(nodenames)) + ' ' + str(nobs) + ' 1 0 ' + str(len(obsweight)),
                            '    1     1 single point   1   0   0',
                            '5.0   2.0   0.3  0.03    10',
                            '3.0   3.0 0.001  0',
                            '0.1',
                            '30  0.01     3     3  0.01     3',
                            '1     1     1',
                            '* parameter groups',
                            'pos           absolute 0.01  0.0  switch  2.0 parabolic',
                            '* parameter data\n']))
                for name in nodenames:
                    parname = 'node_' + name.strip() + '_'
                    for i in range(2):
                        pst.write(parname + str(i) + ' none relative ' + \
                                  '%12.3f' % self.node[name].pos[i] + ' ' + \
                                  '%12.3f' % self.bounds[0][i] + \
                                  ' ' + '%12.3f' % self.bounds[1][i] + ' pos 1.0 0.0 1\n')
                pst.write('* observation groups\n')
                for group in obsweight: pst.write(group + '\n')
                pst.write('* observation data\n')
                for group in obsgroups:
                    if group in obsweight:
                        if group == 'angle':
                            n = len(cons)
                            target = 0.0
                        else:
                            n = len(colnames)
                            target = 1.0
                        for i in range(n):
                            pst.write(group + str(i) + ' %5.2f' % target + \
                                      ' %5.2f' % obsweight[group] + ' ' + group + '\n')
                pst.write('\n'.join([
                            '* model command line', 'python pestmesh_model.py',
                            '* model input/output', 'pestmesh.tpl  pestmesh.in',
                            'pestmesh.ins  pestmesh.out', '* prior information\n']))
                pst.close()
            def write_pest_model_file():
                mod = open('pestmesh_model.py', 'w')
                mod.write('\n'.join([
                        "from mulgrids import *",
                        "geo = mulgrid('" + gridfilename.strip() + "')",
                        "nodenames = np.load('pestmesh_nodes.npy')",
                        "nnodes = len(nodenames)",
                        "dat = np.loadtxt('pestmesh.in').reshape((nnodes,2))",
                        "for pos,nodename in zip(dat,nodenames):",
                        "    geo.node[nodename].pos = pos",
                        "colnames = np.load('pestmesh_columns.npy')",
                        "for colname in colnames:",
                        "    geo.column[colname].centre = geo.column[colname].centroid",
                        "result = []\n"]))
                if 'angle' in obsweight:
                    mod.write("connames = np.load('pestmesh_connections.npy')\n")
                    mod.write("result += [geo.connection[tuple(conname)].angle_cosine" + \
                              " for conname in connames]\n")
                if 'aspect' in obsweight:
                    mod.write("result += [geo.column[colname].side_ratio for" + \
                              " colname in colnames]\n")
                if 'skew' in obsweight:
                    mod.write("result += [geo.column[colname].angle_ratio for" + \
                              " colname in colnames]\n")
                mod.write('\n'.join([
                            "f = open('pestmesh.out', 'w')",
                            "for r in result: f.write('%20.5f\\n'%r)",
                            "f.close()"]))
                mod.close()
            def write_pest_templates():
                tpl = open('pestmesh.tpl', 'w')
                tpl.write("ptf $\n")
                for name in nodenames:
                    parname = 'node_' + name.strip() + '_'
                    for i in range(2):
                        tpl.write("$" + '%18s' % (parname + str(i)) + "$\n")
                tpl.close()
                ins = open('pestmesh.ins', 'w')
                ins.write("pif #\n")
                for group in obsgroups:
                    if group in obsweight:
                        if group == 'angle': n = len(cons)
                        else: n = len(colnames)
                        for i in range(n): ins.write("l1 [" + group + str(i) + "]1:20\n")
                ins.close()
            np.save('pestmesh_nodes.npy', np.array(nodenames))
            np.save('pestmesh_columns.npy', np.array(colnames))
            if 'angle' in obsweight:
                np.save('pestmesh_connections.npy',
                        np.array([[col.name for col in con.column] for con in cons]))
            write_pest_control_file()
            write_pest_templates()
            write_pest_model_file()
            from os import system
            system('pest pestmesh.pst')
            dat = np.loadtxt('pestmesh.in').reshape((len(nodenames), 2))
            for pos, nodename in zip(dat, nodenames): self.node[nodename].pos = pos
            for colname in colnames:
                self.column[colname].centre = self.column[colname].centroid
                self.column[colname].get_area()
            self.filename = original_filename
        else:
            from scipy.optimize import leastsq
            num_nodes = len(nodenames)
            def update_grid(x):
                xpos = x.reshape((num_nodes, 2))
                for nodename, pos in zip(nodenames, xpos):
                    self.node[nodename].pos = pos
                for colname in colnames:
                    self.column[colname].centre = self.column[colname].centroid
            def f(x):
                update_grid(x)
                result = []
                if connection_angle_weight:
                    result += [connection_angle_weight * con.angle_cosine for con in cons]
                if column_aspect_weight:
                    result += [column_aspect_weight *
                               (self.column[colname].side_ratio - 1.) for
                               colname in colnames]
                if column_skewness_weight:
                    result += [column_skewness_weight *
                               (self.column[colname].angle_ratio - 1.) for
                               colname in colnames]
                return np.array(result)
            x0 = np.array([self.node[nodename].pos for
                           nodename in nodenames]).reshape(2 * num_nodes)
            x1, success = leastsq(f, x0)
            if success > 4:
                raise Exception('scipy leastsq() optimization routine did not converge.')
            update_grid(x1)
            for colname in colnames: self.column[colname].get_area()

    def connection_with_nodes(self, nodes):
        """Returns a connection, if one exists, containing the specified two nodes."""
        for con in self.connectionlist:
            if all([node in con.node for node in nodes]): return con
        return None

    def nodes_in_columns(self, columns):
        """Returns a list of all nodes in the specified columns."""
        nodes = set([])
        for col in columns: nodes = nodes | set(col.node)
        return list(nodes)

    def column_boundary_nodes(self, columns):
        """Returns an ordered list of the nodes on the outer boundary of the
        group of specified columns."""
        nodes = self.nodes_in_columns(columns)
        blacklist_connections, blacklist_nodes = [], []
        def next_bdy_node(n):
            for col in [c for c in n.column if c in columns and c.num_nodes > 2]:
                i = col.node.index(n)
                n2 = col.node[(i + 1) % col.num_nodes]
                if n2 not in blacklist_nodes:
                    con = self.connection_with_nodes([n, n2])
                    if not con or any([c.num_nodes <= 2 for c in con.column]):
                        return n2
                    else:
                        if not (con in blacklist_connections) and \
                                not all([(c in columns) for c in con.column]):
                            return n2
            return None
        # look for a starting node along the left-hand edge of the selection (this avoids
        # picking up any interior boundaries):
        startnode = None
        xmin = bounds_of_points([node.pos for node in nodes])[0][0]
        leftnodes = [node for node in nodes if node.pos[0] == xmin]
        for node in leftnodes:
            nextnode = next_bdy_node(node)
            if nextnode:
                startnode = node
                break
        if startnode:
            bdynodes = []
            node = startnode
            back = False
            while not back:
                bdynodes.append(node)
                node = next_bdy_node(node)
                if node is None:
                    raise Exception('Could not detect column boundary nodes.')
                else:
                    back = node.name == startnode.name
                    if (node in bdynodes) and not back : # loop in boundary
                        nodei = bdynodes.index(node)
                        nnodes = len(bdynodes)
                        loopcount = nnodes - nodei - 1
                        for i in range(loopcount):
                            n1, n2 = bdynodes[-2],  bdynodes[-1]
                            lastnode = bdynodes.pop()
                            con = self.connection_with_nodes([n1, n2])
                            if con: blacklist_connections.append(con)
                            else: blacklist_nodes.append(lastnode)
                        node = bdynodes.pop()
            return bdynodes
        else: return []
    def get_boundary_nodes(self): return self.column_boundary_nodes(self.columnlist)
    boundary_nodes = property(get_boundary_nodes)

    def get_boundary_polygon(self):
        """Returns the simplest polygon representing the boundary of the grid."""
        return simplify_polygon([node.pos for node in self.boundary_nodes])
    boundary_polygon = property(get_boundary_polygon)

    def get_boundary_columns(self):
        """Returns a set of columns on the outer boundary of the grid- those
        columns that contain at least two boundary nodes.
        """
        bdynodes = set(self.boundary_nodes)
        return set([col for col in self.columnlist if
                    len(set(col.node) & bdynodes) >= 2])
    boundary_columns = property(get_boundary_columns)

    def grid2d_horizontal(self):
        """Returns 2D horizontal nodes and elements for the grid."""
        z = self.layerlist[0].top
        nodes = []
        node_index = {}
        for i, node in enumerate(self.nodelist):
            nodes.append(np.array(list(node.pos) + [z]))
            node_index[node.name] = i
        elts = []
        for col in self.columnlist:
            elt = [node_index[node.name] for node in col.node]
            elts.append(elt)
        return nodes, elts

    def grid2d_vertical(self, line):
        """Returns 2D vertical slice nodes and elements for the grid. Any
        specified column surface elevations are not taken into
        account."""
        nodes, elts = [], []
        if line is None:
            l = self.bounds
        elif isinstance(line, str):
            axislines = {'x':[np.array([self.bounds[0][0], self.centre[1]]),
                              np.array([self.bounds[1][0], self.centre[1]])],
                         'y':[np.array([self.centre[0], self.bounds[0][1]]),
                              np.array([self.centre[0], self.bounds[1][1]])]}
            if line in axislines: l = axislines[line]
            else: l = self.bounds
        elif isinstance(line, (float, int)):
            r = 0.5 * norm(self.bounds[1] - self.bounds[0])
            from math import radians, cos, sin
            theta = radians(line)
            d = r * np.array([sin(theta), cos(theta)])
            l = [self.centre - d, self.centre + d]
        else: l = line
        if norm(l[1] - l[0]) > 0.0:
            track = self.column_track(l)
            if track:
                node_d = []
                for itrack, trackitem in enumerate(track):
                    points = trackitem[1:]
                    inpoint = points[0]
                    if len(points) > 1: outpoint = points[1]
                    else: outpoint = inpoint
                    if line == 'x': din, dout = inpoint[0], outpoint[0]
                    elif line == 'y': din, dout = inpoint[1], outpoint[1]
                    else: din, dout = norm(inpoint - l[0]), norm(outpoint - l[0])
                    if itrack == 0: node_d.append(din)
                    node_d.append(dout)
                for lay in self.layerlist:
                    for d in node_d:
                        nodes.append(np.array([d, lay.bottom, 0.]))
                n = len(track)
                for ilay, lay in enumerate(self.layerlist[1:]):
                    for icol in range(n):
                        off = ilay * (n + 1)
                        elt = [off + icol + n + 1, off + icol + n + 2,
                               off + icol + 1, off + icol]
                        elts.append(elt)
        return nodes, elts

    def grid3d(self, surface_snap):
        """Returns 3D nodes and elements for the grid. The surface_snap
        parameter is a tolerance determining how close column
        elevations have to be to be considered 'equal'.
        """
        def surf(col):
            return col.surface if col.surface is not None else \
                self.layerlist[0].bottom
        # Create subsurface nodes:
        node3d = []
        node_index = {}
        index = 0
        for lay in self.layerlist[1:]:
            for node in self.nodelist:
                if any([surf(col) > lay.bottom for col in node.column]):
                    pos3d = np.array(list(node.pos) + [lay.bottom])
                    node3d.append(pos3d)
                    node_index[lay.name, node.name] = index
                    index += 1

        def unique_values(a, tol):
            # Returns unique values in array of floats, within given
            # tolerance, and a list of groups of the original values
            # corresponding to each unique value.
            asort = a.copy()
            asort.sort()
            asort = np.array(list(asort) + [asort[-1] + 2. * tol])
            d = np.diff(asort)
            uniques = asort[np.where(d > tol)]
            nearest = np.array([np.argmin(np.abs(uniques - x)) for x in a])
            groups = [list(np.where(nearest == i)[0]) for i in range(len(uniques))]
            return uniques, groups

        # Create surface nodes:
        surface_node_index = {}
        for node in self.nodelist:
            if node.column:
                col_names = [col.name for col in node.column]
                col_surf = np.array([surf(self.column[col]) for col in col_names])
                unique_elevations, columns = unique_values(col_surf, surface_snap)
                for z, col_group in zip(unique_elevations, columns):
                    create_new = True
                    layerdiff = np.array([abs(z - lay.bottom) for
                                          lay in self.layerlist])
                    imin = np.argmin(layerdiff)
                    if layerdiff[imin] <= surface_snap:
                        lay = self.layerlist[imin]
                        if (lay.name, node.name) in node_index:
                            create_new = False
                        else: node_index[lay.name, node.name] = index
                        sindex = node_index[lay.name, node.name]
                    if create_new:
                        pos3d = np.array(list(node.pos) + [z])
                        node3d.append(pos3d)
                        sindex = index
                        index += 1
                    for i in col_group:
                        surface_node_index[col_names[i], node.name] = sindex
        # Create elements:
        elt3d = []
        atmlayer = self.layerlist[0]
        for ilayer, lay in enumerate(self.layerlist[1:]):
            for col in [c for c in self.columnlist if
                        surf(c) > lay.bottom]:
                elt = []
                for node in col.node:
                    elt.append(node_index[lay.name, node.name])
                if ilayer == self.column_surface_layer_index(col) - 1: # top block
                    for node in col.node:
                        elt.append(surface_node_index[col.name, node.name])
                else:
                    above_layer = self.layerlist[ilayer]
                    for node in col.node:
                        elt.append(node_index[above_layer.name, node.name])
                elt3d.append(elt)
        return node3d, elt3d

    def meshio_grid(self, surface_snap = 0.1, dimension = 3, slice = None):
        """Returns mesh in meshio (points, cells) format. If dimension = 3,
        the full 3D mesh is returned. If dimension = 2, the 2D
        horizontal mesh is returned, unless a slice is specified, in
        which case a 2D vertical slice mesh is returned.
        """
        if dimension == 3:
            nodes, elts = self.grid3d(surface_snap)
            cell_types = {
                6: 'wedge',
                8: 'hexahedron'
            }
        elif dimension == 2:
            if slice is None: nodes, elts = self.grid2d_horizontal()
            else: nodes, elts = self.grid2d_vertical(slice)
            cell_types = {
                3: 'triangle',
                4: 'quad'
            }
        else:
            raise Exception("Unrecognised dimension (%d)" % dimension)
        points = np.array(nodes)
        cells = {}
        for elt in elts:
            n = len(elt)
            if n in cell_types:
                cell_type = cell_types[n]
                if cell_type in cells: cells[cell_type].append(elt)
                else: cells[cell_type] = [elt]
            else:
                raise Exception("Unrecognised cell type (nodes: %d)" % n)
        for cell_type in cells.keys():
            cells[cell_type] = np.array(cells[cell_type])
        return points, cells

    def get_vtk_grid(self, arrays = {}, surface_snap = 0.1):
        """Returns a vtkUnstructuredGrid object (for visualisation with VTK)
        corresponding to the grid in 3D.  VTK data arrays may
        optionally be added.
        """
        from vtk import vtkUnstructuredGrid, vtkPoints, vtkIdList, \
            vtkWedge, vtkHexahedron, vtkPentagonalPrism, \
            vtkHexagonalPrism, vtkConvexPointSet
        node3d, elt3d = self.grid3d(surface_snap)
        # Construct the VTK grid:
        grid = vtkUnstructuredGrid()
        pts = vtkPoints()
        pts.SetNumberOfPoints(len(node3d))
        for i, node in enumerate(node3d): pts.SetPoint(i, node)
        grid.SetPoints(pts)
        # Create and add cells:
        celltype = {
            6: vtkWedge,
            8: vtkHexahedron,
            10: vtkPentagonalPrism,
            12: vtkHexagonalPrism}
        for elt in elt3d:
            n = len(elt)
            if n in celltype:
                cell = celltype[n]()
                for i, j in enumerate(elt):
                    cell.GetPointIds().SetId(i, j)
            else:
                cell = vtkConvexPointSet()
                for i, j in enumerate(elt):
                    cell.GetPointIds().InsertId(i, j)
            grid.InsertNextCell(cell.GetCellType(), cell.GetPointIds())
        for array_type, array_dict in arrays.items():
            sortedkeys = sorted(array_dict.keys())
            if array_type=='Block':
                for key in sortedkeys: grid.GetCellData().AddArray(array_dict[key])
            elif array_type=='Node':
                for key in sortedkeys: grid.GetPointData().AddArray(array_dict[key])
        return grid

    def get_vtk_data(self, blockmap = {}):
        """Returns a dictionary of VTK data arrays from the grid (layer and
        column indices (zero-based), column areas, block numbers and
        volumes for each block.
        """
        from vtk import vtkFloatArray, vtkIntArray, vtkCharArray
        arrays = {'Block': {
            'Name': vtkCharArray(),
            'Layer index': vtkIntArray(),
            'Column index': vtkIntArray(),
            'Column area': vtkFloatArray(),
            'Column elevation': vtkFloatArray(),
            'Block number': vtkIntArray(),
            'Volume': vtkFloatArray()},
                  'Node': {}}
        # SetTupleValue() was changed to SetTypedTuple() in VTK 7.1:
        if hasattr(vtkCharArray, 'SetTupleValue'):
            arrays['Block']['Name'].SetTypedTuple = arrays['Block']['Name'].SetTupleValue
        nele = self.num_underground_blocks
        string_properties = ['Name']
        string_length = 5
        array_length = {'Block': nele, 'Node': 0}
        for array_type, array_dict in arrays.items():
            for name, array in array_dict.items():
                array.SetName(name)
                if name in string_properties:
                    array.SetNumberOfComponents(string_length)
                    array.SetNumberOfTuples(array_length[array_type])
                else:
                    array.SetNumberOfValues(array_length[array_type])
                    array.SetNumberOfComponents(1)
        layerindex = self.layer_index
        colindex = self.column_index
        for iblk, blockname in enumerate(self.block_name_list[self.num_atmosphere_blocks:]):
            layname, colname = self.layer_name(blockname), self.column_name(blockname)
            lay, col = self.layer[layname], self.column[colname]
            mapped_name = blockmap[blockname] if blockname in blockmap else blockname
            arrays['Block']['Name'].SetTypedTuple(iblk, mapped_name)
            arrays['Block']['Layer index'].SetValue(iblk, layerindex[layname])
            arrays['Block']['Column index'].SetValue(iblk, colindex[colname])
            arrays['Block']['Column area'].SetValue(iblk, col.area)
            arrays['Block']['Column elevation'].SetValue(iblk, col.surface)
            arrays['Block']['Block number'].SetValue(iblk, self.block_name_index[blockname] + 1)
            arrays['Block']['Volume'].SetValue(iblk, self.block_volume(lay, col))
        return arrays

    def filename_base(self, filename = ''):
        """Returns base of filename (with extension removed).  If specified
        filename is blank, the geometry filename property is used; if
        this is also blank, a default is used.
        """
        from os.path import splitext
        default_filename = 'geometry.dat'
        if filename == '':
            if self.filename == '': filename = default_filename
            else: filename = self.filename
        base, ext = splitext(filename)
        return base

    def write_bna(self, filename = ''):
        """Writes horizontal grid to Atlas BNA file."""
        filename = self.filename_base(filename) + '.bna'
        f = open(filename, 'w')
        headerfmt = '"%3s","",%1d\n'
        nodefmt = '%10.2f,%10.2f\n'
        for col in self.columnlist:
            f.write(headerfmt %  (col.name, col.num_nodes + 1))
            for node in col.node + [col.node[0]]:
                f.write(nodefmt % (node.pos[0], node.pos[1]))
        f.close()

    def write_layer_bln(self, filename = '', aspect = 8.0, left = 0.0):
        """Writes layer grid to Golden Software blanking (BLN) file."""
        filename = self.filename_base(filename) + '_layers.bln'
        f = open(filename, 'w')
        width = (self.layerlist[1].top - self.layerlist[-1].bottom) / aspect
        right = left + width
        nodefmt = '%10.2f,%10.2f\n'
        for layer in self.layerlist:
            f.write('5,1\n')
            pts = [
                (left, layer.top),
                (left, layer.bottom),
                (right, layer.bottom),
                (right, layer.top),
                (left, layer.top)]
            for x, z in pts: f.write(nodefmt % (x, z))
        f.close()

    def write_bna_labels(self, filename = ''):
        """Writes label file for BNA file (containing the column names)."""
        filename = self.filename_base(filename) + '_column_names.csv'
        f = open(filename, 'w')
        fmt = '%.2f,%.2f,"%s"\n'
        for col in self.columnlist:
            f.write(fmt % tuple(list(col.centre) + [col.name]))
        f.close()

    def write_layer_bln_labels(self, filename = '', aspect = 8.0, left = 0.0):
        """Writes label files for layer BLN file (containing the bottom
        elevations, centres and layer names)."""
        base = self.filename_base(filename)
        width = (self.layerlist[1].top - self.layerlist[-1].bottom) / aspect
        right = left + width
        centre = left + 0.5 * width
        labels = ['bottom_elevation', 'centre', 'name']
        filenames = [base + '_layer_' + label + 's.csv' for label in labels]
        files = dict(zip(labels, [open(filename, 'w') for filename in filenames]))
        fmt = dict(zip(labels, ['%.2f,%.2f,%.2f\n'] * 2 + ['%.2f,%.2f,"%s"\n']))
        start = dict(zip(labels, [0, 1, 1]))
        for label in labels: files[label].write('"X", "Y", "' + label + '"\n')
        for i, layer in enumerate(self.layerlist):
            for label in labels:
                vals = dict(zip(labels, [(left, layer.bottom, layer.bottom),
                                         (right, layer.centre, layer.centre),
                                         (centre, layer.centre, layer.name)]))
                if i >= start[label]:
                    files[label].write(fmt[label]%vals[label])
        for f in files.values(): f.close()

    def export_surfer(self, filename = '', aspect = 8.0, left = 0.0):
        """Writes files used for plotting geometry in Surfer."""
        self.write_bna(filename)
        self.write_bna_labels(filename)
        self.write_layer_bln(filename, aspect, left)
        self.write_layer_bln_labels(filename, aspect, left)

    def write_vtk(self, filename = '', arrays = None, wells = False, blockmap = {},
                  surface_snap = 0.1):
        """Writes *.vtu file for a vtkUnstructuredGrid object corresponding to the grid in 3D,
        with the specified filename, for visualisation with VTK."""
        from vtk import vtkXMLUnstructuredGridWriter
        base = self.filename_base(filename)
        filename = base + '.vtu'
        if wells: self.write_well_vtk(filename)
        if arrays is None: arrays = self.get_vtk_data(blockmap)
        vtu = self.get_vtk_grid(arrays, surface_snap)
        writer = vtkXMLUnstructuredGridWriter()
        writer.SetFileName(filename)
        if hasattr(writer, 'SetInput'): writer.SetInput(vtu)
        elif hasattr(writer, 'SetInputData'): writer.SetInputData(vtu)
        writer.Write()

    def get_well_vtk_grid(self):
        """Returns a VTK grid corresponding to the wells in the geometry."""
        from vtk import vtkUnstructuredGrid, vtkPoints, vtkIdList, vtkCharArray
        grid = vtkUnstructuredGrid()
        num_deviations = sum([well.num_deviations + 1 for well in self.welllist])
        pts = vtkPoints()
        pts.SetNumberOfPoints(num_deviations)
        i = 0
        for well in self.welllist:
            for p in well.pos:
                pts.SetPoint(i, list(p))
                i += 1
        grid.SetPoints(pts)
        VTK_POLY_LINE = 4
        i = 0
        for well in self.welllist:
            ids = vtkIdList()
            for p in well.pos:
                ids.InsertNextId(i)
                i += 1
            grid.InsertNextCell(VTK_POLY_LINE, ids)
        namearray = vtkCharArray()
        # SetTupleValue() was changed to SetTypedTuple() in VTK 7.1:
        if hasattr(vtkCharArray, 'SetTupleValue'):
            namearray.SetTypedTuple = namearray.SetTupleValue
        string_length = 5
        namearray.SetName('Name')
        namearray.SetNumberOfComponents(string_length)
        namearray.SetNumberOfTuples(self.num_wells)
        for i, well in enumerate(self.welllist):
            namearray.SetTypedTuple(i, well.name)
        grid.GetCellData().AddArray(namearray)
        return grid

    def write_well_vtk(self, filename = ''):
        from vtk import vtkXMLUnstructuredGridWriter
        filename = self.filename_base(filename) + '_wells.vtu'
        vtu = self.get_well_vtk_grid()
        writer = vtkXMLUnstructuredGridWriter()
        writer.SetFileName(filename)
        if hasattr(writer, 'SetInput'): writer.SetInput(vtu)
        elif hasattr(writer, 'SetInputData'): writer.SetInputData(vtu)
        writer.Write()

    def write_exodusii(self, filename = '', arrays = None,
                       blockmap = {}, surface_snap = 0.1):
        """Writes ExodusII file for a vtkUnstructuredGrid object corresponding
        to the grid in 3D, with the specified filename.
        """
        try:
            from vtk import vtkExodusIIWriter
        except ImportError:
            raise Exception("Either you don't have the Python VTK library installed," + \
                            " or it is too old.\n" + \
                            "On Windows you will need to have version 6.1 or later.")
        base = self.filename_base(filename)
        filename = base + '.exo'
        if arrays is None:
            arrays = self.get_vtk_data(blockmap)
        vtu = self.get_vtk_grid(arrays, surface_snap)
        writer = vtkExodusIIWriter()
        writer.SetFileName(filename)
        if hasattr(writer, 'SetInput'): writer.SetInput(vtu)
        elif hasattr(writer, 'SetInputData'): writer.SetInputData(vtu)
        writer.Write()

    def write_mesh(self, filename, surface_snap = 0.1, dimension = 3,
                   slice = None, file_format = None):
        """Writes mesh file for the grid, with the specified filename. The
        exported mesh file type is determined from the file extension
        of the filename. If dimension = 3, the full 3D mesh is
        written; if dimension = 2, the 2D horizontal mesh is written,
        unless a slice is specified, in which case a 2D vertical slice
        mesh is written.
        """
        try: import meshio
        except ImportError:
            raise Exception("Can't find meshio library- this function " + \
                            "requires it to be installed.")
        points, cells = self.meshio_grid(surface_snap, dimension, slice)
        meshio.write(filename, points, cells, file_format = file_format)

    def snap_columns_to_layers(self, min_thickness = 1.0, columns = []):
        """Snaps column surfaces to the bottom of their layers, if the surface
        block thickness is smaller than a given value.  This can be
        carried out over an optional subset of columns in the grid,
        otherwise over all columns.
        """
        if min_thickness > 0.0:
            if columns == []: columns = self.columnlist
            else:
                if isinstance(columns[0], str):
                    columns = [self.column[col] for col in columns]
            for col in columns:
                toplayer = self.column_surface_layer(col)
                if col.surface - toplayer.bottom < min_thickness:
                    col.surface = toplayer.bottom
                    col.num_layers -= 1
            self.setup_block_name_index()
            self.setup_block_connection_name_index()

    def subdivide_column(self, column_name, i0, colnodelist,
                         chars = ascii_lowercase, spaces = True):
        """Replaces specified column with columns based on subdividing it
        according to the specified local starting node index and list
        of column definitions (each a tuple of either local node
        indices, or 'c' denoting a new node at the centre of the
        column. Returns a list of the names of the new columns.
        """
        justfn = [str.ljust, str.rjust][self.right_justified_names]
        col = self.column[column_name]
        if any(['c' in newcol for newcol in colnodelist]):
            newnodename, nodenumber = self.new_node_name(justfn = justfn,
                                                         chars = chars,
                                                         spaces = spaces)
            centrenode = node(newnodename, col.centre)
            self.add_node(centrenode)
        colnumber = 0
        newcolnames = []
        for colnodes in colnodelist:
            nodes = [centrenode if i == 'c' else
                     col.node[col.index_plus(i0, i)] for
                     i in colnodes]
            name, colnumber = self.new_column_name(colnumber, justfn, chars)
            self.add_column(column(name, nodes, surface = col.surface))
            self.columnlist[-1].num_layers = col.num_layers
            newcolnames.append(name)
        self.delete_column(column_name)
        return newcolnames

    def triangulate_column(self, column_name, replace = True,
                           chars = ascii_lowercase, spaces = True):
        """Replaces specified column with triangulated columns based on a new
        node at its centre, and returns list of new columns created.
        """
        colnodelist = []
        col = self.column[column_name]
        for i, node in enumerate(col.node):
            inext = col.index_plus(i, 1)
            colnodelist.append((i, inext, 'c'))
        colnames = self.subdivide_column(column_name, 0, colnodelist, chars, spaces)
        return colnames

    def decompose_column(self, column_name, chars = ascii_lowercase,
                         spaces = True):
        """Replaces specified column with triangular or quadrilateral columns
        covering the original one, and returns a list of the new
        columns.  There are special cases for columns with lower
        numbers of sides, and 'straight' nodes (i.e. nodes that are on
        a straight line between their adjacent nodes).  Returns a list
        of names of columns that decompose the original column.
        """
        col = self.column[column_name]
        nn = col.num_nodes
        if nn <= 4: return [column_name]
        elif nn <= 8:
            angles = col.interior_angles
            tol = 1.e-3
            straight = [i for i, angle in enumerate(angles) if angle > np.pi - tol]
            ns = len(straight)
            if (nn, ns) == (5, 1):
                return self.subdivide_column(column_name, straight[0],
                                             [(0, 1, 2), (0, 2, 3), (0, 3, 4)],
                                             chars, spaces)
            elif (nn, ns) == (6, 2):
                d = col.index_dist(straight[0], straight[1])
                if d == 2:
                    last2 = [col.index_minus(i, 2) for i in straight]
                    start = [s for s, l in zip(straight, last2) if l not in straight][0]
                    return self.subdivide_column(column_name, start,
                                                 [(0, 1, 2, 'c'), (2, 3, 'c'),
                                                  (3, 4, 'c'), (4, 5, 'c'),
                                                  (5, 0, 'c')], chars, spaces)
                elif d == 3:
                    return self.subdivide_column(column_name, straight[0],
                                                 [(0, 1, 2, 3), (3, 4, 5, 0)],
                                                 chars, spaces)
                else: return self.triangulate_column(column_name, chars, spaces)
            elif (nn, ns) == (7, 3):
                last2 = [col.index_minus(i, 2) for i in straight]
                start = [s for s, l in zip(straight, last2) if l not in straight][0]
                return self.subdivide_column(column_name, start,
                                             [(0, 1, 2), (2, 3, 4),
                                              (0, 2, 4), (4, 5, 6, 0)],
                                             chars, spaces)
            elif (nn, ns) == (8, 4):
                return self.subdivide_column(column_name, straight[0],
                                             [(1, 2, 'c', 0), (2, 3, 4, 'c'),
                                              (4, 5, 6, 'c'), (6, 7, 0, 'c')],
                                             chars, spaces)
            else: return self.triangulate_column(column_name, chars, spaces)
        else: return self.triangulate_column(column_name, chars, spaces)

    def decompose_columns(self, columns = [], mapping = False,
                          chars = ascii_lowercase, spaces = True):
        """Decomposes columns with more than four sides to triangles and
        quadrilaterals. Optionally returns a dictionary mapping column
        names in the original geometry to lists of corresponding
        column names in the reduced geometry.
        """
        if columns == []: columns = self.columnlist
        else:
            if isinstance(columns[0], str): columns = [self.column[col] for col in columns]
        colmap = dict([(col.name, self.decompose_column(col.name, chars, spaces))
                       for col in columns])
        for c in self.missing_connections: self.add_connection(c)
        self.setup_block_name_index()
        self.setup_block_connection_name_index()
        if mapping: return colmap

    def fit_columns(self, data, alpha = 0.1, beta = 0.1, columns = [],
                    min_columns = [], grid_boundary = False,
                    silent = False, output_dict = False):
        """Fits scattered data to the column centres of the geometry, using
        least-squares bilinear finite element fitting with Sobolev
        smoothing.  The parameter data should be in the form of a
        3-column array with x,y,z data in each row.  The smoothing
        parameters alpha and beta control the first and second
        derivatives of the surface.  If the parameter columns is
        specified, data will only be fitted for the specified column
        names.  For columns with names in min_columns, column centre
        values will be calculated as the minimum of the fitted nodal
        values.  For all other columns, the average of the nodal
        values is used.  If grid_boundary is True, only data inside
        the bounding polygon of the grid are used- this can speed up
        the fitting if there are many data outside the grid, and the
        grid has a simply-shaped boundary.  The result is by default
        an array of fitted values corresponding to each of the
        columns. If output_dict is True, the result is a dictionary of
        fitted values indexed by column names.
        """

        if columns == []: columns = self.columnlist
        else: 
            if isinstance(columns[0], str):
                columns = [self.column[col] for col in columns]
        if min_columns != []:
            if not isinstance(min_columns[0], str):
                min_columns = [col.name for col in min_columns]

        # make copy of geometry and decompose into 3, 4 sided columns:
        geo = mulgrid(convention = self.convention)
        for n in self.nodelist: geo.add_node(node(n.name, n.pos))
        for col in self.columnlist:
            geo.add_column(column(col.name, [geo.node[n.name] for n in col.node]))
        colmap = geo.decompose_columns(columns, mapping = True,
                                       chars = ascii_lowercase + ascii_uppercase,
                                       spaces = True)
        geo_columns = []
        for col in columns: geo_columns += [geo.column[geocol] for
                                            geocol in colmap[col.name]]

        nodes = geo.nodes_in_columns(geo_columns)
        node_index = dict([(n.name, i) for i, n in enumerate(nodes)])
        num_nodes = len(nodes)
        # assemble least squares FEM fitting system:
        from scipy import sparse
        A = sparse.lil_matrix((num_nodes, num_nodes))
        b = np.zeros(num_nodes)
        guess = None
        if grid_boundary: bounds = geo.boundary_polygon
        else: bounds = None
        qtree = geo.column_quadtree(columns)
        import sys
        nd = len(data)
        for idata, d in enumerate(data):
            col = geo.column_containing_point(d[0:2], geo_columns,
                                              guess, bounds, qtree)
            percent = 100. * idata / nd
            if not silent:
                ps = 'fit_columns %3.0f%% done'% percent
                sys.stdout.write('%s\r' % ps)
                sys.stdout.flush()
            if col:
                xi = col.local_pos(d[0:2])
                if xi is not None:
                    guess = col
                    psi = col.basis(xi)
                    for i, nodei in enumerate(col.node):
                        I = node_index[nodei.name]
                        for j, nodej in enumerate(col.node):
                            J = node_index[nodej.name]
                            A[I, J] += psi[i] * psi[j]
                        b[I] += psi[i] * d[2]

        # add smoothing:
        smooth = {3: 0.5 * alpha * np.array([[1., 0., -1.], 
                                             [0., 1., -1.],
                                             [-1., -1., 2.]]),
                  4: alpha / 6. * np.array([[4., -1., -2., -1.],
                                            [-1., 4., -1., -2.],
                                            [-2., -1., 4., -1.],
                                            [-1., -2., -1., 4.]]) + \
                      beta * np.array([[1., -1., 1., -1.],
                                       [-1., 1., -1., 1.],
                                       [1., -1., 1., -1.],
                                       [-1., 1., -1., 1.]])}
        for col in geo_columns:
            for i, nodei in enumerate(col.node):
                I = node_index[nodei.name]
                for j, nodej in enumerate(col.node):
                    J = node_index[nodej.name]
                    A[I, J] += smooth[col.num_nodes][i, j]

        A = A.tocsr()
        from scipy.sparse.linalg import spsolve, use_solver
        use_solver(useUmfpack = False)
        z = spsolve(A, b)

        column_values = []
        def colnodez(col): return [z[node_index[node.name]] for node in col.node]
        for col in columns:
            if col.name in min_columns:
                geocol_min = None
                for geocolname in colmap[col.name]:
                    geocol = geo.column[geocolname]
                    nodez = colnodez(geocol)
                    minnodez = min(nodez)
                    if geocol_min is None: geocol_min = minnodez
                    else: geocol_min = min(geocol_min, minnodez)
                column_values.append(geocol_min)
            else:
                geocol_area, geocol_values = [], []
                for geocolname in colmap[col.name]:
                    geocol = geo.column[geocolname]
                    nodez = colnodez(geocol)
                    geocol_area.append(geocol.area)
                    geocol_values.append(sum(nodez) / geocol.num_nodes)
                val = sum([area * val for area, val in
                           zip(geocol_area, geocol_values)]) / sum(geocol_area)
                column_values.append(val)
        if output_dict:
            return dict(zip([col.name for col in columns], column_values))
        else: return np.array(column_values)

    def fit_surface(self, data, alpha = 0.1, beta = 0.1, columns = [],
                    min_columns = [], grid_boundary = False,
                    layer_snap = 0.0, silent = False):
        """Fits column surface elevations to the grid from the data, using the
        fit_columns() method (see documentation for that method for
        more detail). The layer_snap parameter can be specified as a
        positive number to avoid the creation of very thin top surface
        layers, if the fitted elevation is very close to the bottom of
        a layer.  In this case the value of layer_snap is a tolerance
        representing the smallest permissible layer thickness.
        """

        if columns == []: columns = self.columnlist
        else:
            if isinstance(columns[0], str):
                columns = [self.column[col] for col in columns]

        col_elevations = self.fit_columns(data, alpha, beta, columns,
                                          min_columns, grid_boundary, silent)

        for col, elev in zip(columns, col_elevations):
            col.surface = elev
            self.set_column_num_layers(col)

        self.snap_columns_to_layers(layer_snap, columns)
        self.setup_block_name_index()
        self.setup_block_connection_name_index()

    def refine(self, columns = [], bisect = False,
               bisect_edge_columns = [], chars = ascii_lowercase,
               spaces = True):
        """Refines selected columns in the grid.  If no columns are specified,
        all columns are refined.  Refinement is carried out by
        splitting: each column is divided into four, unless the bisect
        parameter is 'x' or 'y', in which case they are divided in the
        specified direction, or unless bisect is True, in which case
        they are divided into two between their longest sides.
        Triangular transition columns are added around the edge of the
        refinement region as needed.  Only 3 and 4-sided columns are
        supported.  The parameter bisect_edge_columns can contain a
        list of columns outside the edge of the refinement area (as
        specified by the columns parameter) which should be bisected
        prior to the refinement.  This is useful for columns with
        larger aspect ratios just outside the refinement area, whose
        aspect ratios would become even greater from simple
        refinement.
        """
        if columns == []: columns = self.columnlist
        else:
            if isinstance(columns[0], str):
                columns = [self.column[col] for col in columns]
        connections = set([])
        sidenodes = {}
        chars = uniqstring(chars)
        nodenumber, colnumber  = 0, 0
        justfn = [str.ljust, str.rjust][self.right_justified_names]
        def create_mid_node(node1, node2, sidenodes, nodenumber):
            midpos = 0.5 * (node1.pos + node2.pos)
            nodenames = frozenset((node1.name, node2.name))
            name, nodenumber = self.new_node_name(nodenumber, justfn, chars, spaces)
            self.add_node(node(name, midpos))
            sidenodes[nodenames] = self.nodelist[-1]
            return sidenodes, nodenumber
        if bisect:
            if bisect == True: direction = None
            else: direction = bisect
            for col in columns:
                for i in col.bisection_sides(direction):
                    n1, n2 = col.node[i], col.node[(i + 1) % col.num_nodes]
                    con = self.connection_with_nodes([n1, n2])
                    if con: connections.add(con)
                    else:
                        sidenodes, nodenumber = create_mid_node(n1, n2,
                                                                sidenodes, nodenumber)
        else:
            for col in columns: connections = connections | col.connection
        if bisect_edge_columns != []:
            if isinstance(bisect_edge_columns[0], str):
                bisect_edge_columns = [self.column[col] for col in bisect_edge_columns]
        columns_plus_edge = set(columns) | set(bisect_edge_columns)
        for con in connections: columns_plus_edge = columns_plus_edge | set(con.column)
        if all([col.num_nodes in [3, 4] for col in columns_plus_edge]):
            # bisect edge columns if required:
            for col in bisect_edge_columns:
                for con in col.connection:
                    if all([concol in bisect_edge_columns for concol in con.column]):
                        connections.add(con)
            # create midside nodes at connections:
            for con in connections:
                sidenodes, nodenumber = create_mid_node(con.node[0], con.node[1],
                                                        sidenodes, nodenumber)
            if not bisect:
                # create midside nodes on grid boundaries in the refinement area:
                bdy = self.boundary_nodes
                for col in columns:
                    nn = col.num_nodes
                    for i, corner in enumerate(col.node):
                        next_corner = col.node[(i + 1) % nn]
                        if (corner in bdy) and (next_corner in bdy):
                            sidenodes, nodenumber = create_mid_node(corner, next_corner,
                                                                    sidenodes, nodenumber)
            def transition_type(nn, sides):
                # returns transition type- classified by how many
                # refined sides, starting side, and range
                nref = len(sides)
                missing = list(set(range(nn)) - set(sides))
                nunref = len(missing)
                if nref == 1: return 1, sides[0], 0
                elif nref == nn: return nn, 0, nn - 1
                elif nunref == 1: return nref, (missing[0] + 1) % nn, nn - 2
                elif nn == 4 and nref == 2:
                    diff=sides[1] - sides[0]
                    if diff < 3: return nref, sides[0], diff
                    else: return nref, sides[1], 1
                else:
                    print('Error in refine()- unrecognised transition type:',
                          nref, 'of', nn, 'sides refined.')
            # How to subdivide, based on no. of nodes, no. of
            # refined sides and range of refined sides
            transition_column = {3: {(1, 0): ((0, (0, 1), 2), ((0, 1), 1, 2)),  
                                     (2, 1): ((0, (0, 1), (1, 2), 2), ((0, 1), 1, (1, 2))), 
                                     (3, 2): ((0, (0, 1), (2, 0)), ((0, 1), 1, (1, 2)),
                                              ((1, 2), 2, (2, 0)), ((0, 1), (1, 2), (2, 0)))},
                                 4: {(1, 0): ((0, (0, 1), 3), ((0, 1), 1, 2),
                                              ((0, 1), 2, 3)),
                                     (2, 1): ((0, (0, 1), 'c'), ((0, 1), 1, 'c'),
                                              (1, (1, 2), 'c'), ((1, 2), 2, 'c'),
                                              (2, 3, 'c'), (0, 'c', 3)),
                                     (2, 2): ((0, (0, 1), (2, 3), 3),
                                              ((0, 1), 1, 2, (2, 3))),
                                     (3, 2): ((0, (0, 1), (2, 3), 3),
                                              ((0, 1), 1, (1, 2)), ((1, 2), 2, (2, 3)),
                                              ((0, 1), (1, 2), (2, 3))),
                                     (4, 3): ((0, (0, 1), 'c', (3, 0)),
                                              ((0, 1), 1, (1, 2), 'c'),
                                              ((1, 2), 2, (2, 3), 'c'),
                                              ((2, 3), 3, (3, 0), 'c'))}}
            # create refined columns (and centre nodes for quadrilaterals that need them):
            centrenodes = {}
            for col in columns_plus_edge:
                nn = col.num_nodes
                refined_sides = []
                for i, corner in enumerate(col.node):
                    if frozenset((corner.name, col.node[(i + 1) % nn].name)) in sidenodes:
                        refined_sides.append(i)
                nrefined, istart, irange = transition_type(nn, refined_sides)
                if (col.num_nodes == 4) and ((nrefined == 4) or
                                             ((nrefined == 2) and (irange == 1))):
                    # create quadrilateral centre node:
                    name, nodenumber = self.new_node_name(nodenumber, justfn, chars, spaces)
                    self.add_node(node(name, col.centre))
                    centrenodes[col.name] = self.nodelist[-1]
                for subcol in transition_column[nn][nrefined, irange]:
                    name, colnumber = self.new_column_name(colnumber, justfn, chars)
                    nodes = []
                    for vert in subcol:
                        if isinstance(vert, int): n = col.node[(istart + vert) % nn]
                        elif vert == 'c': n = centrenodes[col.name]
                        else:
                            n = sidenodes[frozenset([col.node[(istart + i) % nn].name for
                                                     i in vert])]
                        nodes.append(n)
                    self.add_column(column(name, nodes, surface = col.surface))
                    self.columnlist[-1].num_layers = col.num_layers
            # clean up:
            for col in columns_plus_edge: self.delete_column(col.name)
            for con in self.missing_connections: self.add_connection(con)
            self.identify_neighbours()
            self.setup_block_name_index()
            self.setup_block_connection_name_index()
        else: print('Grid selection contains columns with more than 4 nodes: not supported.')

    def refine_layers(self, layers = [], factor = 2, chars = ascii_lowercase,
                      spaces = True):
        """Refines selected layers in the grid.  If no layers are specified,
        all layers are refined.  Each layer is refined by the
        specified factor.  Layer names for all subsurface layers in
        the grid are regenerated in sequence.
        """
        chars = uniqstring(chars)
        if layers == []: layers = self.layerlist
        else:
            if isinstance(layers[0], str):
                layers = [self.layer[lay] for lay in layers]
        factor = int(factor)
        top_elevation = self.layerlist[0].top
        atm_name = self.layerlist[0].name
        thicknesses = []
        for lay in self.layerlist[1:]:
            if lay in layers: thicknesses += [lay.thickness / factor] * factor
            else: thicknesses.append(lay.thickness)
        self.clear_layers()
        justify = ['l', 'r'][self.right_justified_names]
        self.add_layers(thicknesses, top_elevation, justify, chars, spaces)
        # Preserve old atmosphere layer name:
        self.rename_layer(self.layerlist[0].name, atm_name)
        for col in self.columnlist: self.set_column_num_layers(col)
        self.setup_block_name_index()
        self.setup_block_connection_name_index()

    def column_neighbour_groups(self, columns):
        """Given a list or set of columns, finds sets of columns that are
        connected together, and returns a list of them.
        """
        columns = list(set(columns))
        groups = []
        for col in columns: groups.append(set([col]))
        from copy import copy
        done = False
        while not done:
            done = True
            for i, g in enumerate(groups):
                ng = copy(g)
                for col in g: ng = ng | col.neighbour
                if i < len(groups) - 1:
                    for g2 in groups[i + 1:]:
                        if ng & g2:
                            g.update(g2)
                            groups.remove(g2)
                            done = False
                            break
                    if not done: break
        return groups

    def reduce(self, columns):
        """Reduce the geometry so that it contains only the specified columns."""
        if len(columns) > 0:
            if isinstance(columns[0], str):
                columns = [self.column[col] for col in columns]
        delete_columns = set(self.columnlist) - set(columns)
        colnames = [col.name for col in delete_columns]
        for colname in colnames: self.delete_column(colname)
        self.check(fix = True, silent = True)
        self.setup_block_name_index()
        self.setup_block_connection_name_index()

    def minc_array(self, vals, minc_indices, level = 0, outside = 0.0):
        """Returns an array for all blocks, with values taken from the vals array
        for the given MINC level, based on the index array minc_indices (which
        can be obtained from the output of the t2grid minc() method). For partial MINC
        grids, blocks outside the MINC area are assigned the value given by the
        parameter outside, unless outside is True, in which case the porous medium
        values are assigned."""

        if outside is True: vals_level = vals[: self.num_blocks]
        else:
            if outside is False: outside = 0.0
            vals_level = np.ones(self.num_blocks) * outside
        vals_level[minc_indices[0]] = vals[minc_indices[level]]
        return vals_level

    def from_amesh(self, input_filename = 'in', segment_filename = 'segmt',
              convention = 0, node_tolerance = None,
              justify = 'r', chars = ascii_lowercase, spaces = True):
        """Reads in AMESH input and segment files for a Voronoi mesh and
        returns a corresponding mulgrid object and block mapping. The
        block mapping dictionary maps block names in the geometry to
        block names in the AMESH grid. The atmosphere type is assumed
        to be 2 (no atmosphere blocks)."""

        from scipy.spatial import cKDTree
        thickness_tol = 1.e-3

        def parse_layers(filename):
            """Parse AMESH input to identify layer structure."""
            f = open(filename, 'r')
            found_locat = False
            for line in f:
                found_locat = line[:5].lower() == 'locat'
                if found_locat: break
            if not found_locat:
                raise Exception('Could not find locat block in AMESH input file: ' +
                                input_filename)
            layers = {}
            for line in f:
                if line.strip():
                    blkname = line[:5]
                    vals = line[5:].split()
                    index = int(vals[0])
                    x, y, z = float(vals[1]), float(vals[2]), float(vals[3])
                    pos = np.array([x, y])
                    thickness = float(vals[4])
                    if index in layers:
                        layers[index]['block_name'].append(blkname)
                        layers[index]['column_centre'][blkname] = pos
                        thickness_diff = thickness - layers[index]['thickness']
                        thickness_err = thickness_diff / layers[index]['thickness']
                        if abs(thickness_err) > thickness_tol:
                            raise Exception('Non-constant thickness ' +
                                            'for layer containing block: ' + blkname)
                    else:
                        layers[index] = {'block_name': [blkname],
                                         'column_centre': {blkname: pos},
                                         'thickness': thickness,
                                         'elevation': z}
                else: break
            layer_names = list(layers.keys())
            elevations = np.array([layers[name]['elevation'] for name in layer_names])
            isort = np.argsort(elevations)[::-1]
            layers = [layers[layer_names[i]] for i in isort]
            return layers

        def parse_segments(filename, bottom_layer):
            """Parse AMESH segment file and return list of 2-D segments, together
            with the minimum segment length."""
            from sys import float_info
            segment_data = []
            min_segment_length = float_info.max
            for line in open(filename):
                x1, y1, x2, y2 = float(line[0: 15]), float(line[15: 30]), \
                                 float(line[30: 45]), float(line[45: 60])
                points = (np.array([x1, y1]), np.array([x2, y2]))
                idx = int(line[60: 63])
                blocknames = (line[65: 70], line[70: 75])
                if all([blkname in bottom_layer['block_name'] or blkname.startswith('*')
                        for blkname in blocknames]):
                    segment_data.append({'points': points,
                                         'index': idx, 'blocknames': blocknames})
                    min_segment_length = min(min_segment_length,
                                             np.linalg.norm(points[0] - points[1]))
            return segment_data, min_segment_length

        layers = parse_layers(input_filename)
        bottom_layer = layers[-1]
        segment_data, min_segment_length = parse_segments(segment_filename, bottom_layer)

        justfn = [str.rjust, str.ljust][justify == 'l']
        geo = mulgrid(convention = convention, atmos_type = 2)

        # Add nodes:
        nodeindex = 1
        segments = {}
        for blkname in bottom_layer['block_name']: segments[blkname] = []
        if node_tolerance is None: node_tolerance = 0.9 * min_segment_length
        for seg in segment_data:
            nodes = []
            for point in seg['points']:
                new = True
                if geo.num_nodes > 1:
                    kdt = cKDTree([n.pos for n in geo.nodelist])
                    r,i = kdt.query(point)
                    if r < node_tolerance:
                        nodes.append(geo.nodelist[i]) # existing node
                        new = False
                if new: # new node
                    name = geo.node_name_from_number(nodeindex, justfn, chars, spaces)
                    newnode = node(name, point)
                    geo.add_node(newnode)
                    nodes.append(newnode)
                    nodeindex += 1
            for iname, blockname in enumerate(seg['blocknames']):
                if not blockname.startswith('*'):
                    segnodes = nodes if iname == 0 else nodes[::-1]
                    if segnodes not in segments[blockname]:
                        segments[blockname].append(segnodes)

        # Add columns:
        colindex = 1
        for blockname in bottom_layer['block_name']:
            segs = segments[blockname]
            if segs:
                colnodes = segs[0]
                done = False
                while not done:
                    nextsegs = [seg for seg in segs if seg[0] == colnodes[-1]]
                    try:
                        nextseg = nextsegs[0]
                        if nextseg[-1] == colnodes[0]: done = True
                        else: colnodes.append(nextseg[-1])
                    except: raise Exception(
                            "Could not identify column nodes for block:" + blockname)
                colname = geo.column_name_from_number(colindex, justfn, chars, spaces)
                colindex += 1
                pos = bottom_layer['column_centre'][blockname]
                geo.add_column(column(colname, colnodes, pos))
            else: raise Exception(
                    "No line segments found for block: " + blockname)
        # Add layers:
        top_elevation = layers[0]['elevation'] + 0.5 * layers[0]['thickness']
        geo.add_layers([lay['thickness'] for lay in layers], top_elevation,
                       justify, chars, spaces)
        for geolayer, lay in zip(geo.layerlist[1:], layers):
            geolayer.centre = lay['elevation']

        geo.set_default_surface()
        geo.identify_neighbours()
        geo.check(fix = True, silent = True)
        geo.setup_block_name_index()
        geo.setup_block_connection_name_index()

        # compute block mapping:
        orig_block_names = []
        orig_centres = []
        for lay in layers:
            orig_block_names += lay['block_name']
            for blkname in lay['block_name']:
                pos = np.hstack((lay['column_centre'][blkname],
                                 np.array([lay['elevation']])))
                orig_centres.append(pos)
        kdt = cKDTree(orig_centres)
        blockmap = {}
        for blkname in geo.block_name_list:
            layname = geo.layer_name(blkname)
            colname = geo.column_name(blkname)
            pos = geo.block_centre(layname, colname)
            r, i = kdt.query(pos)
            blockmap[blkname] = orig_block_names[i]

        return geo, blockmap

