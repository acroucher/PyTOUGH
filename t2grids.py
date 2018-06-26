"""For manipulating TOUGH2 grids.

Copyright 2011 University of Auckland.

This file is part of PyTOUGH.

PyTOUGH is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

PyTOUGH is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with PyTOUGH.  If not, see <http://www.gnu.org/licenses/>."""

from __future__ import print_function

from mulgrids import *
from t2incons import *

class rocktype(object):
    """Rock type"""
    def __init__(self, name = "dfalt", nad = 0, density = 2600.0, porosity = 0.1,
                 permeability = None, conductivity = 1.5,
                 specific_heat = 900.0):
        self.name = name
        self.nad = nad
        self.density = density
        self.porosity = porosity
        if permeability is None: permeability = np.ones(3) * 1.e-15
        if isinstance(permeability, list): permeability = np.array(permeability)
        self.permeability = permeability
        self.conductivity = conductivity
        self.specific_heat = specific_heat
        self.compressibility = 0.0
        self.expansivity = 0.0
        self.dry_conductivity = 0.0
        self.tortuosity = 0.0
        self.relative_permeability = {}
        self.capillarity = {}
    def __repr__(self):
        return self.name

class t2block(object):
    """Grid block"""
    def __init__(self, name = '     ', volume = 1.0, blockrocktype = None,
                 centre = None, atmosphere = False, ahtx = None, pmx = None,
                 nseq = None, nadd = None):
        if blockrocktype is None: blockrocktype = rocktype()
        self.name = name
        self.volume = volume
        self.rocktype = blockrocktype
        if isinstance(centre, list): centre = np.array(centre)
        self.centre = centre
        self.atmosphere = atmosphere
        self.ahtx = ahtx
        self.pmx = pmx
        self.nseq, self.nadd = nseq, nadd
        self.connection_name = set([])
    def __repr__(self): return self.name
    def get_num_connections(self): return len(self.connection_name)
    num_connections = property(get_num_connections)
    def get_neighbour_names(self):
        """Returns a set of neighbouring block names- those connected to this one."""
        return set([[blkname for blkname in cn if blkname != self.name][0] for
                    cn in self.connection_name])
    neighbour_name = property(get_neighbour_names)

class t2connection(object):
    """Connection between two blocks"""
    def __init__(self, blocks = None, direction = 1,
                 distance = None, area = 1.0, dircos = 0.0, sigma = None,
                 nseq = None, nad1 = None, nad2 = None):
        if blocks is None: blocks = [t2block(), t2block()]
        if distance is None: distance = [0.0, 0.0]
        self.block = blocks
        self.direction = direction # permeability direction
        self.distance = distance
        self.area = area
        self.dircos = dircos # direction cosine
        self.sigma = sigma # radiant emittance factor (TOUGH2)
        self.nseq, self.nad1, self.nad2 = nseq, nad1, nad2
        self.centre = None
        self.midpoint = None
        self.normal = None
    def __repr__(self):
        return self.block[0].name + ':' + self.block[1].name

class t2grid(object):
    """TOUGH2 grid"""
    def __init__(self): self.empty()

    def get_num_rocktypes(self):
        return len(self.rocktypelist)
    num_rocktypes = property(get_num_rocktypes)

    def get_num_blocks(self):
        return len(self.blocklist)
    num_blocks = property(get_num_blocks)

    def get_num_connections(self):
        return len(self.connectionlist)
    num_connections = property(get_num_connections)

    def get_num_atmosphere_blocks(self):
        return len(self.atmosphere_blocks)
    num_atmosphere_blocks = property(get_num_atmosphere_blocks)

    def get_num_underground_blocks(self):
        return self.num_blocks - self.num_atmosphere_blocks
    num_underground_blocks = property(get_num_underground_blocks)

    def get_atmosphere_blocks(self):
        return [blk for blk in self.blocklist if blk.atmosphere]
    atmosphere_blocks = property(get_atmosphere_blocks)

    def get_block_centres_defined(self):
        if self.num_atmosphere_blocks == 1: istart = 1
        else: istart = 0
        return all([blk.centre is not None for blk in self.blocklist[istart:]])
    block_centres_defined = property(get_block_centres_defined)

    def get_connection_centres_defined(self):
        return all([con.centre is not None for con in self.connectionlist])
    connection_centres_defined = property(get_connection_centres_defined)

    def calculate_block_centres(self, geo, blockmap = {}):
        """Calculates block centres from geometry object."""
        def mapname(blk): return blockmap[blk] if blk in blockmap else blk
        if geo.atmosphere_type == 0:
            istart = 1
            atmblockname = mapname(geo.block_name_list[0])
            # Centre not well defined for single atmosphere block:
            self.block[atmblockname].centre = None
        else: istart = 0
        for blkname in geo.block_name_list[istart:]:
            layername = geo.layer_name(blkname)
            colname = geo.column_name(blkname)
            self.block[mapname(blkname)].centre = geo.block_centre(layername, colname)

    def calculate_connection_centres(self, geo, blockmap = {}):
        """Calculates centre (and normal vector) for each connection face.  Note that the 'centre'
        depends on which block the connection is approached from- in the case where the connection
        face is not orthogonal to the line between the block centres.  Hence there are two 'centres'.
        The mipoint is just midway between the connection face nodes."""
        layindex = dict([(lay.name, i) for i, lay in enumerate(geo.layerlist)])
        iblockmap = dict([(v, k) for k, v in blockmap.items()])
        def imapblock(blk): return iblockmap[blk] if blk in iblockmap else blk
        for con in self.connectionlist:
            con.centre = {}
            geoblks = [imapblock(blk.name) for blk in con.block]
            layernames = [geo.layer_name(blk) for blk in geoblks]
            if layernames[0] == layernames[1]: # horizontal connection
                lay = geo.layer[layernames[0]]
                colnames = tuple([geo.column_name(blk) for blk in geoblks])
                geocon = geo.connection[colnames]
                nodepos = [node.pos for node in geocon.node]
                dpos = nodepos[1] - nodepos[0]
                for colname, blk in zip(colnames, con.block):
                    col = geo.column[colname]
                    if blk.centre is not None: vcentre = blk.centre[2]
                    else: vcentre = geo.block_centre(lay, col)[2]
                    hcentre = line_projection(col.centre, nodepos)
                    con.centre[blk.name] = np.hstack((hcentre, np.array([vcentre])))
                con.normal = np.array([dpos[1], -dpos[0], 0.]) / np.linalg.norm(dpos)
                con.midpoint = np.hstack((0.5 * sum(nodepos),
                                          min([pos[2] for name, pos in con.centre.items()])))
            else: # vertical connection
                layindices = np.array([layindex[layname] for layname in layernames])
                ilower = np.argmax(layindices)
                colname = geo.column_name(geoblks[ilower])
                col = geo.column[colname]
                hcentre = col.centre
                lay = geo.layer[layernames[ilower]]
                vcentre = geo.block_surface(lay, col)
                sgn = [1., -1.][ilower]
                con.normal = np.array([0., 0., sgn])
                con.midpoint = np.hstack((hcentre, np.array([vcentre])))
                for blk in con.block: con.centre[blk.name] = con.midpoint

    def rocktype_frequency(self, rockname):
        """Returns how many times the rocktype with given name is used in the grid."""
        return [blk.rocktype.name for blk in self.blocklist].count(rockname)

    def get_rocktype_frequencies(self):
        """Returns a list of tuples of occurring frequencies of rock types in
        the grid and the names of rocktypes with that frequency,
        ordered by increasing frequency.
        """
        freq = [(rt.name, self.rocktype_frequency(rt.name)) for rt in self.rocktypelist]
        occurring_freqs = list(set([item[1] for item in freq]))
        occurring_freqs.sort()
        frocks = dict([(f, []) for f in occurring_freqs])
        for item in freq: frocks[item[1]].append(item[0])
        return [(f, frocks[f]) for f in occurring_freqs]
    rocktype_frequencies = property(get_rocktype_frequencies)

    def sort_rocktypes(self):
        """Sorts rocktype list in alphabetical order by name."""
        rocknames = [rt.name for rt in self.rocktypelist]
        rocknames.sort()
        self.rocktypelist = [self.rocktype[name] for name in rocknames]

    def rename_rocktype(self, rockname, newrockname):
        """Renames rocktype with specified name. If that rocktype does not exist,
        or the target rocktype name has already been used, an exception is raised."""
        if rockname in self.rocktype:
            if newrockname in self.rocktype:
                raise Exception("Target rocktype name " + newrockname + " already exists.")
            else:
                rock = self.rocktype[rockname]
                del self.rocktype[rockname]
                rock.name = newrockname
                self.rocktype[newrockname] = rock
        else: raise Exception("Rocktype " + rockname + " not found.")

    def __repr__(self):
        return str(self.num_rocktypes) + ' rock types; ' + str(self.num_blocks) + \
            ' blocks; ' + str(self.num_connections) + ' connections'

    def __add__(self, other):
        """Adds two grids together."""
        result = t2grid()
        for grid in [self, other]:
            for rt in grid.rocktypelist: result.add_rocktype(rt)
            for blk in grid.blocklist: result.add_block(blk)
            for con in grid.connectionlist: result.add_connection(con)
        return result

    def embed(self, subgrid, connection):
        """Returns a grid with a subgrid embedded inside one of its blocks.
        The connection specifies how the two grids are to be
        connected: the blocks to be connected and the connection
        distances, area etc. between them.  The first block should be
        the host block, the second the connecting block in the
        subgrid.
        """
        result = None
        subvol = sum([blk.volume for blk in subgrid.blocklist])
        hostblock = connection.block[0]
        if subvol<hostblock.volume:
            dupblks = set([blk.name for blk in self.blocklist]) & \
                      set([blk.name for blk in subgrid.blocklist])
            if len(dupblks) == 0:
                result = self + subgrid
                connection.block = [result.block[blk.name] for blk in connection.block]
                result.add_connection(connection)
                result.block[hostblock.name].volume -= subvol # remove subgrid volume from host block
            else: print('Grid embedding error: the following blocks are in both grids:', dupblks)
        else: print('Grid embedding error: the host block is not big enough to contain the subgrid.')
        return result

    def empty(self):
        """Empties a TOUGH2 grid"""
        self.rocktypelist = []
        self.blocklist = []
        self.connectionlist = []
        self.rocktype = {}
        self.block = {}
        self.connection = {}

    def add_rocktype(self, newrocktype = None):
        """Adds a rock type to the grid.  Any existing rocktype of the same name is replaced."""
        if newrocktype is None: newrocktype = rocktype()
        if newrocktype.name in self.rocktype:
            i = self.rocktypelist.index(self.rocktype[newrocktype.name])
            self.rocktypelist[i] = newrocktype
        else: self.rocktypelist.append(newrocktype)
        self.rocktype[newrocktype.name] = newrocktype

    def delete_rocktype(self, rocktypename):
        """Deletes a rock type from the grid"""
        if rocktypename in self.rocktype:
            rt = self.rocktype[rocktypename]
            del self.rocktype[rocktypename]
            self.rocktypelist.remove(rt)

    def clean_rocktypes(self):
        """Deletes any unused rock types from the grid"""
        unused_rocktypes = []
        for rt in self.rocktypelist:
            if self.rocktype_frequency(rt.name) == 0: unused_rocktypes.append(rt.name)
        for name in unused_rocktypes: self.delete_rocktype(name)

    def add_block(self, newblock = None):
        """Adds a block to the grid"""
        if newblock is None: newblock = t2block()
        if newblock.name in self.block:
            i = self.blocklist.index(self.block[newblock.name])
            self.blocklist[i] = newblock
        else: self.blocklist.append(newblock)
        self.block[newblock.name] = newblock

    def delete_block(self, blockname):
        """Deletes a block from the grid"""
        if blockname in self.block:
            blk = self.block[blockname]
            from copy import copy
            connames = copy(blk.connection_name)
            for conname in connames: self.delete_connection(conname)
            del self.block[blockname]
            if blk in self.blocklist: self.blocklist.remove(blk)

    def demote_block(self, blockname):
        """Shifts blocks with specified names to the end of the block list.  The blockname
        parameter can be a block name or a list of them."""
        if isinstance(blockname, str): blockname = [blockname]
        for name in blockname:
            i = self.block_index(name)
            self.blocklist.append(self.blocklist.pop(i))

    def add_connection(self, newconnection = None):
        """Adds a connection to the grid"""
        if newconnection is None: newconnection = t2connection()
        conname = tuple([blk.name for blk in newconnection.block])
        if conname in self.connection:
            i = self.connectionlist.index(self.connection[conname])
            self.connectionlist[i] = newconnection
        else: self.connectionlist.append(newconnection)
        self.connection[conname] = newconnection
        for block in newconnection.block: block.connection_name.add(conname)

    def delete_connection(self, connectionname):
        """Deletes a connection from the grid"""
        if connectionname in self.connection:
            con = self.connection[connectionname]
            for block in con.block: block.connection_name.remove(connectionname)
            del self.connection[connectionname]
            self.connectionlist.remove(con)

    def block_index(self, blockname):
        """Returns index of block with specified name in the block list of the grid"""
        if blockname in self.block:
            return self.blocklist.index(self.block[blockname])
        else: return None

    def connection_index(self, connectionnames):
        """Returns index of connection with specified pair of names in the
        connection list of the grid"""
        if connectionnames in self.connection:
            return self.connectionlist.index(self.connection[connectionnames])
        else: return None

    def fromgeo(self, geo, blockmap = {}):
        """Converts a MULgraph grid to a TOUGH2 grid. The blockmap parameter
        applies an optional mapping to the block names from the geometry.
        """
        self.empty()
        self.add_rocktype(rocktype()) # add default rock type
        self.add_blocks(geo, blockmap)
        self.add_connections(geo, blockmap)
        return self

    def add_blocks(self, geo, blockmap = {}):
        """Adds blocks to grid from MULgraph geometry file"""
        self.add_atmosphereblocks(geo, blockmap)
        self.add_underground_blocks(geo, blockmap)

    def add_atmosphereblocks(self, geo, blockmap = {}):
        """Adds atmosphere blocks from geometry"""
        atmosrocktype = self.rocktypelist[0]
        if geo.atmosphere_type == 0: # one atmosphere block
            atmblockname = geo.block_name(geo.layerlist[0].name,
                                          geo.atmosphere_column_name, blockmap)
            centre = None
            self.add_block(t2block(atmblockname, geo.atmosphere_volume,
                                   atmosrocktype, centre = centre, atmosphere = True))
        elif geo.atmosphere_type == 1: # one atmosphere block per column
            for col in geo.columnlist:
                atmblockname = geo.block_name(geo.layerlist[0].name, col.name, blockmap)
                centre = geo.block_centre(geo.layerlist[0], col)
                self.add_block(t2block(atmblockname, geo.atmosphere_volume,
                                       atmosrocktype, centre = centre, atmosphere = True))

    def add_underground_blocks(self, geo, blockmap = {}):
        """Add underground blocks from geometry"""
        for lay in geo.layerlist[1:]:
            for col in [col for col in geo.columnlist if col.surface > lay.bottom]:
                name = geo.block_name(lay.name, col.name, blockmap)
                centre = geo.block_centre(lay, col)
                self.add_block(t2block(name, geo.block_volume(lay, col),
                                       self.rocktypelist[0], centre = centre))

    def add_connections(self, geo, blockmap = {}):
        """Add connections from geometry"""
        tilt = geo.tilt_vector
        for lay in geo.layerlist[1:]:
            layercols = [col for col in geo.columnlist if col.surface > lay.bottom]
            self.add_vertical_layer_connections(geo, lay, layercols, tilt, blockmap)
            self.add_horizontal_layer_connections(geo, lay, layercols, tilt, blockmap)

    def add_vertical_layer_connections(self, geo, lay, layercols = [],
                                       tilt = None, blockmap = {}):
        """Add vertical connections in layer"""
        if tilt is None: tilt = np.array([0., 0., -1.])
        for col in layercols:
            thisblk = self.block[geo.block_name(lay.name, col.name, blockmap)]
            if (geo.layerlist.index(lay) == 1) or (col.surface <= lay.top):
                # connection to atmosphere:
                abovelayer = geo.layerlist[0]
                abovedist = geo.atmosphere_connection
                belowdist = col.surface - thisblk.centre[2]
                if geo.atmosphere_type == 0:
                    aboveblk = self.blocklist[0]
                elif geo.atmosphere_type == 1:
                    aboveblk = self.block[geo.block_name(abovelayer.name, col.name, blockmap)]
                else: # no atmosphere blocks
                    continue
            else:
                ilayer = geo.layerlist.index(lay)
                abovelayer = geo.layerlist[ilayer - 1]
                aboveblk = self.block[geo.block_name(abovelayer.name, col.name, blockmap)]
                abovedist = aboveblk.centre[2] - abovelayer.bottom
                belowdist = lay.top - lay.centre
            con = t2connection([thisblk, aboveblk], 3, [belowdist, abovedist], col.area, tilt[2])
            self.add_connection(con)

    def add_horizontal_layer_connections(self, geo, lay, layercols = [],
                                         tilt = None, blockmap = {}):
        """Add horizontal connections in layer"""
        if tilt is None: tilt = np.array([0., 0., -1.])
        from math import cos, sin, radians
        layercolset = set(layercols)
        anglerad = radians(geo.permeability_angle)
        c, s = cos(anglerad), sin(anglerad)
        rotation = np.array([[c, s], [-s, c]])
        for con in [con for con in geo.connectionlist if set(con.column).issubset(layercolset)]:
            conblocks = [self.block[geo.block_name(lay.name, concol.name, blockmap)]
                         for concol in con.column]
            [dist, area] = geo.connection_params(con, lay)
            d = conblocks[1].centre - conblocks[0].centre
            d2 = np.dot(rotation, d[0:2])
            direction = np.argmax(abs(d2)) + 1
            dircos = np.dot(d, tilt) / np.linalg.norm(d)
            self.add_connection(t2connection(conblocks, direction, dist, area, dircos))

    def copy_connection_directions(self, geo, grid):
        """Copies connection permeability directions from another grid.  It is
        assumed the two grids have the same column structure.  The geo
        argument is the geometry file corresponding to grid.
        """
        nlayercons = len(geo.connectionlist)
        noldcons = len(grid.connectionlist)
        # create dictionary of permeability directions within a layer
        # (from bottom layer, assumed complete):
        dirn = {}
        for i, con in enumerate(geo.connectionlist):
            dirn[(con.column[0].name,
                  con.column[1].name)] = grid.connectionlist[noldcons - nlayercons + i].direction
        # transfer permeability directions to horizontal connections:
        for con in self.connectionlist:
            if con.direction<3:
                colnames = tuple([geo.column_name(blk.name) for blk in con.block])
                con.direction = dirn[colnames]

    def get_unconnected_blocks(self):
        """Returns a set of blocks in the grid that are not connected to any other blocks."""
        return set([blk.name for blk in self.blocklist if len(blk.connection_name) == 0])
    unconnected_blocks = property(get_unconnected_blocks)

    def get_isolated_rocktype_blocks(self):
        """Returns a list of blocks with isolated rocktypes- that is, blocks
        with a rocktype different from that of all other blocks they are connected to.
        """
        # blocks with volume greater than this are considered boundary
        # condition blocks and not counted:
        bc_volume = 1.e20
        return set([blk.name for blk in self.blocklist if (blk.volume<bc_volume) and
                    not (blk.rocktype.name in [self.block[nbr].rocktype.name for
                                               nbr in blk.neighbour_name])])
    isolated_rocktype_blocks = property(get_isolated_rocktype_blocks)

    def check(self, fix = False, silent = False):
        """Checks a grid for errors, and optionally fixes them.  Errors checked for are:
        - blocks not connected to any other blocks
        - blocks with isolated rocktypes
        Returns True if no errors were found, and False otherwise.  If silent is True,
        there is no printout. Unconnected blocks are fixed by deleting them.
        Isolated rocktype blocks are fixed by assigning them the most popular rocktype
        of their neighbours."""
        ok = True
        ub = self.unconnected_blocks
        if len(ub)>0:
            ok = False
            if not silent: print('Unconnected blocks:', list(ub))
            if fix:
                for blk in ub: self.delete_block(blk)
                if not silent: print('Unconnected blocks fixed.')
        ib = self.isolated_rocktype_blocks
        if len(ib)>0:
            ok = False
            if not silent: print('Isolated rocktype blocks:', list(ib))
            if fix:
                for blk in ib:
                    nbr_rocktype = [self.block[nbr].rocktype.name for
                                    nbr in self.block[blk].neighbour_name]
                    pop_rocktype = max(set(nbr_rocktype), key = nbr_rocktype.count)
                    self.block[blk].rocktype = self.rocktype[pop_rocktype]
                if not silent: print('Isolated rocktype blocks fixed.')
        if ok and not silent: print('No problems found.')
        return ok

    def get_rocktype_indices(self, geo = None, blockmap = {}):
        """Returns an integer array containing the rocktype index for each block in the grid."""
        rocknames = [rt.name for rt in self.rocktypelist]
        rockdict = dict([(name, i) for i, name in enumerate(rocknames)])
        if geo is None:
            return np.array([rockdict[blk.rocktype.name] for blk in self.blocklist])
        else:
            return np.array([rockdict[self.block[
                            blockmap[blkname] if blkname in blockmap else blkname].rocktype.name]
                             for blkname in geo.block_name_list])
    rocktype_indices = property(get_rocktype_indices)

    def get_vtk_data(self, geo, blockmap = {}):
        """Returns dictionary of VTK data arrays from rock types.  The
        geometry object geo must be passed in."""
        from vtk import vtkIntArray, vtkFloatArray, vtkCharArray
        arrays = {'Block':{'Rock type index': vtkIntArray(), 'Porosity': vtkFloatArray(),
                         'Permeability': vtkFloatArray(), 'Name': vtkCharArray()}, 'Node': {}}
        # SetTupleValue() was changed to SetTypedTuple() in VTK 7.1:
        if hasattr(vtkCharArray, 'SetTupleValue'):
            arrays['Block']['Name'].SetTypedTuple = arrays['Block']['Name'].SetTupleValue
        vector_properties = ['Permeability']
        string_properties = ['Name']
        string_length = 5
        nele = geo.num_underground_blocks
        array_length = {'Block': nele, 'Node': 0}
        for array_type, array_dict in arrays.items():
            for name, array in array_dict.items():
                array.SetName(name)
                if name in vector_properties:
                    array.SetNumberOfComponents(3)
                    array.SetNumberOfTuples(array_length[array_type])
                elif name in string_properties:
                    array.SetNumberOfComponents(string_length)
                    array.SetNumberOfTuples(array_length[array_type])
                else:
                    array.SetNumberOfComponents(1)
                    array.SetNumberOfValues(array_length[array_type])
        natm = geo.num_atmosphere_blocks
        rindex = self.get_rocktype_indices(geo, blockmap)
        rockdict = dict(zip([blk.name for blk in self.blocklist], rindex))
        for i, blkname in enumerate(geo.block_name_list[natm:]):
            mapped_name = blockmap[blkname] if blkname in blockmap else blkname
            arrays['Block']['Name'].SetTypedTuple(i, mapped_name)
            ri = rockdict[mapped_name]
            arrays['Block']['Rock type index'].SetValue(i, ri)
            rt = self.rocktypelist[ri]
            arrays['Block']['Porosity'].SetValue(i, rt.porosity)
            k = rt.permeability
            arrays['Block']['Permeability'].SetTuple3(i, k[0], k[1], k[2])
        return arrays

    def write_vtk(self, geo, filename, wells = False, blockmap = {}, surface_snap = 0.1):
        """Writes *.vtu file for a vtkUnstructuredGrid object corresponding to
        the grid in 3D, with the specified filename, for visualisation
        with VTK.
        """
        from vtk import vtkXMLUnstructuredGridWriter
        if wells: geo.write_well_vtk()
        arrays = geo.get_vtk_data(blockmap)
        grid_arrays = self.get_vtk_data(geo, blockmap)
        for array_type, array_dict in arrays.items():
            array_dict.update(grid_arrays[array_type])
        vtu = geo.get_vtk_grid(arrays, surface_snap)
        writer = vtkXMLUnstructuredGridWriter()
        writer.SetFileName(filename)
        if hasattr(writer, 'SetInput'): writer.SetInput(vtu)
        elif hasattr(writer, 'SetInputData'): writer.SetInputData(vtu)
        writer.Write()

    def flux_matrix(self, geo, blockmap = {}):
        """Returns a sparse matrix which can be used to multiply a vector of
        connection table values for underground blocks, to give
        approximate average fluxes of those values at the block centres.
        """
        natm = geo.num_atmosphere_blocks
        nele = geo.num_underground_blocks
        conindex = dict([((c.block[0].name, c.block[1].name), i) for
                         i, c in enumerate(self.connectionlist)])
        from scipy import sparse
        A = sparse.lil_matrix((3 * nele, self.num_connections))
        if not self.block_centres_defined: self.calculate_block_centres(geo, blockmap)
        if not self.connection_centres_defined: self.calculate_connection_centres(geo, blockmap)
        def mname(blk): return blockmap[blk] if blk in blockmap else blk
        for iblk, geoblkname in enumerate(geo.block_name_list[natm:]):
            blkname = mname(geoblkname)
            blk = self.block[blkname]
            ncons = blk.num_connections
            if ncons > 0:
                M, icons = [], []
                for conname in blk.connection_name:
                    con = self.connection[conname]
                    anormal = con.normal * con.area
                    row = list(anormal) # fit constant flows
                    if ncons >= 6: # fit linear flows
                        row += list((con.centre[blk.name] - blk.centre) * anormal)
                    M.append(row)
                    icons.append(conindex[conname])
                Ablk = -np.linalg.pinv(np.array(M))
                ib3 = iblk * 3
                for ic in range(ncons):
                    for ip in range(3): A[ib3 + ip, icons[ic]] = Ablk[ip, ic]
        return A

    def radial(self, rblocks, zblocks, convention = 0, atmos_type = 2,
               origin = None, justify = 'r',
               case = None, dimension = 2, blockmap = {},
               chars = ascii_lowercase, spaces = True):
        """Returns a radial TOUGH2 grid with the specified radial and vertical
        block sizes.  The arguments are arrays of the block sizes in
        each dimension (r,z).  Naming convention, atmosphere type and
        grid origin can optionally be specified.  The origin is in
        (r,z) coordinates, so origin[0] is the starting radius of the
        grid.  (The origin can also be specified with three
        components, in which case the second one is ignored.)  The
        optional justify and case parameters specify the format of the
        character part of the block names (whether they are right or
        left justified, and lower or upper case). The case parameter
        is now deprecated- the more flexible chars parameter should be
        used instead.  Specifying dimension<>2 (between 1 and 3)
        simulates flow in fractured rock using the "generalized radial
        flow" concept of Barker, J.A. (1988), "A generalized radial
        flow model for hydraulic tests in fractured rock", Water
        Resources Research 24(10), 1796-1804.  In this case it
        probably doesn't make much sense to have more than one block
        in the z direction.
        """

        if origin is None: origin = np.array([0.,0.])
        if isinstance(rblocks, list): rblocks = np.array(rblocks)
        if isinstance(zblocks, list): zblocks = np.array(zblocks)
        if isinstance(origin, list): origin = np.array(origin)
        if len(origin) > 2: origin = origin[[0, 2]]

        justfn = [str.rjust, str.ljust][justify == 'l']
        if case is not None:
            casefn = [str.upper, str.lower][case == 'l']
            chars = casefn(chars)
        chars = uniqstring(chars)

        n2 = 0.5 * dimension
        if dimension != 2: # need gamma function
            try:
                from math import gamma # included for Python 2.7 or later
            except ImportError:
                from scipy.special import gamma
            gamman2 = gamma(n2)
        else: gamman2 = 1.0
        alpha = 2. / gamman2 * np.pi ** n2

        b = np.sum(zblocks) # total thickness
        b2n = b ** (2 - dimension)
        r = origin[0] + np.concatenate((np.zeros(1), np.cumsum(rblocks)))
        rin, rout=r[:-1], r[1:] # inner and outer radii
        rc = 0.5 * (rin + rout) # centre radius
        A = alpha * b2n / dimension * np.diff(r ** dimension) # "top area"
        c = alpha * b2n * rout ** (dimension - 1) # "outer circumference"
        ncols = len(rblocks)

        grid = t2grid()
        grid.add_rocktype(rocktype()) # add default rock type

        # dummy geometry for creating block names etc:
        geo = mulgrid(type = 'GENER', convention = convention, atmos_type = atmos_type)
        for ir, dr in enumerate(rblocks):
            colname = geo.column_name_from_number(ir + 1, justfn, chars, spaces)
            geo.add_column(column(colname, [], centre = np.array([rc[ir], 0.])))
        geo.add_layers(zblocks, origin[1], justify, chars, spaces)
        grid.add_atmosphereblocks(geo)

        for lay in geo.layerlist[1:]: # add blocks
            V = A * lay.thickness
            for col, rcentre, vol in zip(geo.columnlist, rc, V):
                name = geo.block_name(lay.name, col.name, blockmap)
                centre = np.array([rcentre, 0., lay.centre])
                grid.add_block(t2block(name, vol, grid.rocktypelist[0], centre = centre))

        for ilay, lay in enumerate(geo.layerlist[1:]):
            Ar = c * lay.thickness
            for icol, col in enumerate(geo.columnlist): # vertical connections
                top_area = A[icol]
                blkindex = ilay * ncols + icol + geo.num_atmosphere_blocks
                thisblk = grid.blocklist[blkindex]
                if ilay == 0: # atmosphere connections
                    abovedist = geo.atmosphere_connection
                    belowdist = 0.5 * geo.layerlist[1].thickness
                    if atmos_type == 0: aboveblk = grid.blocklist[0]
                    elif atmos_type == 1: aboveblk = grid.blocklist[icol]
                    else: continue
                else:
                    abovelayer = geo.layerlist[ilay]
                    aboveblk = grid.blocklist[blkindex - ncols]
                    abovedist = aboveblk.centre[2] - abovelayer.bottom
                    belowdist = lay.top - lay.centre
                con = t2connection([thisblk, aboveblk], 3, [belowdist, abovedist], top_area, -1.0)
                grid.add_connection(con)
            for icol, col in enumerate(geo.columnlist[:-1]): # radial connections
                nextcol = geo.columnlist[icol + 1]
                conblocks = [grid.block[geo.block_name(lay.name, acol.name, blockmap)] for
                             acol in [col, nextcol]]
                dist, area = [0.5 * rblocks[icol], 0.5 * rblocks[icol + 1]], Ar[icol]
                direction, dircos = 1, 0.0
                grid.add_connection(t2connection(conblocks, direction, dist, area, dircos))

        return grid

    def incons(self, values = (101.3e3, 20.)):
        """Creates a t2incon initial condtions object corresponding to the
        grid from the given values.  If initial conditions are given
        for one block only, these are applied to all blocks.
        """
        inc = t2incon()
        values = np.array(values)
        if len(np.shape(values)) == 1:
            from copy import copy
            for blk in self.blocklist: inc[blk.name] = copy(tuple(values))
        else:
            for blk, val in zip(self.blocklist, values): inc[blk.name] = tuple(val)
        return inc

    def neargroups(self, blocknames):
        """Given a list or set of block names, finds groups of 'near' blocks.
        Blocks are assigned the same group if they are neighbours, or
        share a neighbour.
        """
        blocknames = list(set(blocknames))
        groups = []
        for blk in blocknames: groups.append(set([blk]))
        from copy import copy
        done = False
        while not done:
            done = True
            for i, g in enumerate(groups):
                ng = copy(g)
                for blk in g: ng = ng | self.block[blk].neighbour_name
                if i<len(groups) - 1:
                    for g2 in groups[i + 1:]:
                        ng2 = copy(g2)
                        for blk in g2: ng2 = ng2 | self.block[blk].neighbour_name
                        if ng & ng2:
                            g.update(g2)
                            groups.remove(g2)
                            done = False
                            break
                    if not done: break
        return groups

    def rectgeo(self, origin_block = None, atmos_volume = 1.e25, remove_inactive = False,
                convention = 0, atmos_type = 2, justify = 'r',
                chars = ascii_lowercase, spaces = True,
                layer_snap = 0.1):
        """For a rectangular grid, returns a mulgrid object representing the
        geometry. The 'origin block' (the block on the bottom layer of
        the grid, at the origin of the permeability direction 1 & 2
        axes) can optionally be specified manually. If the
        remove_inactive parameter is True, the mulgrid surface
        property will be determined by the positions of inactive
        blocks in the grid. The atmos_volume parameter specifies the
        maximum block volume considered to be part of the geometrical
        grid. The layer_snap parameter can be used to eliminate blocks
        with very small volumes at the ground surface.  The method
        also returns a block mapping dictionary, mapping geometry
        block names into grid block names.
        """

        def blockelevs(grid, max_volume = None):
            """Returns array of block elevations."""
            def elev(blk, max_volume = None):
                if blk.centre is None: return np.nan
                elif max_volume is None: return blk.centre[2]
                elif 0. < blk.volume <= max_volume: return blk.centre[2]
                else: return np.nan
            return np.array([elev(blk, max_volume) for blk in grid.blocklist])

        def topmost_block(grid, max_volume = None):
            """"Returns block with highest centre in the grid. If max_volume is
            given, only blocks with volume greater than zero and less
            than max_volume are considered.
            """
            itop = np.nanargmax(blockelevs(grid, max_volume))
            return grid.blocklist[itop]

        def find_origin_block(grid):
            """Returns block in bottom layer of rectangular mesh at horizontal
            position corresponding to the origin of the coordinate
            system, defined by the permeability directions assigned to
            the connections. It's assumed that the origin is the first
            block in the bottom layer.
            """
            iorigin = np.nanargmin(blockelevs(grid))
            return grid.blocklist[iorigin]

        def con_name_index(con, blkname):
            """Index of block name in connection block names."""
            if con[0] == blkname: return 0
            elif con[1] == blkname: return 1
            else: return None

        def next_block_in_direction(blk, last, direction, grid, max_volume = None):
            """Finds next block in specified direction, and connection connecting them."""
            cons = [con for con in blk.connection_name if
                    grid.connection[con].direction == direction]
            if last: cons = [con for con in cons if last.name not in con]
            nextblk, nextcon, found = None, None, None
            for con in cons:
                i = con_name_index(con, blk.name)
                i2 = (i + 1) % 2
                nextblk = grid.block[con[i2]]
                nextcon = con
                if max_volume is None: found = True
                else:
                    if 0. < nextblk.volume < max_volume: found = True
                if found: break
            if found: return nextblk, nextcon
            else: return None, None

        def block_direction_track(grid, start_block, dirn, max_volume = None):
            """Returns list of blocks and block sizes found by following specified
            direction from the starting block. Specify max_volume as a float to set
            maximum allowed block volume- used for excluding boundary condition
            blocks."""
            blks, sizes = [], []
            con, last, last_con = None, None, None
            blk = start_block
            done = False
            while not done:
                if max_volume is not None: ok = 0. < blk.volume < max_volume
                else: ok = True
                if ok: blks.append(blk)
                last_con = con
                next_blk,con = next_block_in_direction(blk, last, dirn, grid, max_volume)
                if next_blk:
                    if ok:
                        i = con_name_index(con, blk.name)
                        sizes.append(2. * grid.connection[con].distance[i])
                    last = blk
                    blk = next_blk
                else:
                    done = True
                    if last_con:
                        i = con_name_index(last_con, blk.name)
                        sizes.append(2. * grid.connection[last_con].distance[i])
            return blks, sizes

        def block_spacings(grid, ob, max_volume):
            """Calculate block spacings in each direction"""
            spacings = {}
            for dirn in range(1, 4):
                if dirn < 3: blk = ob
                else: blk = topmost_block(grid, max_volume)
                blks, spacings[dirn] = block_direction_track(grid, blk, dirn, max_volume)
                if dirn == 3 and blks[0].centre[2] < blks[-1].centre[2]:
                    spacings[dirn].reverse()
            # calculate missing spacing for 2D meshes:
            missing_dirns = [i for i in range(1,4) if len(spacings[i]) == 0]
            num_missing = len(missing_dirns)
            if num_missing == 1:
                dirn = missing_dirns[0]
                present_dirns = [1,2,3]; present_dirns.remove(dirn)
                d = ob.volume
                for pd in present_dirns:
                    con = [con for con in ob.connection_name if
                           grid.connection[con].direction == pd][0]
                    i = con_name_index(con, ob.name)
                    d /= (2. * grid.connection[con].distance[i])
                spacings[dirn].append(d)
            elif num_missing == 2:
                raise Exception("Mesh appears to be 1-D: can't reconstruct geometry.")
            return spacings

        def find_surface(geo, grid, blockmap, remove_inactive, max_volume):
            """Determines surface of mulgrid from the positions of missing or inactive
            blocks in grid. The blockmap parameter maps mulgrid block names to grid
            block names."""
            remove_blocks = []
            inactive = False
            for blk in grid.blocklist:
                if remove_inactive and blk.volume <= 0.: inactive = True
                if inactive: remove_blocks.append(blk)
            bottom_layer = geo.layerlist[-1]
            for col in geo.columnlist:
                geoblkname = geo.block_name(bottom_layer.name, col.name)
                blkname = blockmap[geoblkname]
                bottom_block = grid.block[blkname]
                colblocks, layerthicks = block_direction_track(grid, bottom_block, 3, max_volume)
                if colblocks:
                    remove_col = [blk for blk in colblocks if blk in remove_blocks]
                    for blk in remove_col:
                        i = colblocks.index(blk)
                        del colblocks[i]; del layerthicks[i]
                    topblock = colblocks[-1]
                    zc = topblock.centre[2]
                    if topblock.volume > 0.:
                        block_height = colblocks[-1].volume / col.area
                        if layerthicks: lt = layerthicks[-1]
                        else: lt = block_height
                        if block_height <= lt: surface = zc + 0.5 * block_height
                        else: surface = zc - 0.5 * lt + block_height
                    else: # top block is inactive
                        surface = zc + 0.5 * layerthicks[-1]
                    col.surface = surface
                    geo.set_column_num_layers(col)
            geo.setup_block_name_index()
            geo.setup_block_connection_name_index()
            return geo

        def match_position(geo, grid, ob):
            """Rotate and translate geometry as needed."""
            blks, sp = block_direction_track(grid, ob, 1)
            angle =  0.5 * np.pi - vector_heading(blks[-1].centre[0:2] - ob.centre[:2])
            from math import degrees
            angle = degrees(angle)
            geo.rotate(-angle, np.zeros(2))
            geo.permeability_angle = angle
            col = geo.columnlist[0]
            lay = geo.layerlist[-1]
            origin = geo.block_centre(lay, col)
            trans = ob.centre - origin
            geo.translate(trans)
            return geo

        def block_mapping(geo, grid, ob, nblks, max_volume):
            """Generates a mapping from geometry block names to grid block names."""
            mapping = {}
            icol = 0
            start2, last2 = ob, None
            for i2 in range(nblks[2]):
                start1, last1 = start2, None
                for i1 in range(nblks[1]):
                    col = geo.columnlist[icol]
                    blk, last3 = start1, None
                    for lay in geo.layerlist[::-1]:
                        geoblkname = geo.block_name(lay.name, col.name)
                        mapping[geoblkname] = blk.name
                        next_blk,con = next_block_in_direction(blk, last3, 3, grid, max_volume)
                        if next_blk is None: # incomplete column
                            atm_blk,con = next_block_in_direction(blk, last3, 3, grid)
                            if atm_blk:
                                if geo.atmosphere_type == 0:
                                    atmblockname = geo.block_name(geo.layerlist[0].name,
                                                                  geo.atmosphere_column_name)
                                    mapping[atmblockname] = atm_blk.name
                                elif geo.atmosphere_type == 1:
                                    atmblockname = geo.block_name(geo.layerlist[0].name, col.name)
                                    mapping[atmblockname] = atm_blk.name
                            break
                        last3 = blk; blk = next_blk
                    icol += 1
                    next1, con = next_block_in_direction(start1, last1, 1, grid, max_volume)
                    last1 = start1; start1 = next1
                next2, con = next_block_in_direction(start2, last2, 2, grid, max_volume)
                last2 = start2; start2 = next2
            return mapping

        def required_centres_present(grid, max_volume):
            """Returns True if all required block centres are specified."""
            return all([blk.centre is not None
                        for blk in grid.blocklist if 0. < blk.volume < max_volume])

        if not required_centres_present(self, atmos_volume):
            raise Exception("Not all required grid block centres are defined.")
        else:
            if isinstance(origin_block, str): ob = self.block[origin_block]
            elif isinstance(origin_block, t2block): ob = origin_block
            else: ob = find_origin_block(self)
            if ob is None: raise Exception("Can't find origin block for grid.")
            else:
                spacings = block_spacings(self, ob, atmos_volume)
                nblks = dict([(i, len(spacings[i])) for i in range(1,4)])
                geo = mulgrid().rectangular(spacings[1], spacings[2], spacings[3],
                                            convention = convention, atmos_type = atmos_type,
                                            justify = justify, chars = chars, spaces = spaces)
                blockmap = block_mapping(geo, self, ob, nblks, atmos_volume)
                geo = match_position(geo, self, ob)
                geo = find_surface(geo, self, blockmap, remove_inactive, atmos_volume)
                geo.snap_columns_to_layers(layer_snap)
                prune = list(set(blockmap.keys()) - set(geo.block_name_list))
                for blk in prune: del blockmap[blk]
                return geo, blockmap

    def minc(self, volume_fractions, spacing = 50., num_fracture_planes = 1,
             blocks = None, matrix_blockname = None, minc_rockname = None,
             proximity = None, atmos_volume = 1.e25, incon = None,
             fracture_connection_distance = 0.0):
        """Adds MINC blocks to the grid, and returns an array of block
        indices for the different MINC levels. The first three parameters define
        the fracture geometry. The blocks parameter is a list of blocks or block
        names specifying where MINC is to be applied. The matrix_blockname, 
        minc_rockname and proximity parameters are optional functions for 
        determining the names of matrix blocks for a given level and rocktype names
        for MINC blocks, and the proximity function. If these are not specified,
        defaults will be used. The atmos_volume parameter specifies the minimum
        volume of large blocks used for boundary conditions, which will be excluded
        from MINC processing. If incon is specifed, a new t2incon is also returned
        with initial conditions in the MINC blocks copied from the original."""

        num_levels = len(volume_fractions)
        if num_levels < 2:
            raise Exception("Need at least two volume fractions specified " +
                            "for MINC.")
        else:

            from scipy.optimize import bisect
            from scipy.misc import derivative
            from numbers import Number

            volume_fractions = np.array(volume_fractions, dtype = float64)

            if isinstance(spacing, Number): spacing = [spacing]
            missing = num_fracture_planes - len(spacing)
            if missing > 0: spacing = list(spacing) + [spacing[0]] * missing
            spacing = np.array(spacing)
            if blocks is None or blocks == []:
                blocks = [blk.name for blk in self.blocklist]
            if isinstance(blocks[0], t2block):
                blocks = [blk.name for blk in blocks]

            def default_proximity(x):
                """Default proximity functions."""
                if num_fracture_planes == 1:
                    if x >= 0.5 * spacing[0]: return 1.
                    else: return 2. * x / spacing[0]
                elif num_fracture_planes == 2:
                    if any([0.5 * fs <= x for fs in spacing[:2]]):
                        return 1.
                    else:
                        f01 = spacing[0] * spacing[1]
                        return 2. * x *(spacing[0] + spacing[1] - 2. * x) / f01
                elif num_fracture_planes == 3:
                    u = 2. * x / spacing[0:3]
                    if any([ui >= 1. for ui in u]): return 1.
                    else: return u[0] * u[1] * u[2] \
                            - (u[0] * u[1] + u[1] * u[2] + u[0] * u[2]) \
                            + u[0] + u[1] + u[2]
                else:
                    raise Exception("Invalid number of MINC fracture planes" +
                                    "(" + str(num_fracture_planes) + ").")

            if proximity is None: proximity = default_proximity

            def invert_proximity(vf, xl, xr):
                def proxv(x): return proximity(x) - vf
                while proxv(xr) < 0.: xr *= 2.
                x,r = bisect(proxv, xl, xr, full_output = True)
                if r.converged: return x, xr
                else: return None, None

            def inner_dist(x):
                """Returns innermost connection distance, at distance x."""
                if num_fracture_planes == 1:
                    return (spacing[0] - 2. * x) / 6.
                elif num_fracture_planes == 2:
                    u = spacing[0:2] - 2. * x
                    return 0.25 * u[0] * u[1] / (u[0] + u[1])
                elif num_fracture_planes == 3:
                    u = spacing[0:3] - 2. * x
                    return 0.3 * u[0] * u[1] * u[2] / (u[0] * u[1] + u[1] * u[2] + u[0] * u[2])
                else:
                    raise Exception("Invalid number of MINC fracture planes" +
                                    "(" + str(num_fracture_planes) + ").")

            # Calculate MINC geometry parameters:
            volume_fractions /= np.sum(volume_fractions)
            vf0 = 1. - volume_fractions[0]
            x, d = [0.], [fracture_connection_distance]
            z, delta = 1.e-10, 1.e-8
            a = [vf0 * proximity(z) / z]
            volsum = np.cumsum(volume_fractions[1:]) / vf0
            volsum[-1] = 1. - delta
            xl, xr = 0., volume_fractions[1] / a[0]

            for vs in volsum:
                xm, xr = invert_proximity(vs, xl, xr)
                if xm is None:
                    raise Exception("Could not invert MINC proximity function.")
                else:
                    x.append(xm)
                    a.append(vf0 * derivative(proximity, xm, xm * delta))
                    d.append(0.5 * (xm - xl))
                    xl = xm
            d[-1] = inner_dist(x[-2])

            def default_matrix_blockname(blkname, level):
                """Returns default matrix block name, given the original
                block name and the MINC level (> 0)."""
                levelstr = str(level)
                return levelstr + blkname[len(levelstr):]

            def default_minc_rockname(rockname, level):
                """Returns default MINC rocktype name, given the
                original rocktype name and MINC level (>= 0)."""
                if level == 0: return rockname
                else: return 'X' + rockname[1:]

            def duplicate_rock(newrockname, r):
                """Adds new rocktype (if necessary) based on the given one r."""
                if newrockname not in self.rocktype:
                    newrock = rocktype(newrockname, 0, r.density, r.porosity,
                                     r.permeability, r.conductivity, r.specific_heat)
                    self.add_rocktype(newrock)

            # Add MINC blocks and connections:
            if matrix_blockname is None: matrix_blockname = default_matrix_blockname
            if minc_rockname is None: minc_rockname = default_minc_rockname
            blkidict = dict([(blk.name, i) for i, blk in enumerate(self.blocklist)])
            iblk = self.num_blocks - 1
            blockindex = np.zeros((num_levels, len(blocks)), np.int)
            if incon is not None:
                template_vars = incon[0].variable
                newincon = self.incons(template_vars)
                from copy import copy
                for blkinc in incon:
                    newincon[blkinc.block] = copy(incon[blkinc.block])

            for blk_index, blkname in enumerate(blocks):

                blk = self.block[blkname]
                original_vol = blk.volume

                if 0. < original_vol < atmos_volume:

                    blk.volume *= volume_fractions[0]
                    blockindex[0, blk_index] = blkidict[blkname]
                    original_rock = blk.rocktype

                    m, lastblk = 0, blk
                    for vf in volume_fractions[1:]:
                        m += 1
                        mrockname = minc_rockname(original_rock.name, m)
                        duplicate_rock(mrockname, original_rock)
                        mblockname = matrix_blockname(blkname, m)
                        if mblockname in self.block:
                            raise Exception("Duplicate MINC matrix block name: " + mblockname)
                        else:
                            mincblk = t2block(mblockname, original_vol * vf,
                                              self.rocktype[mrockname], centre = blk.centre)
                            self.add_block(mincblk)
                            if incon is not None:
                                inc = copy(incon[blkname])
                                inc.block = mblockname
                                newincon[mblockname] = inc
                            iblk += 1
                            blockindex[m, blk_index] = iblk
                            con = t2connection([lastblk, mincblk], 1, [d[m - 1], d[m]],
                                               original_vol * a[m - 1], None)
                            self.add_connection(con)
                            lastblk = mincblk

                    fract_rockname = minc_rockname(original_rock.name, 0)
                    duplicate_rock(fract_rockname, original_rock)
                    blk.rocktype = self.rocktype[fract_rockname]

            if incon is None: return blockindex
            else:
                newincon.porosity = None
                return blockindex, newincon

    def blockmap(self, geo, index = None):
        """Returns a block mapping from the block name list of the specified
        geometry to the block names in the grid. If the index parameter is
        present (a list of integer indices), the mapping will be to the blocks
        with the specified indices in the grid."""

        if index is None:
            gridblknames = [blk.name for blk in self.blocklist]
        else:
            gridblknames = [self.blocklist[i].name for i in index]
        return dict(zip(geo.block_name_list, gridblknames))

    def reorder(self, block_names = None, connection_names = None, geo = None):
        """Reorders the block (and optionally connection) list to have the
        specified ordering. The block_names parameter should be a list
        of block name strings, while the connection_names parameter
        should be a list of two-element tuples of block name
        strings. If a mulgrid geometry object is passed in using the
        geo parameter, the block and connection name lists from that
        are used (and the block_names and connection_names parameters
        are ignored.)"""

        if geo is not None:
            block_names = geo.block_name_list
            connection_names = geo.block_connection_name_list

        if block_names:
            self.blocklist = [self.block[name] for name in block_names]

        if connection_names:
            connectionlist = []
            for names in connection_names:
                if names in self.connection:
                    connectionlist.append(self.connection[names])
                else:
                    orignames = names[::-1]
                    if orignames in self.connection:
                        con = self.connection[orignames]
                        con.block = con.block[::-1]
                        for blk in con.block:
                            blk.connection_name.remove(orignames)
                            blk.connection_name.add(names)
                        del self.connection[orignames]
                        self.connection[names] = con
                        connectionlist.append(con)
                    else:
                        raise Exception("Unknown connection name: " + names)
            self.connectionlist = connectionlist

    def rename_blocks(self, blockmap = {}, fix_blocknames = True):
        """Rename blocks according to the specified block mapping. The
        connections involving the renamed blocks must also be renamed."""

        if fix_blocknames: mapping = fix_block_mapping(blockmap)

        for blk in self.blocklist:
            name = blk.name
            if name in blockmap:
                del self.block[name]
                mapped_name = blockmap[name]
                self.block[mapped_name] = blk
                blk.name = mapped_name
            cons = set()
            for names in list(blk.connection_name):
                con = []
                for name in names:
                    mapped_name = blockmap[name] if name in blockmap else name
                    con.append(mapped_name)
                cons.add(tuple(con))
            blk.connection_name = cons

        self.connection = {}
        for con in self.connectionlist:
            names = tuple([blk.name for blk in con.block])
            self.connection[names] = con
