"""For manipulating TOUGH2 grids."""

"""
Copyright 2011 University of Auckland.

This file is part of PyTOUGH.

PyTOUGH is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

PyTOUGH is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with PyTOUGH.  If not, see <http://www.gnu.org/licenses/>."""

from mulgrids import *
from t2incons import *

class rocktype(object):
    """Rock type"""
    def __init__(self,name="dfalt",nad=0,density=2600.0,porosity=0.1,permeability=1.0e-15*np.ones(3),conductivity=1.5,specific_heat=900.0):
        self.name=name
        self.nad=nad
        self.density=density
        self.porosity=porosity
        if isinstance(permeability,list): permeability=np.array(permeability)
        self.permeability=permeability
        self.conductivity=conductivity
        self.specific_heat=specific_heat
        self.compressibility=0.0
        self.expansivity=0.0
        self.dry_conductivity=0.0
        self.tortuosity=0.0
        self.relative_permeability={}
        self.capillarity={}
    def __repr__(self):
        return self.name

class t2block(object):
    """Grid block"""
    def __init__(self,name='     ',volume=1.0,blockrocktype=rocktype(),centre=None,atmosphere=False,ahtx=None,pmx=None,nseq=None,nadd=None):
        self.name=name
        self.volume=volume
        self.rocktype=blockrocktype
        if isinstance(centre,list): centre=np.array(centre)
        self.centre=centre
        self.atmosphere=atmosphere
        self.ahtx=ahtx
        self.pmx=pmx
        self.nseq,self.nadd=nseq,nadd
        self.connection_name=set([])
    def __repr__(self): return self.name
    def get_num_connections(self): return len(self.connection_name)
    num_connections=property(get_num_connections)
    def get_neighbour_names(self):
        """Returns a set of neighbouring block names- those connected to this one."""
        return set([[blkname for blkname in cn if blkname<>self.name][0] for cn in self.connection_name])
    neighbour_name=property(get_neighbour_names)

class t2connection(object):
    """Connection between two blocks"""
    def __init__(self,blocks=[t2block(),t2block()],direction=1,distance=[0.0,0.0],area=1.0,dircos=0.0,sigma=None,nseq=None,nad1=None,nad2=None):
        self.block=blocks
        self.direction=direction # permeability direction
        self.distance=distance
        self.area=area
        self.dircos=dircos # direction cosine
        self.sigma=sigma # radiant emittance factor (TOUGH2)
        self.nseq,self.nad1,self.nad2=nseq,nad1,nad2
        self.centre = None
        self.midpoint = None
        self.normal = None
    def __repr__(self):
        return self.block[0].name+':'+self.block[1].name

class t2grid(object):
    """TOUGH2 grid"""
    def __init__(self): self.empty()

    def get_num_rocktypes(self):
        return len(self.rocktypelist)
    num_rocktypes=property(get_num_rocktypes)
        
    def get_num_blocks(self):
        return len(self.blocklist)
    num_blocks=property(get_num_blocks)
        
    def get_num_connections(self):
        return len(self.connectionlist)
    num_connections=property(get_num_connections)

    def get_num_atmosphere_blocks(self):
        return len(self.atmosphere_blocks)
    num_atmosphere_blocks=property(get_num_atmosphere_blocks)

    def get_num_underground_blocks(self):
        return self.num_blocks-self.num_atmosphere_blocks
    num_underground_blocks=property(get_num_underground_blocks)
    
    def get_atmosphere_blocks(self):
        return [blk for blk in self.blocklist if blk.atmosphere]
    atmosphere_blocks=property(get_atmosphere_blocks)

    def get_block_centres_defined(self):
        if self.num_atmosphere_blocks==1: istart=1
        else: istart=0
        return any([blk.centre is not None for blk in self.blocklist[istart:]])
    block_centres_defined=property(get_block_centres_defined)

    def get_connection_centres_defined(self):
        return any([con.centre is not None for con in self.connectionlist])
    connection_centres_defined = property(get_connection_centres_defined)

    def calculate_block_centres(self,geo):
        """Calculates block centres from geometry object."""
        if geo.atmosphere_type==0:
            istart=1
            self.blocklist[0].centre=None  # centre not well defined for single atmosphere block
        else: istart=0
        for blk in self.blocklist[istart:]:
            layername=geo.layer_name(blk.name)
            lay=geo.layer[layername]
            colname=geo.column_name(blk.name)
            col=geo.column[colname]
            blk.centre=geo.block_centre(lay,col)

    def calculate_connection_centres(self, geo):
        """Calculates centre (and normal vector) for each connection face.  Note that the 'centre'
        depends on which block the connection is approached from- in the case where the connection
        face is not orthogonal to the line between the block centres.  Hence there are two 'centres'.
        The mipoint is just midway between the connection face nodes."""
        layindex = dict([(lay.name,i) for i,lay in enumerate(geo.layerlist)])
        for con in self.connectionlist:
            con.centre = {}
            layernames = [geo.layer_name(blk.name) for blk in con.block]
            if layernames[0] == layernames[1]: # horizontal connection
                lay = geo.layer[layernames[0]]
                colnames = tuple([geo.column_name(blk.name) for blk in con.block])
                geocon = geo.connection[colnames]
                nodepos = [node.pos for node in geocon.node]
                dpos = nodepos[1] - nodepos[0]
                for colname,blk in zip(colnames, con.block):
                    col = geo.column[colname]
                    if blk.centre is not None: vcentre = blk.centre[2]
                    else: vcentre = geo.block_centre(lay,col)[2]
                    hcentre = line_projection(col.centre, nodepos)
                    con.centre[blk.name] = np.hstack((hcentre, np.array([vcentre])))
                con.normal = np.array([dpos[1], -dpos[0], 0.]) / np.linalg.norm(dpos)
                con.midpoint = np.hstack((0.5 * sum(nodepos), min([centre[2] for centre in con.centre])))
            else: # vertical connection
                layindices = np.array([layindex[layname] for layname in layernames])
                ilower = np.argmax(layindices)
                colname = geo.column_name(con.block[ilower].name)
                col = geo.column[colname]
                hcentre = col.centre
                lay = geo.layer[layernames[ilower]]
                vcentre = geo.block_surface(lay,col)
                sgn = [1.,-1.][ilower]
                con.normal = np.array([0., 0., sgn])
                con.midpoint = np.hstack((hcentre, np.array([vcentre])))
                for blk in con.block: con.centre[blk.name] = con.midpoint

    def rocktype_frequency(self,rockname):
        """Returns how many times the rocktype with given name is used in the grid."""
        return [blk.rocktype.name for blk in self.blocklist].count(rockname)
    def get_rocktype_frequencies(self):
        """Returns a list of tuples of occurring frequencies of rock types in the grid and the names of rocktypes with that frequency, 
        ordered by increasing frequency."""
        freq=[(rt.name,self.rocktype_frequency(rt.name)) for rt in self.rocktypelist]
        occurring_freqs=list(set([item[1] for item in freq]))
        occurring_freqs.sort()
        frocks=dict([(f,[]) for f in occurring_freqs])
        for item in freq: frocks[item[1]].append(item[0])
        return [(f,frocks[f]) for f in occurring_freqs]
    rocktype_frequencies=property(get_rocktype_frequencies)

    def sort_rocktypes(self):
        """Sorts rocktype list in alphabetical order by name."""
        rocknames=[rt.name for rt in self.rocktypelist]
        rocknames.sort()
        self.rocktypelist=[self.rocktype[name] for name in rocknames]

    def __repr__(self):
        return str(self.num_rocktypes)+' rock types; '+str(self.num_blocks)+' blocks; '+str(self.num_connections)+' connections'

    def __add__(self,other):
        """Adds two grids together."""
        result=t2grid()
        for grid in [self,other]:
            for rt in grid.rocktypelist: result.add_rocktype(rt)
            for blk in grid.blocklist: result.add_block(blk)
            for con in grid.connectionlist: result.add_connection(con)
        return result

    def embed(self,subgrid,connection):
        """Returns a grid with a subgrid embedded inside one of its blocks.  The connection specifies how the two grids
        are to be connected: the blocks to be connected and the connection distances, area etc. between them.  The first 
        block should be the host block, the second the connecting block in the subgrid."""
        result=None
        subvol=sum([blk.volume for blk in subgrid.blocklist])
        hostblock=connection.block[0]
        if subvol<hostblock.volume:
            dupblks=set([blk.name for blk in self.blocklist]) & set([blk.name for blk in subgrid.blocklist])
            if len(dupblks)==0:
                result=self+subgrid
                connection.block=[result.block[blk.name] for blk in connection.block]
                result.add_connection(connection)
                result.block[hostblock.name].volume-=subvol # remove subgrid volume from host block
            else: print 'Grid embedding error: the following blocks are in both grids:',dupblks
        else: print 'Grid embedding error: the host block is not big enough to contain the subgrid.'
        return result

    def empty(self):
        """Empties a TOUGH2 grid"""
        self.rocktypelist=[]
        self.blocklist=[]
        self.connectionlist=[]
        self.rocktype={}
        self.block={}
        self.connection={}
        
    def add_rocktype(self,newrocktype=rocktype()):
        """Adds a rock type to the grid.  Any existing rocktype of the same name is replaced."""
        if newrocktype.name in self.rocktype:
            i=self.rocktypelist.index(self.rocktype[newrocktype.name])
            self.rocktypelist[i]=newrocktype
        else: self.rocktypelist.append(newrocktype)
        self.rocktype[newrocktype.name]=newrocktype

    def delete_rocktype(self,rocktypename):
        """Deletes a rock type from the grid"""
        if rocktypename in self.rocktype:
            rt=self.rocktype[rocktypename]
            del self.rocktype[rocktypename]
            self.rocktypelist.remove(rt)

    def clean_rocktypes(self):
        """Deletes any unused rock types from the grid"""
        unused_rocktypes=[]
        for rt in self.rocktypelist:
            if self.rocktype_frequency(rt.name)==0: unused_rocktypes.append(rt.name)
        for name in unused_rocktypes: self.delete_rocktype(name)

    def add_block(self,newblock=t2block()):
        """Adds a block to the grid"""
        if newblock.name in self.block:
            i=self.blocklist.index(self.block[newblock.name])
            self.blocklist[i]=newblock
        else: self.blocklist.append(newblock)
        self.block[newblock.name]=newblock
    
    def delete_block(self,blockname):
        """Deletes a block from the grid"""
        if blockname in self.block:
            blk=self.block[blockname]
            from copy import copy
            connames=copy(blk.connection_name)
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

    def add_connection(self,newconnection=t2connection()):
        """Adds a connection to the grid"""
        conname=tuple([blk.name for blk in newconnection.block])
        if conname in self.connection: 
            i=self.connectionlist.index(self.connection[conname])
            self.connectionlist[i]=newconnection
        else: self.connectionlist.append(newconnection)
        self.connection[conname]=newconnection
        for block in newconnection.block: block.connection_name.add(conname)

    def delete_connection(self,connectionname):
        """Deletes a connection from the grid"""
        if connectionname in self.connection:
            con=self.connection[connectionname]
            for block in con.block: block.connection_name.remove(connectionname)
            del self.connection[connectionname]
            self.connectionlist.remove(con)

    def block_index(self,blockname):
        """Returns index of block with specified name in the block list of the grid"""
        if blockname in self.block:
            return self.blocklist.index(self.block[blockname])
        else: return None

    def connection_index(self,connectionnames):
        """Returns index of connection with specified pair of names in the connection list of the grid"""
        if connectionnames in self.connection:
            return self.connectionlist.index(self.connection[connectionnames])
        else: return None

    def fromgeo(self,geo):
        """Converts a MULgraph grid to a TOUGH2 grid"""
        self.empty()
        self.add_rocktype(rocktype()) # add default rock type
        self.add_blocks(geo)
        self.add_connections(geo)
        return self

    def add_blocks(self,geo):
        """Adds blocks to grid from MULgraph geometry file"""
        self.add_atmosphereblocks(geo)
        self.add_underground_blocks(geo)

    def add_atmosphereblocks(self,geo):
        """Adds atmosphere blocks from geometry"""
        atmosrocktype=self.rocktypelist[0]
        if geo.atmosphere_type==0: # one atmosphere block
            atmblockname=geo.block_name(geo.layerlist[0].name,geo.atmosphere_column_name)
            centre=None
            self.add_block(t2block(atmblockname,geo.atmosphere_volume,atmosrocktype,centre=centre,atmosphere=True))
        elif geo.atmosphere_type==1: # one atmosphere block per column
            for col in geo.columnlist:
                atmblockname=geo.block_name(geo.layerlist[0].name,col.name)
                centre=geo.block_centre(geo.layerlist[0],col)
                self.add_block(t2block(atmblockname,geo.atmosphere_volume,atmosrocktype,centre=centre,atmosphere=True))

    def add_underground_blocks(self,geo):
        """Add underground blocks from geometry"""
        for lay in geo.layerlist[1:]:
            for col in [col for col in geo.columnlist if col.surface>lay.bottom]:
                name=geo.block_name(lay.name,col.name)
                centre=geo.block_centre(lay,col)
                self.add_block(t2block(name,geo.block_volume(lay,col),self.rocktypelist[0],centre=centre))

    def add_connections(self,geo):
        """Add connections from geometry"""
        for lay in geo.layerlist[1:]:
            layercols=[col for col in geo.columnlist if col.surface>lay.bottom]
            self.add_vertical_layer_connections(geo,lay,layercols)
            self.add_horizontal_layer_connections(geo,lay,layercols)

    def add_vertical_layer_connections(self,geo,lay,layercols=[]):
        """Add vertical connections in layer"""
        for col in layercols:
            thisblk=self.block[geo.block_name(lay.name,col.name)]
            if (geo.layerlist.index(lay)==1) or (col.surface<=lay.top): # connection to atmosphere
                abovelayer=geo.layerlist[0]
                abovedist=geo.atmosphere_connection
                belowdist=col.surface-thisblk.centre[2]
                if geo.atmosphere_type==0:
                    aboveblk=self.blocklist[0]
                elif geo.atmosphere_type==1:
                    aboveblk=self.block[geo.block_name(abovelayer.name,col.name)]
                else: # no atmosphere blocks
                    continue
            else:
                ilayer=geo.layerlist.index(lay)
                abovelayer=geo.layerlist[ilayer-1]
                aboveblk=self.block[geo.block_name(abovelayer.name,col.name)]
                abovedist=aboveblk.centre[2]-abovelayer.bottom
                belowdist=lay.top-lay.centre
            con=t2connection([thisblk,aboveblk],3,[belowdist,abovedist],col.area,-1.0)
            self.add_connection(con)

    def add_horizontal_layer_connections(self,geo,lay,layercols=[]):
        """Add horizontal connections in layer"""
        from math import cos,sin
        layercolset=set(layercols)
        anglerad=geo.permeability_angle*np.pi/180.
        c,s=cos(anglerad),sin(anglerad)
        rotation=np.array([[c,s],[-s,c]])
        for con in [con for con in geo.connectionlist if set(con.column).issubset(layercolset)]:
            conblocks=[self.block[geo.block_name(lay.name,concol.name)] for concol in con.column]
            [dist,area]=geo.connection_params(con,lay)
            d=conblocks[1].centre-conblocks[0].centre
            d2=np.dot(rotation,d[0:2])
            direction=np.argmax(abs(d2))+1
            dircos=-d[2]/np.linalg.norm(d)
            self.add_connection(t2connection(conblocks,direction,dist,area,dircos))

    def copy_connection_directions(self,geo,grid):
        """Copies connection permeability directions from another grid.  It is assumed the two grids have
        the same column structure.  The geo argument is the geometry file corresponding to grid."""
        nlayercons=len(geo.connectionlist)
        noldcons=len(grid.connectionlist)
        # create dictionary of permeability directions within a layer (from bottom layer, assumed complete):
        dirn={}
        for i,con in enumerate(geo.connectionlist):
            dirn[(con.column[0].name,con.column[1].name)]=grid.connectionlist[noldcons-nlayercons+i].direction
        # transfer permeability directions to horizontal connections:
        for con in self.connectionlist:
            if con.direction<3:
                colnames=tuple([geo.column_name(blk.name) for blk in con.block])
                con.direction=dirn[colnames]

    def get_unconnected_blocks(self):
        """Returns a set of blocks in the grid that are not connected to any other blocks."""
        return set([blk.name for blk in self.blocklist if len(blk.connection_name)==0])
    unconnected_blocks=property(get_unconnected_blocks)
    
    def get_isolated_rocktype_blocks(self):
        """Returns a list of blocks with isolated rocktypes- that is, blocks with a rocktype different from that of
        all other blocks they are connected to."""
        bc_volume=1.e20  # blocks with volume greater than this are considered boundary condition blocks and not counted
        return set([blk.name for blk in self.blocklist if (blk.volume<bc_volume) and not (blk.rocktype.name in [self.block[nbr].rocktype.name for nbr in blk.neighbour_name])])
    isolated_rocktype_blocks=property(get_isolated_rocktype_blocks)

    def check(self,fix=False,silent=False):
        """Checks a grid for errors, and optionally fixes them.  Errors checked for are:
        - blocks not connected to any other blocks
        - blocks with isolated rocktypes
        Returns True if no errors were found, and False otherwise.  If silent is True, there is no printout.
        Unconnected blocks are fixed by deleting them.  Isolated rocktype blocks are fixed by assigning them the
        most popular rocktype of their neighbours."""
        ok=True
        ub=self.unconnected_blocks
        if len(ub)>0:
            ok=False
            if not silent: print 'Unconnected blocks:',list(ub)
            if fix:
                for blk in ub: self.delete_block(blk)
                if not silent: print 'Unconnected blocks fixed.'
        ib=self.isolated_rocktype_blocks
        if len(ib)>0:
            ok=False
            if not silent: print 'Isolated rocktype blocks:',list(ib)
            if fix:
                for blk in ib:
                    nbr_rocktype=[self.block[nbr].rocktype.name for nbr in self.block[blk].neighbour_name]
                    pop_rocktype=max(set(nbr_rocktype), key=nbr_rocktype.count)
                    self.block[blk].rocktype=self.rocktype[pop_rocktype]
                if not silent: print 'Isolated rocktype blocks fixed.'
        if ok and not silent: print 'No problems found.'
        return ok

    def get_rocktype_indices(self):
        """Returns an integer array containing the rocktype index for each block in the grid."""
        rocknames=[rt.name for rt in self.rocktypelist]
        rockdict=dict([(name,i) for i,name in enumerate(rocknames)])
        return np.array([rockdict[blk.rocktype.name] for blk in self.blocklist])
    rocktype_indices=property(get_rocktype_indices)

    def get_vtk_data(self,geo):
        """Returns dictionary of VTK data arrays from rock types.  The geometry object geo must be passed in."""
        from vtk import vtkIntArray,vtkFloatArray,vtkCharArray
        arrays={'Block':{'Rock type index':vtkIntArray(),'Porosity':vtkFloatArray(),
                         'Permeability':vtkFloatArray(),'Name':vtkCharArray()},'Node':{}}
        vector_properties=['Permeability']
        string_properties=['Name']
        string_length=5
        nele=geo.num_underground_blocks
        array_length={'Block':nele,'Node':0}
        for array_type,array_dict in arrays.items():
            for name,array in array_dict.items():
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
        natm=geo.num_atmosphere_blocks
        rindex=self.rocktype_indices[natm:]
        for i,ri in enumerate(rindex):
            arrays['Block']['Rock type index'].SetValue(i,ri)
            rt=self.rocktypelist[ri]
            arrays['Block']['Porosity'].SetValue(i,rt.porosity)
            k=rt.permeability
            arrays['Block']['Permeability'].SetTuple3(i,k[0],k[1],k[2])
        for i,blk in enumerate(self.blocklist[natm:]):
            arrays['Block']['Name'].SetTupleValue(i,blk.name)
        return arrays

    def write_vtk(self,geo,filename,wells=False):
        """Writes *.vtu file for a vtkUnstructuredGrid object corresponding to the grid in 3D, with the specified filename,
        for visualisation with VTK."""
        from vtk import vtkXMLUnstructuredGridWriter
        if wells: geo.write_well_vtk()
        arrays=geo.vtk_data
        grid_arrays=self.get_vtk_data(geo)
        for array_type,array_dict in arrays.items():
            array_dict.update(grid_arrays[array_type])
        vtu=geo.get_vtk_grid(arrays)
        writer=vtkXMLUnstructuredGridWriter()
        writer.SetFileName(filename)
        writer.SetInput(vtu)
        writer.Write()

    def flux_matrix(self,geo):
        """Returns a sparse matrix which can be used to multiply a vector of connection table values for underground
        blocks, to give approximate average fluxes of those values at the block centres."""
        natm = geo.num_atmosphere_blocks
        nele = geo.num_underground_blocks
        conindex = dict([((c.block[0].name,c.block[1].name),i) for i,c in enumerate(self.connectionlist)])
        from scipy import sparse
        A = sparse.lil_matrix((3*nele, self.num_connections))
        if not self.block_centres_defined: self.calculate_block_centres(geo)
        if not self.connection_centres_defined: self.calculate_connection_centres(geo)
        for iblk,blk in enumerate(self.blocklist[natm:]):
            ncons = blk.num_connections
            if ncons > 0:
                M,icons = [],[]
                for conname in blk.connection_name:
                    con = self.connection[conname]
                    anormal = con.normal * con.area
                    row = list(anormal) # fit constant flows
                    if ncons >= 6: row += list((con.centre[blk.name] - blk.centre) * anormal) # fit linear flows
                    M.append(row)
                    icons.append(conindex[conname])
                Ablk = -np.linalg.pinv(np.array(M))
                ib3 = iblk*3
                for ic in xrange(ncons):
                    for ip in xrange(3): A[ib3+ip,icons[ic]] = Ablk[ip,ic]
        return A

    def radial(self,rblocks,zblocks,convention=0,atmos_type=2,origin=np.array([0.,0.]),justify='r',case='l',dimension=2):
        """Returns a radial TOUGH2 grid with the specified radial and vertical block sizes.
        The arguments are arrays of the block sizes in each dimension (r,z).
        Naming convention, atmosphere type and grid origin can optionally be specified.  The origin is in 
        (r,z) coordinates, so origin[0] is the starting radius of the grid.  (The origin can also be specified
        with three components, in which case the second one is ignored.)
        The optional justify and case parameters specify the format of the character part of the block names
        (whether they are right or left justified, and lower or upper case).
        Specifying dimension<>2 (between 1 and 3) simulates flow in fractured rock using the
        "generalized radial flow" concept of Barker, J.A. (1988), "A generalized radial flow model for hydraulic
        tests in fractured rock", Water Resources Research 24(10), 1796-1804.  In this case it probably doesn't
        make much sense to have more than one block in the z direction. """

        if isinstance(rblocks,list): rblocks=np.array(rblocks)
        if isinstance(zblocks,list): zblocks=np.array(zblocks)
        if isinstance(origin,list): origin=np.array(origin)
        if len(origin)>2: origin=origin[[0,2]]

        from string import ljust,rjust,lowercase,uppercase
        justfn=[rjust,ljust][justify=='l']
        casefn=[uppercase,lowercase][case=='l']

        n2=0.5*dimension
        if dimension<>2: # need gamma function
            try:
                from math import gamma # included for Python 2.7 or later
            except ImportError:
                from scipy.special import gamma
            gamman2=gamma(n2)
        else: gamman2=1.0
        alpha=2./gamman2*np.pi**n2

        b=np.sum(zblocks) # total thickness
        b2n=b**(2-dimension)
        r=origin[0]+np.concatenate((np.zeros(1),np.cumsum(rblocks)))
        rin,rout=r[:-1],r[1:] # inner and outer radii
        rc=0.5*(rin+rout) # centre radius
        A=alpha*b2n/dimension*np.diff(r**dimension) # "top area"
        c=alpha*b2n*rout**(dimension-1) # "outer circumference"
        ncols=len(rblocks)

        grid=t2grid()
        grid.add_rocktype(rocktype()) # add default rock type

        # dummy geometry for creating block names etc:
        geo=mulgrid(type='GENER',convention=convention,atmos_type=atmos_type)
        for ir,dr in enumerate(rblocks):
            colname=geo.column_name_from_number(ir+1,justfn,casefn)
            geo.add_column(column(colname,[],centre=np.array([rc[ir],0.])))
        geo.add_layers(zblocks,origin[1],justify,case)
        grid.add_atmosphereblocks(geo)

        for lay in geo.layerlist[1:]: # add blocks
            V=A*lay.thickness
            for col,rcentre,vol in zip(geo.columnlist,rc,V):
                name=geo.block_name(lay.name,col.name)
                centre=np.array([rcentre,0.,lay.centre])
                grid.add_block(t2block(name,vol,grid.rocktypelist[0],centre=centre))

        for ilay,lay in enumerate(geo.layerlist[1:]):
            Ar=c*lay.thickness
            for icol,col in enumerate(geo.columnlist): # vertical connections
                top_area=A[icol]
                blkindex=ilay*ncols+icol+geo.num_atmosphere_blocks
                thisblk=grid.blocklist[blkindex]
                if ilay==0: # atmosphere connections
                    abovedist=geo.atmosphere_connection
                    belowdist=0.5*geo.layerlist[1].thickness
                    if atmos_type==0: aboveblk=grid.blocklist[0]
                    elif atmos_type==1: aboveblk=grid.blocklist[icol]
                    else: continue
                else:
                    abovelayer=geo.layerlist[ilay]
                    aboveblk=grid.blocklist[blkindex-ncols]
                    abovedist=aboveblk.centre[2]-abovelayer.bottom
                    belowdist=lay.top-lay.centre
                con=t2connection([thisblk,aboveblk],3,[belowdist,abovedist],top_area,-1.0)
                grid.add_connection(con)
            for icol,col in enumerate(geo.columnlist[:-1]): # radial connections
                nextcol=geo.columnlist[icol+1]
                conblocks=[grid.block[geo.block_name(lay.name,acol.name)] for acol in [col,nextcol]]
                dist,area=[0.5*rblocks[icol],0.5*rblocks[icol+1]],Ar[icol]
                direction,dircos=1,0.0
                grid.add_connection(t2connection(conblocks,direction,dist,area,dircos))

        return grid

    def incons(self,values=(101.3e3,20.)):
        """Creates a t2incon initial condtions object corresponding to the grid from the given values.  If initial
        conditions are given for one block only, these are applied to all blocks."""
        inc=t2incon()
        values=np.array(values)
        if len(np.shape(values))==1:
            from copy import copy
            for blk in self.blocklist: inc[blk.name]=copy(tuple(values))
        else:
            for blk,val in zip(self.blocklist,values): inc[blk.name]=tuple(val)
        return inc

    def neargroups(self,blocknames):
        """Given a list or set of block names, finds groups of 'near' blocks.  Blocks are assigned the same group
        if they are neighbours, or share a neighbour."""
        blocknames=list(set(blocknames))
        groups=[]
        for blk in blocknames: groups.append(set([blk]))
        from copy import copy
        done=False
        while not done:
            done=True
            for i,g in enumerate(groups):
                ng=copy(g)
                for blk in g: ng = ng | self.block[blk].neighbour_name
                if i<len(groups)-1:
                    for g2 in groups[i+1:]:
                        ng2=copy(g2)
                        for blk in g2: ng2 = ng2 | self.block[blk].neighbour_name
                        if ng & ng2:
                            g.update(g2)
                            groups.remove(g2)
                            done=False
                            break
                    if not done: break
        return groups
