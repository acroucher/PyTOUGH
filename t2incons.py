"""For reading and writing TOUGH2 initial conditions files."""

"""
Copyright 2012 University of Auckland.

This file is part of PyTOUGH.

PyTOUGH is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

PyTOUGH is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with PyTOUGH.  If not, see <http://www.gnu.org/licenses/>."""

from mulgrids import *

class t2blockincon(object):
    """Class for a single set of initial conditions at one block."""
    def __init__(self,variable,block='',porosity=None):
        self.block=block
        self.variable=list(variable)
        self.porosity=porosity
    def __getitem__(self,key): return self.variable[key]
    def __setitem__(self,key,value): self.variable[key]=value
    def __repr__(self):
        result=self.block+':'+str(self.variable)
        if self.porosity: result+=' '+str(self.porosity)
        return result

class t2incon(object):
    """Class for a set of initial conditions over a TOUGH2 grid."""
    def __init__(self,filename=''):
        self.empty()
        if filename: self.read(filename)

    def __getitem__(self,key):
        if isinstance(key,(int,slice)): return self._blocklist[key]
        elif isinstance(key,str): return self._block[key]
        else: return None
    def __setitem__(self,key,value):
        if isinstance(value,(list,tuple)): value=t2blockincon(value,key)
        if value.block<>key: value.block=key
        self.add_incon(value)

    def __repr__(self): return 'Initial conditions for '+str(self.num_blocks)+' blocks'

    def empty(self):
        self._block={}
        self._blocklist=[]
        self.timing=None

    def get_num_blocks(self): return len(self._blocklist)
    num_blocks=property(get_num_blocks)

    def get_num_variables(self):
        if self.num_blocks>0: return len(self._blocklist[0].variable)
        else: return 0
    num_variables=property(get_num_variables)

    def get_variable(self):
        """Returns an array of initial condition variables."""
        return np.array([inc.variable for inc in self._blocklist])
    def set_variable(self,val):
        """Sets all initial condition variables to values in an array."""
        for i,inc in enumerate(self._blocklist): inc.variable=val[i]
    variable=property(get_variable,set_variable)

    def get_porosity(self):
        """Returns an array of porosities for each block."""
        return np.array([inc.porosity for inc in self._blocklist])
    def set_porosity(self,por):
        if isinstance(por,np.ndarray):
            if len(por)==self.num_blocks: 
                for i,p in enumerate(por): self._blocklist[i].porosity=p
            else: print 'Porosity array is the wrong length (',len(por),').'
        elif isinstance(por,float) or (por==None):
            for blk in self._blocklist: blk.porosity=por
    porosity=property(get_porosity,set_porosity)

    def get_blocklist(self):
        """Returns an ordered list of blocks."""
        return [inc.block for inc in self._blocklist]
    blocklist=property(get_blocklist)

    def add_incon(self,incon):
        """Adds a t2blockincon."""
        if incon.block in self._block:
            i=self._blocklist.index(self._block[incon.block])
            self._blocklist[i]=incon
        else: self._blocklist.append(incon)
        self._block[incon.block]=incon

    def insert_incon(self,index,incon):
        """Inserts a t2blockincon at the specified index."""
        self._blocklist.insert(index,incon)
        self._block[incon.block]=incon

    def delete_incon(self,block):
        """Deletes a t2blockincon."""
        if block in self._block:
            incon=self._block[block]
            del self._block[block]
            self._blocklist.remove(incon)

    def read(self,filename):
        """Reads initial conditions from file."""
        self.empty()
        f=open(filename,'rU')
        line=f.readline() # header
        finished=False
        timing=False
        while not finished:
            line=f.readline()
            if line.strip():
                if line.startswith('+++'):
                    finished=True
                    timing=True
                else:
                    blkname=fix_blockname(line[0:5])
                    if len(line)>=30: porosity=fortran_float(line[15:30])
                    else: porosity=None
                    line=f.readline()
                    linestrs=[line[i*20:(i+1)*20] for i in xrange(len(line)/20)]
                    vals=tuple([fortran_float(s) for s in linestrs])
                    incon=t2blockincon(vals,blkname,porosity)
                    self.add_incon(incon)
            else: finished=True
        self.timing=None
        if timing:
            line=f.readline().rstrip()
            if len(line)>=45: self.timing={'kcyc':int(line[0:5]),'iter':int(line[5:10]),'nm':int(line[10:15]),
                                  'tstart':float(line[15:30]),'sumtim':float(line[30:45])}

    def write(self,filename,reset=True):
        """Writes initial conditions to file."""
        def list_fmt(x): return ['%20.14e','%20.13e'][abs(x)<1.e-99]
        f=open(filename,'w')
        if (self.timing is None) or reset:
            f.write('INCON\n')
        else:
            f.write('%s%5d%s%12.6e\n' % ('INCON -- INITIAL CONDITIONS FOR',self.num_blocks,' ELEMENTS AT TIME  ',self.timing['sumtim']))
        for incon in self._blocklist:
            f.write('%5s'%unfix_blockname(incon.block))  # TOUGH2 needs mangled block names in the incon file
            if incon.porosity is not None: f.write('          %15.9e'%incon.porosity)
            f.write('\n')
            for v in incon.variable: f.write(list_fmt(v)%v)
            f.write('\n')
        if (self.timing is None) or reset: f.write('\n\n')
        else: f.write('%5d%5d%5d%15.9e%15.9e\n'%(self.timing['kcyc'],self.timing['iter'],self.timing['nm'],
                                                 self.timing['tstart'],self.timing['sumtim']))
        f.close()

    def transfer_from(self,sourceinc,sourcegeo,geo,mapping={},colmapping={}):
        """Transfers initial conditions from another t2incon object, using the two corresponding geometry objects, and the
        optionally specified block and column mappings between them (these are created if not specified).  If there are 
        no atmosphere blocks in the source initial conditions, default atmosphere conditions are assigned if necessary."""
        self.empty()
        if (colmapping=={}) or (mapping=={}): mapping,colmapping=geo.block_mapping(sourcegeo,True)
        from copy import copy
        # atmosphere blocks:
        default_atm_incons=t2blockincon([1.013e5,20.])
        if geo.atmosphere_type==0: # single atmosphere block
            atmblk=geo.block_name(geo.layerlist[0].name,geo.atmosphere_column_name)
            if sourcegeo.atmosphere_type==0: self[atmblk]=copy(sourceinc[0])
            elif sourcegeo.atmosphere_type==1: # take average over source column atmosphere incons
                varsum=np.zeros(len(sourceinc[0].variable))
                for col in sourcegeo.columnlist:
                    blk=sourcegeo.block_name(sourcegeo.layerlist[0].name,col.name)
                    varsum+=np.array(sourceinc[blk].variable)
                self[atmblk]=t2blockincon(varsum/sourcegeo.num_columns)
            else: self[atmblk]=copy(default_atm_incons)
        elif geo.atmosphere_type==1: # atmosphere block over each column
            if sourcegeo.atmosphere_type==0: # broadcast single source atmosphere incons to each column
                for col in geo.columnlist:
                    blk=geo.block_name(geo.layerlist[0].name,col.name)
                    self[blk]=copy(sourceinc[0])
            elif sourcegeo.atmosphere_type==1: # atmosphere over each column in both source and destination
                for col in geo.columnlist:
                    mappedcol=colmapping[col.name]
                    oldatmosblockname=sourcegeo.block_name(sourcegeo.layerlist[0].name,mappedcol)
                    blk=geo.block_name(geo.layerlist[0].name,col.name)
                    self[blk]=copy(sourceinc[oldatmosblockname])
            else: # use default
                for col in geo.columnlist:
                    blk=geo.block_name(geo.layerlist[0].name,col.name)
                    self[blk]=copy(default_atm_incons)
        # underground blocks:
        for blk in geo.block_name_list[geo.num_atmosphere_blocks:]: self[blk]=copy(sourceinc[mapping[blk]])
