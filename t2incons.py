"""For reading and writing TOUGH2 initial conditions files.

Copyright 2012 University of Auckland.

This file is part of PyTOUGH.

PyTOUGH is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

PyTOUGH is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with PyTOUGH.  If not, see <http://www.gnu.org/licenses/>."""

from mulgrids import *
from fixed_format_file import *

t2incon_format_specification = {
    'header_short': [['title'], ['5s']],
    'header_long': [['s1', 'nele', 's2', 'sumtim'], ['31s', '5d', '19s', '12.6e']],
    'incon1': [['name', 'nseq', 'nadd', 'porx'],
               ['5s', '5d', '5d', '15.9e']],
    'incon1_toughreact': [['name', 'nseq', 'nadd', 'porx', 'k1', 'k2', 'k3'],
               ['5s', '5d', '5d'] + ['15.9e'] * 4],
    'incon2': [['x1', 'x2', 'x3', 'x4'], ['20.13e'] * 4],
    'timing': [['kcyc', 'iter', 'nm', 'tstart', 'sumtim'],
               ['5d'] * 3 + ['15.9e'] * 2],
    'timing_toughreact': [['kcyc', 'iter', 'nm', 'tstart', 'sumtim'],
               ['6d'] * 2 + ['3d'] + ['15.9e'] * 2]
}

class t2incon_parser(fixed_format_file):
    """Class for parsing TOUGH2 incon file."""
    def __init__(self, filename, mode, read_function = fortran_read_function):
        super(t2incon_parser,self).__init__(filename, mode,
                                            t2incon_format_specification, read_function)

class t2blockincon(object):
    """Class for a single set of initial conditions at one block."""
    def __init__(self, variable, block = '', porosity = None, permeability = None,
                 nseq = None, nadd = None):
        self.block = block
        self.variable = list(variable)
        self.porosity = porosity
        self.permeability = permeability
        self.nseq, self.nadd = nseq, nadd
    def __getitem__(self, key): return self.variable[key]
    def __setitem__(self, key, value): self.variable[key] = value
    def __repr__(self):
        result = self.block + ':' + str(self.variable)
        if self.porosity is not None:
            result += ' ' + str(self.porosity)
        if self.permeability is not None:
            result += ' ' + str(self.permeability)
        if self.nseq is not None:
            result += ' (' + str(self.nseq) + ', ' + str(self.nadd) + ')'
        return result

class t2incon(object):
    """Class for a set of initial conditions over a TOUGH2 grid."""
    def __init__(self, filename = '',
                 read_function = fortran_read_function, num_variables = None):
        self.simulator = 'TOUGH2'
        self.read_function = read_function
        self.empty()
        if filename: self.read(filename, num_variables)

    def __getitem__(self, key):
        if isinstance(key, (int, slice)): return self._blocklist[key]
        elif isinstance(key, str): return self._block[key]
        else: return None
    def __setitem__(self, key, value):
        if isinstance(value,(list,tuple)):
            value = t2blockincon(value,key)
        if value.block != key: value.block = key
        self.add_incon(value)

    def __repr__(self):
        return self.simulator + ' initial conditions for ' + \
            str(self.num_blocks) + ' blocks'

    def empty(self):
        self._block = {}
        self._blocklist = []
        self.timing = None

    def get_num_blocks(self): return len(self._blocklist)
    num_blocks = property(get_num_blocks)

    def get_num_variables(self):
        if self.num_blocks > 0: return len(self._blocklist[0].variable)
        else: return 0
    num_variables = property(get_num_variables)

    def get_variable(self):
        """Returns an array of initial condition variables."""
        return np.array([inc.variable for inc in self._blocklist])
    def set_variable(self, val):
        """Sets all initial condition variables to values in an array."""
        for i, inc in enumerate(self._blocklist): inc.variable = val[i]
    variable = property(get_variable, set_variable)

    def get_porosity(self):
        """Returns an array of porosities for each block."""
        return np.array([inc.porosity for inc in self._blocklist])
    def set_porosity(self, por):
        if isinstance(por, np.ndarray):
            if len(por) == self.num_blocks:
                for i, p in enumerate(por): self._blocklist[i].porosity = p
            else:
                raise Exception('Porosity array is the wrong length (' + \
                                str(len(por)) + ').')
        elif isinstance(por, float) or (por == None):
            for blk in self._blocklist: blk.porosity = por
    porosity = property(get_porosity, set_porosity)

    def get_permeability(self):
        """Returns an array of permeabilities for each block."""
        return np.array([inc.permeability for inc in self._blocklist])
    def set_permeability(self, perm):
        if isinstance(perm, np.ndarray):
            shape = np.shape(perm)
            if shape == (3,):
                from copy import copy
                for blk in self._blocklist: blk.permeability = copy(perm)
            elif shape == (self.num_blocks, 3):
                for i, p in enumerate(perm): self._blocklist[i].permeability = p
            else:
                raise Exception('Permeability array is the wrong shape (' + \
                                str(shape) + ').')
        elif isinstance(perm, float):
            for blk in self._blocklist:
                blk.permeability = perm * np.ones(3)
        elif perm == None:
            for blk in self._blocklist:
                blk.permeability = None

    porosity = property(get_porosity, set_porosity)
    permeability = property(get_permeability, set_permeability)

    def get_blocklist(self):
        """Returns an ordered list of blocks."""
        return [inc.block for inc in self._blocklist]
    blocklist = property(get_blocklist)

    def add_incon(self, incon):
        """Adds a t2blockincon."""
        if incon.block in self._block:
            i = self._blocklist.index(self._block[incon.block])
            self._blocklist[i] = incon
        else: self._blocklist.append(incon)
        self._block[incon.block] = incon

    def insert_incon(self, index, incon):
        """Inserts a t2blockincon at the specified index."""
        self._blocklist.insert(index, incon)
        self._block[incon.block] = incon

    def delete_incon(self, block):
        """Deletes a t2blockincon."""
        if block in self._block:
            incon = self._block[block]
            del self._block[block]
            self._blocklist.remove(incon)

    def read(self, filename, num_variables = None):
        """Reads initial conditions from file."""
        self.empty()
        infile = t2incon_parser(filename, 'rU', read_function = self.read_function)
        infile.readline() # skip header
        finished = False
        timing = False
        while not finished:
            line = infile.readline()
            if line.strip():
                if line.startswith('+++'):
                    finished = True
                    timing = True
                else:
                    line = padstring(line)
                    [blkname, nseq, nadd, porosity, k1, k2, k3] = \
                        infile.parse_string(line, 'incon1_toughreact')
                    if valid_blockname(blkname):
                        blkname = fix_blockname(blkname)
                        if (k1 is None or k2 is None or k3 is None): permeability = None
                        else:
                            permeability = np.array([k1, k2, k3])
                            self.simulator = 'TOUGHREACT'
                        vals, more = [], True
                        while more:
                            linevals = infile.read_values('incon2')
                            while linevals and linevals[-1] is None: linevals.pop()
                            vals += linevals
                            more = False if num_variables is None else \
                                   len(vals) < num_variables
                        incon = t2blockincon(vals, blkname, porosity, permeability,
                                             nseq, nadd)
                        self.add_incon(incon)
                    else:
                        raise Exception('Invalid block name (' + blkname + \
                                        ') in incon file: ' + filename)
            else: finished = True
        self.timing = None
        if timing:
            line = infile.readline()
            if line.strip():
                line = padstring(line)
                timing_fmt = 'timing'
                if self.simulator == 'TOUGHREACT': timing_fmt += '_toughreact'
                [kcyc, itr, nm, tstart, sumtim] = infile.parse_string(line, timing_fmt)
                self.timing = {'kcyc': kcyc, 'iter': itr, 'nm': nm,
                               'tstart': tstart, 'sumtim': sumtim}

    def write(self, filename, reset = True):
        """Writes initial conditions to file."""
        outfile = t2incon_parser(filename, 'w')
        if (self.timing is None) or reset:
            outfile.write_values(['INCON'], 'header_short')
        else:
            outfile.write_values(['INCON -- INITIAL CONDITIONS FOR', self.num_blocks,
                                  ' ELEMENTS AT TIME  ', self.timing['sumtim']], 'header_long')
        for incon in self._blocklist:
            blkname = unfix_blockname(incon.block)
            if self.simulator == 'TOUGHREACT' and incon.permeability is not None:
                outfile.write_values([blkname, incon.nseq, incon.nadd, incon.porosity] +
                                      list(incon.permeability), 'incon1_toughreact')
            else:
                outfile.write_values([blkname, incon.nseq, incon.nadd, incon.porosity], 'incon1')
            vals = list(incon.variable)
            while vals:
                linelen = min(len(vals), 4)
                linevals = vals[:linelen]
                outfile.write_values(linevals, 'incon2')
                vals = vals[linelen:]
        if (self.timing is None) or reset: outfile.write('\n\n')
        else:
            outfile.write('+++\n')
            timing_fmt = 'timing'
            if self.simulator == 'TOUGHREACT': timing_fmt += '_toughreact'
            outfile.write_value_line(self.timing, timing_fmt)
        outfile.close()

    def transfer_from(self, sourceinc, sourcegeo, geo, mapping = {}, colmapping = {}):
        """Transfers initial conditions from another t2incon object, using the
        two corresponding geometry objects, and the optionally
        specified block and column mappings between them (these are
        created if not specified).  If there are no atmosphere blocks
        in the source initial conditions, default atmosphere
        conditions are assigned if necessary.
        """
        self.empty()
        if (colmapping == {}) or (mapping == {}):
            mapping, colmapping = sourcegeo.block_mapping(geo, True)
        from copy import copy
        # atmosphere blocks:
        default_atm_incons = t2blockincon([1.013e5, 20.])
        if geo.atmosphere_type == 0: # single atmosphere block
            atmblk = geo.block_name(geo.layerlist[0].name, geo.atmosphere_column_name)
            if sourcegeo.atmosphere_type == 0: self[atmblk] = copy(sourceinc[0])
            elif sourcegeo.atmosphere_type == 1:
                # take average over source column atmosphere incons
                varsum = np.zeros(len(sourceinc[0].variable))
                for col in sourcegeo.columnlist:
                    blk = sourcegeo.block_name(sourcegeo.layerlist[0].name, col.name)
                    varsum+=np.array(sourceinc[blk].variable)
                self[atmblk] = t2blockincon(varsum / sourcegeo.num_columns)
            else: self[atmblk] = copy(default_atm_incons)
        elif geo.atmosphere_type == 1:
            # atmosphere block over each column
            if sourcegeo.atmosphere_type == 0:
                # broadcast single source atmosphere incons to each column
                for col in geo.columnlist:
                    blk = geo.block_name(geo.layerlist[0].name, col.name)
                    self[blk] = copy(sourceinc[0])
            elif sourcegeo.atmosphere_type == 1:
                # atmosphere over each column in both source and destination
                for col in geo.columnlist:
                    mappedcol = colmapping[col.name]
                    oldatmosblockname = sourcegeo.block_name(sourcegeo.layerlist[0].name, mappedcol)
                    blk = geo.block_name(geo.layerlist[0].name, col.name)
                    self[blk] = copy(sourceinc[oldatmosblockname])
            else: # use default
                for col in geo.columnlist:
                    blk = geo.block_name(geo.layerlist[0].name, col.name)
                    self[blk] = copy(default_atm_incons)
        # underground blocks:
        for blk in geo.block_name_list[geo.num_atmosphere_blocks:]:
            self[blk] = copy(sourceinc[mapping[blk]])
