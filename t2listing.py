"""For reading TOUGH2 listing files.

Copyright 2012 University of Auckland.

This file is part of PyTOUGH.

PyTOUGH is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

PyTOUGH is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with PyTOUGH.  If not, see <http://www.gnu.org/licenses/>."""

try:
    import numpy as np
    from numpy import float64
except ImportError: # try importing Numeric on old installs
    import Numeric as np
    from Numeric import Float64 as float64
from mulgrids import fix_blockname, valid_blockname
from fixed_format_file import fortran_float, fortran_int
import io
from sys import version_info

class listingtable(object):

    """Class for table in listing file, with values addressable by index
    (0-based) or row name, and column name: e.g. table[i] returns the
    ith row (as a dictionary), table[rowname] returns the row with the
    specified name, and table[colname] returns the column with the
    specified name."""

    def __init__(self, cols, rows, row_format = None, row_line = None, num_keys = 1,
                 allow_reverse_keys = False, header_skiplines = 0, skiplines = []):
        """The row_format parameter is a dictionary with three keys,
        'key','index' and 'values'.  These contain the positions, in
        each row of the table, of the start of the keys, index and
        data fields.  The row_line parameter is a list containing, for
        each row of the table, the number of lines before it in the
        listing file, from the start of the table.  This is needed for
        TOUGH2_MP listing files, in which the rows are not in index
        order and can also be duplicated.
        """
        self.column_name = cols
        self.row_name = rows
        self.row_format = row_format
        self.row_line = row_line
        self.num_keys = num_keys
        self.allow_reverse_keys = allow_reverse_keys
        self.header_skiplines = header_skiplines
        self.skiplines = skiplines
        self._col = dict([(c,i) for i,c in enumerate(cols)])
        self._row = dict([(r,i) for i,r in enumerate(rows)])
        self._data = np.zeros((len(rows), len(cols)), float64)

    def __repr__(self):
        return repr(self.column_name) + '\n' + repr(self._data)

    def __getitem__(self,key):
        if isinstance(key,int):
            return dict(zip(['key'] + self.column_name, [self.row_name[key]] +
                            list(self._data[key,:])))
        else:
            if key in self.column_name:
                return self._data[:,self._col[key]]
            elif key in self.row_name:
                rowindex = self._row[key]
                return dict(zip(['key'] + self.column_name,
                                [self.row_name[rowindex]] +
                                list(self._data[rowindex,:])))
            elif len(key) > 1 and self.allow_reverse_keys:
                revkey = key[::-1] # try reversed key for multi-key tables
                if revkey in self.row_name:
                    rowindex = self._row[revkey]
                    return dict(zip(['key'] + self.column_name,
                                    [self.row_name[rowindex][::-1]] +
                                    list(-self._data[rowindex,:])))
            else: return None

    def __setitem__(self, key, value):
        if isinstance(key,int): self._data[key,:] = value
        else: self._data[self._row[key],:] = value

    def get_num_columns(self):
        return len(self.column_name)
    num_columns=property(get_num_columns)

    def get_num_rows(self):
        return len(self.row_name)
    num_rows=property(get_num_rows)

    def key_from_line(self, line):
        key = [fix_blockname(line[pos: pos + 5]) for pos in self.row_format['key']]
        if len(key) == 1: return key[0]
        else: return tuple(key)

    def rows_matching(self, pattern, index = 0, match_any = False):
        """Returns rows in the table with keys matching the specified regular expression pattern
        string.
        For tables with multiple keys, pattern can be a list or tuple of regular expressions.  If
        a single string pattern is given for a multiple-key table, the pattern is matched on the
        index'th key (and any value of the other key- unless match_any is used; see below).
        If match_any is set to True, rows are returned with keys matching any of the specified
        patterns (instead of all of them).  If this option is used in conjunction with a single
        string pattern, the specified pattern is applied to all keys."""
        from re import search
        if self.num_keys == 1:
            return [self[key] for key in self.row_name if search(pattern,key)]
        else:
            if isinstance(pattern, str): pattern = [pattern]
            else: pattern = list(pattern)
            if len(pattern) < self.num_keys:
                if match_any: default = [pattern[0]]
                else: default = ['.*']
                if 0 <= index <= self.num_keys:
                    pattern = default * index + pattern + default * (self.num_keys - 1 - index)
                else: return []
            combine = [all, any][match_any]
            return [self[key] for key in self.row_name if
                    combine([search(p, n) for p, n in zip(pattern,key)])]

    def __add__(self, other):
        """Adds two listing tables together."""
        if self.column_name == other.column_name and self.row_name == other.row_name:
            from copy import copy
            result = listingtable(copy(self.column_name), copy(self.row_name), num_keys = self.num_keys,
                                  allow_reverse_keys = self.allow_reverse_keys)
            result._data = self._data + other._data
            return result
        else: raise Exception("Incompatible tables: can't be added together.")

    def __sub__(self, other):
        """Subtracts one listing table from another."""
        if self.column_name == other.column_name and self.row_name == other.row_name:
            from copy import copy
            result = listingtable(copy(self.column_name), copy(self.row_name), num_keys = self.num_keys,
                                  allow_reverse_keys = self.allow_reverse_keys)
            result._data = self._data - other._data
            return result
        else: raise Exception("Incompatible tables: can't be subtracted.")

    def is_header(self, line):
        """Returns True if all column headers in the table are present in the given line.  Used for
        detecting internal table headers in TOUGH2 listings.  Some internal header lines have different
        spacings from the top header, so a simple string comparison doesn't always work."""
        return all([col in line for col in self.column_name])

    def get_DataFrame(self):
        """Returns data in the table as a pandas DataFrame object."""
        import pandas as pd
        row_header = 'row'
        datadict = {row_header: self.row_name}
        for colname in self.column_name: datadict[colname] = self[colname]
        return pd.DataFrame(datadict, columns = [row_header] + self.column_name)
    DataFrame = property(get_DataFrame)

class t2listing(object):
    """Class for TOUGH2 listing file.  The element, connection and
       generation tables can be accessed via the element, connection
       and generation fields.  (For example, the pressure in block
       'aa100' is given by element['aa100']['Pressure'].)  It is
       possible to navigate through time in the listing by using the
       next() and prev() functions to step through, or using the
       first() and last() functions to go to the start or end, or to
       set the index, step (model time step number) or time properties
       directly."""

    def __init__(self, filename = None, skip_tables = None, encoding = 'latin-1'):
        if skip_tables is None: skip_tables = []
        self._table = {}
        self._tablenames = []
        self.filename = filename
        self.skip_tables = skip_tables
        self.encoding = encoding
        self._file = io.open(filename, 'rb', newline = None)
        self.detect_simulator()
        if self.simulator is None:
            raise Exception('Could not detect simulator type.')
        else:
            self.setup_short_types()
            self.setup_pos()
            if self.num_fulltimes > 0:
                self._index = 0
                self.setup_tables()
                self.set_table_attributes()
                self.setup_short_indices()
                self.first()
            else:
                raise Exception('No full results found in listing file.')

    def __repr__(self): return self.title
    def close(self):
        self._file.close()

    if version_info[0] < 3:
        def readline(self):
            """Reads next line and returns it as a string. In Python 2.x,
            readline() already returns a string and doesn't need decoding.
            Skipping the decoding step speeds it up."""
            return self._file.readline()
    else:
        def readline(self):
            """Reads next line, decodes it and returns it as a string."""
            return self._file.readline().decode(self.encoding)

    def get_index(self): return self._index
    def set_index(self, i):
        self._file.seek(self._fullpos[i])
        self._index = i
        if self._index < 0: self._index += self.num_fulltimes
        self.read_tables()
    index = property(get_index,set_index)

    def get_time(self): return self._time
    def set_time(self, t):
        if t < self.fulltimes[0]: self.index = 0
        elif t > self.fulltimes[-1]: self.index = -1
        else:
            dt = np.abs(self.fulltimes - t)
            self.index = np.argmin(dt)
    time = property(get_time, set_time)

    def get_num_times(self): return len(self.times)
    num_times = property(get_num_times)
    def get_num_fulltimes(self): return len(self.fulltimes)
    num_fulltimes = property(get_num_fulltimes)

    def get_step(self): return self._step
    def set_step(self, step):
        if step < self.fullsteps[0]: self.index = 0
        elif step > self.fullsteps[-1]: self.index = -1
        else:
            dstep = np.abs(self.fullsteps - step)
            self.index = np.argmin(dstep)
    step = property(get_step,set_step)

    def get_table_names(self):
        return sorted(self._table.keys())
    table_names = property(get_table_names)

    def rewind(self):
        """Rewinds to start of listing (without reading any results)"""
        self._file.seek(0)
        self._index = -1

    def first(self): self.index = 0
    def last(self): self.index = -1
    def next(self):
        """Find and read next set of results; returns false if at end of listing"""
        more = self.index < self.num_fulltimes - 1
        if more: self.index += 1
        return more
    def prev(self):
        """Find and read previous set of results; returns false if at start of listing"""
        more = self.index > 0
        if more: self.index -= 1
        return more

    def skiplines(self, number = 1):
        """Skips specified number of lines in listing file"""
        for i in range(number):  self._file.readline()

    def skipto(self, keyword = '', start = 1):
        """Skips to line starting with keyword.  keyword can be either a
        string or a list of strings, in which case it skips to a line
        starting with any of the specified strings.  Returns the
        keyword found, or false if it can't find any of them.  The
        start parameter specifies which character in the line is to be
        considered the first to search."""
        line = ''
        if isinstance(keyword, list): keywords = keyword
        else: keywords = [keyword]
        while not any([line[start:].startswith(kw) for kw in keywords]):
            line = self.readline()
            if line == '': return False
        return [kw for kw in keywords if line[start:].startswith(kw)][0]

    def skip_to_nonblank(self):
        """Skips to start of next non-blank line."""
        pos = self._file.tell()
        while not self.readline().strip(): pos = self._file.tell()
        self._file.seek(pos)

    def skip_to_blank(self):
        """Skips to start of next blank line."""
        pos = self._file.tell()
        while self.readline().strip(): pos = self._file.tell()
        self._file.seek(pos)

    def skip_over_next_blank(self):
        """Skips past next blank line."""
        while self.readline().strip(): pass

    def detect_simulator(self):
        """Detects whether the listing has been produced by AUTOUGH2,
        TOUGH2/TOUGH2_MP or TOUGH+, and sets some internal methods
        according to the simulator type."""
        self._file.seek(0)
        simulator = {'EEEEE':'AUTOUGH2','ESHORT':'AUTOUGH2','BBBBB':'AUTOUGH2',
                     '@@@@@':'TOUGH2','=====':'TOUGH+'}
        MP = self.filename.endswith('OUTPUT_DATA') and self.readline().startswith('\f') \
             and not ('@@@@@' in self.readline())
        line = ' '
        while not ('output data after' in line or 'output after' in line or line == ''):
            line = self.readline().lower()
        if line == '': self.simulator = None
        else:
            self.skip_to_nonblank()
            line = self.readline()
            if line[1:].startswith('THE TIME IS'): line = self.readline() # AUTOUGH2
            linechars = line[1:6]
            if linechars in simulator.keys():
                self.simulator = simulator[linechars]
                if self.simulator == 'TOUGH2' and MP: self.simulator += '_MP'
            else: self.simulator = None
            if self.simulator:
                # Set internal methods according to simulator type:
                simname = self.simulator.replace('+','plus')
                internal_fns = ['setup_pos','table_type','setup_table',
                                'setup_tables','read_header','read_table','next_table',
                                'read_tables','skip_to_table','read_table_line','read_title',
                                'skip_table']
                for fname in internal_fns:
                    fname_sim = fname + '_' + simname
                    # use TOUGH2 methods for TOUGH2_MP/TOUGH+ unless there are customized methods
                    # for these simulators:
                    if simname == 'TOUGH2_MP' and not hasattr(self, fname_sim):
                        fname_sim = fname_sim.replace('_MP','')
                    if simname == 'TOUGHplus' and not hasattr(self, fname_sim):
                        fname_sim = fname_sim.replace('plus','2')
                    setattr(self, fname, getattr(self, fname_sim))

    def table_type_AUTOUGH2(self, keyword):
        """Returns AUTOUGH2 table name based on the 5-character keyword read at the top of the table."""
        keytable = {'EEEEE': 'element', 'CCCCC': 'connection', 'GGGGG': 'generation'}
        if keyword in keytable: return keytable[keyword]
        else: return None

    def table_type_TOUGH2(self, headers):
        """Returns TOUGH2 table name based on a tuple of the first three column headings."""
        if headers[0:2] in [('ELEM.','INDEX'),('ELEM.','IND.')]:
            if headers[2] == 'P': return 'element'
            elif headers[2] == 'X1': return 'primary'
        else:
            keytable = {('ELEM1','ELEM2','INDEX'): 'connection',
                        ('ELEMENT','SOURCE','INDEX'): 'generation'}
            if headers in keytable: return keytable[headers]
        return None

    def table_type_TOUGHplus(self, headers):
        """Returns TOUGH+ table name based on a tuple of the first three column headings."""
        if headers[0:2] == ('ELEM','INDEX'):
            if headers[2] == 'X1': return 'primary'
            else: return 'element'
        else:
            keytable = {('ELEM1','ELEM2','INDEX'): 'connection',
                        ('ELEMENT','SOURCE','INDEX'): 'generation'}
            if headers in keytable: return keytable[headers]
        return None

    def setup_short_types(self):
        """Sets up short_types, for handling short output."""
        self.short_types = []
        if self.simulator == 'AUTOUGH2':
            startpos = self._file.tell()
            self._file.seek(0)
            done = False
            while not done:
                shortkw = self.skipto(['ESHORT','CSHORT','GSHORT'])
                if (shortkw in self.short_types) or not shortkw: done = True
                else:
                    self.short_types.append(shortkw)
                    self.skipto(shortkw)
                    self.skipto(shortkw) # to end of short table
            self._file.seek(startpos)
        # (no SHORT output for TOUGH2 or TOUGH+)

    def setup_short_indices(self):
        """Sets up short_indices (indices of main table items in short tables)."""
        self.short_indices = {}
        if self.simulator == 'AUTOUGH2':
            startpos = self._file.tell()
            shortpos = [pos for pos,short in zip(self._pos, self._short) if short]
            if len(shortpos) > 0:
                self._file.seek(shortpos[0])
                for itable,table in enumerate(self.short_types):
                    fulltable_keyword = table[0].upper()*5
                    fulltable = self.table_type(fulltable_keyword)
                    if itable > 0: self.skipto(table)
                    self.short_indices[table] = {}
                    self.skipto(table)
                    self.skip_to_blank()
                    self.skip_to_nonblank()
                    if fulltable in self._table:
                        def rowindex(line): # can get row indices from full table
                            key = self._table[fulltable].key_from_line(line)
                            return self._table[fulltable]._row[key]
                        self.skip_to_blank()
                    else:
                        # have to get row index from index in table (possibly not necessarily reliable)
                        indexpos = self.readline().index('INDEX')
                        def rowindex(line): return fortran_int(line[indexpos: indexpos + 5]) - 1
                    self.skip_to_nonblank()
                    endtable = False
                    lineindex = 0
                    while not endtable:
                        line = self.readline()
                        if line[1:].startswith(table): endtable = True
                        else:
                            index = rowindex(line)
                            self.short_indices[table][index] = lineindex
                        lineindex += 1
            self._file.seek(startpos)

    def setup_pos_AUTOUGH2(self):
        """Sets up _pos list for AUTOUGH2 listings, containing file position
        at the start of each set of results.  Also sets up the times
        and steps arrays.
        """
        self._file.seek(0)
        # set up pos,times, steps and short arrays:
        self._fullpos, self._pos, self._short = [], [], []
        fullt, fulls, t, s = [], [], [], []
        keywords = ['EEEEE']
        if len(self.short_types) > 0: keywords.append(self.short_types[0])
        endfile = False
        while not endfile:
            kwfound = self.skipto(keywords)
            if kwfound:
                self._pos.append(self._file.tell())
                self.read_header_AUTOUGH2()
                if kwfound == 'EEEEE': # full results
                    self._fullpos.append(self._pos[-1])
                    fullt.append(self.time)
                    fulls.append(self.step)
                    self._short.append(False)
                else: self._short.append(True)
                t.append(self.time)
                s.append(self.step)
                self._file.readline()
                self.skipto(kwfound) # to end of table
            else: endfile = True
        self.times = np.array(t)
        self.steps = np.array(s)
        self.fulltimes = np.array(fullt)
        self.fullsteps = np.array(fulls)

    def setup_pos_TOUGH2(self):
        """Sets up _pos list for TOUGH2 listings, containing file position at
        the start of each set of results.  Also sets up the times and steps arrays."""
        self._file.seek(0)
        # set up pos,times, steps and short arrays:
        self._fullpos,self._pos = [],[]
        t,s = [],[]
        endfile = False
        while not endfile:
            line = ' '
            while not (line.lstrip().lower().startswith('output data after') or line == ''):
                line = self.readline()
            if line != '':
                while not ('total time' in self.readline().lower()): pass
                pos = self._file.tell()
                self._pos.append(pos)
                self._fullpos.append(pos)
                self.read_header()
                t.append(self.time)
                s.append(self.step)
                self.skipto('@@@@@')
            else: endfile = True
        self.times = np.array(t)
        self.steps = np.array(s)
        self.fulltimes = np.array(t)
        self.fullsteps = np.array(s)
        self._short = [False for p in self._pos]

    def set_table_attributes(self):
        """Makes tables in self._table accessible as attributes."""
        for key,table in self._table.items(): setattr(self, key, table)

    def setup_tables_AUTOUGH2(self):
        """Sets up configuration of element, connection and generation tables."""
        tablename = 'element'
        self._file.seek(self._fullpos[0])
        while tablename:
            self.read_header()
            if tablename in self.skip_tables: self.skip_table(tablename)
            else: self.setup_table(tablename)
            tablename = self.next_table()

    def setup_tables_TOUGH2(self):
        self.read_title()
        tablename = 'element'
        self._file.seek(self._fullpos[0])
        self.read_header() # only one header at each time
        while tablename:
            if tablename in self.skip_tables: self.skip_table(tablename)
            else: self.setup_table(tablename)
            tablename = self.next_table()

    def setup_tables_TOUGHplus(self):
        self.read_title()
        tablename = 'element'
        self._file.seek(self._fullpos[0])
        self.read_header() # only one header at each time
        nelt_tables = 0 # can have multiple element tables
        while tablename:
            if tablename in self.skip_tables: self.skip_table(tablename)
            else: self.setup_table(tablename)
            tablename = self.next_table()
            if tablename == 'element':
                nelt_tables += 1
                tablename += str(nelt_tables)

    def next_table_AUTOUGH2(self):
        """Goes to start of next table at current time and returns its type,
        or None if there are no more."""
        keyword = self.readline()[1:6]
        return self.table_type(keyword)

    def next_table_TOUGH2(self):
        found = False
        while not found:
            line = '\n'
            while not ((line.strip().startswith('KCYC') and 'ITER' in line) or line == ''):
                line = self.readline()
            if line == '': return None
            else:
                pos = self._file.tell()
                if (self.num_fulltimes > 1) and (self.index < self.num_fulltimes-1):
                    if pos >= self._fullpos[self.index+1]: return None
                self.skip_to_nonblank()
                headpos = self._file.tell()
                line = self.readline().strip()
                if line == 'MASS FLOW RATES (KG/S) FROM DIFFUSION':
                    # skip over extra mass flow rate table in EOS7c listings:
                    self.skipto('@@@@@')
                else:
                    headers = tuple(line.strip().split()[0:3])
                    self._file.seek(headpos)
                    found = True
                    return self.table_type(headers)

    def next_table_TOUGHplus(self):
        if self.skipto('_____',0):
            self._file.readline()
            pos = self._file.tell()
            if (self.num_fulltimes>1) and (self.index<self.num_fulltimes-1):
                if pos >= self._fullpos[self.index+1]: return None
            headpos = self._file.tell()
            headers = tuple(self.readline().strip().split()[0:3])
            self._file.seek(headpos)
            return self.table_type(headers)
        else: return None

    def skip_to_table_AUTOUGH2(self, tablename, last_tablename, nelt_tables):
        """Skips forwards to headers of table with specified name at the current time."""
        tablechar = tablename[0].upper()
        if self._short[self._index]:
            keyword, first_tablechar = tablechar + 'SHORT', self.short_types[0][0]
        else: keyword, first_tablechar = tablechar*5, 'E'
        if tablechar != first_tablechar: self.skipto(keyword)
        self.skipto('OUTPUT')
        self.skipto(keyword)
        self.skip_to_blank()
        self.skip_to_nonblank()

    def skip_to_table_TOUGH2(self, tablename, last_tablename, nelt_tables):
        if last_tablename is None:
            self.skipto('@@@@@')
            self.skip_to_nonblank()
            tname = 'element'
        else: tname = last_tablename
        while tname != tablename:
            self.skipto('@@@@@')
            tname = self.next_table_TOUGH2()

    def skip_to_table_TOUGHplus(self, tablename, last_tablename, nelt_tables):
        if last_tablename is None:
            self.skipto('=====',0)
            self.skip_to_nonblank()
            tname = 'element'
            nelt_tables = 0
        else: tname = last_tablename
        while tname != tablename:
            if tname == 'primary': keyword='_____'
            else: keyword = '@@@@@'
            self.skipto(keyword,0)
            tname = self.next_table_TOUGHplus()
            if tname == 'element':
                nelt_tables += 1
                tname += str(nelt_tables)

    def start_of_values(self, line):
        """Returns start index of values in a table line.  Characters before
        this start index are taken to contain the key(s) and row index."""
        pt = line.find('.')
        if pt >= 2:
            nextpt = line.find('.', pt + 1)
            if nextpt < 0 : nextpt = len(line)
            s = line[pt + 1: nextpt - 1].lower()
            exponential = s.find('e') >= 0 or s.find('+') >= 0 or s.find('-') >= 0
            if exponential:
                c = line[pt - 2]
                if c in ['-',' ']: return pt - 2
                elif c.isdigit(): return pt - 1
            else:
                pos = pt - 1
                while line[pos] != ' ' and pos > 0 : pos -= 1
                while line[pos] == ' ' and pos > 0 : pos -= 1
                if pos > 0 : return pos + 1
        return None

    def key_positions(self, line, nkeys):
        """Returns detected positions of keys in the start of a table line.
        This works on the assumption that key values must have a digit
        present in their last character.  It searches backwards from
        the end of the position of INDEX (or IND.) in the header line-
        sometimes keys can overlap into this area. """
        keylength = 5
        keypos = []
        pos = len(line)-1
        while line[pos] == ' ': pos -= 1
        while line[pos] != ' ': pos -= 1
        for k in range(nkeys):
            while not line[pos].isdigit() and pos >= keylength:
                pos -= 1
            pos -= (keylength - 1)
            if valid_blockname(line[pos: pos + keylength]):
                keypos.append(pos)
                pos -= 1
            else: return None
        keypos.reverse()
        if len(keypos) != nkeys: return None
        else: return keypos

    def parse_table_header_AUTOUGH2(self):
        """Parses table header line for AUTOUGH2, returning the number of keys
        and the column names."""
        cols = []
        headline = self.readline()
        headstrs = headline.strip().split()
        indexstr = 'INDEX'
        nkeys = headstrs.index(indexstr)
        for s in headstrs[nkeys+1:]:
            if s[0] == s[0].upper(): cols.append(s)
            else: cols[-1] += ' ' + s
        return nkeys,cols

    def setup_table_AUTOUGH2(self, tablename):
        """Sets up table from AUTOUGH2 listing file."""
        keyword = tablename[0].upper()*5
        self.skiplines(3)
        nkeys, cols = self.parse_table_header_AUTOUGH2()
        self._file.readline()
        line = self.readline()
        start = self.start_of_values(line)
        rows = []
        # Double-check number of columns:
        nvalues = len([s for s in line[start:].strip().split()])
        if (len(cols) == nvalues):
            keypos = self.key_positions(line[:start], nkeys)
            if keypos:
                # determine row names:
                while line[1:6] != keyword:
                    keyval = [fix_blockname(line[kp: kp + 5]) for kp in keypos]
                    if len(keyval) > 1: keyval = tuple(keyval)
                    else: keyval = keyval[0]
                    rows.append(keyval)
                    line = self.readline()
                row_format = {'key': keypos, 'values': [start]}
                allow_rev = tablename == 'connection'
                self._table[tablename] = listingtable(cols, rows, row_format, num_keys = nkeys,
                                                      allow_reverse_keys = allow_rev)
                self._tablenames.append(tablename)
                self._file.readline()
            else: raise Exception('Error parsing '+tablename+' table keys: table not created.')
        else:
            raise Exception('Error parsing '+tablename+' table columns: table not created.')

    def parse_table_line(self, line, start):
        """Parses line of a table and returns starting indices of each column"""
        numpos = [start]
        from re import finditer,escape
        # find all decimal points:
        pts = [match.start() for match in finditer(escape('.'), line)]
        for i, pt in enumerate(pts[: -1]):
            nextpt = pts[i+1]
            pstart, pend = pt+1, nextpt-1
            spacepos = line.find(' ', pstart, pend)
            if spacepos > 0: next_start = spacepos
            else: # no space at end
                exppos = line.find('E', pstart, pend)
                if exppos > 0:
                    endpos = exppos + 3
                    next_start = endpos + 1
                else: raise Exception("Unable to parse table line:\n" + line)
            numpos.append(next_start)
        numpos.append(len(line))
        return numpos

    def parse_table_header_TOUGH2(self):
        """Parses table header line for TOUGH2, returning the number of keys and the column names."""
        cols = []
        if self.simulator in ['TOUGH2','TOUGH2_MP']: flow_headers = ['RATE']
        else: flow_headers = ['Flow','Veloc']
        headline = self.readline().strip()
        headstrs = headline.split()
        indexstrs = ['INDEX','IND.'] # for EWASG
        for s in indexstrs:
            if s in headstrs:
                nkeys = headstrs.index(s)
                break
        for s in headstrs[nkeys+1:]:
            if s in flow_headers: cols[-1] += ' ' + s
            else: cols.append(s)
        return nkeys,cols

    def is_results_line(self, line, expected_floats):
        """Detects whether the given string represents a line of results values, by seeing if
        it contains at least expected_floats floating point numbers."""
        from re import findall
        return len(findall('\.[0-9]+', line)) >= expected_floats

    def skip_to_results_line(self, expected_floats):
        """Skips to the start of the next line of results values, returning the number
        of lines skipped."""
        found, num_lines = False, 1
        while not found:
            pos = self._file.tell()
            line = self.readline().strip()
            found = self.is_results_line(line, expected_floats)
            if found: self._file.seek(pos)
            else: num_lines += 1
        return num_lines

    def table_expected_floats(self, tablename, num_columns):
        """Returns number of floating point numbers expected in each line of the table.
        For most tables this is the number of table columns, but generation tables can
        have incomplete lines in them and sometimes have as few as one value per line."""
        return 1 if tablename == 'generation' else num_columns

    def setup_table_TOUGH2(self, tablename):
        """Sets up table from TOUGH2 (or TOUGH+) listing file."""
        nkeys,cols = self.parse_table_header_TOUGH2()
        ncols = len(cols)
        expected_floats = self.table_expected_floats(tablename, ncols)
        header_skiplines = self.skip_to_results_line(expected_floats)
        line = self.readline()
        start = self.start_of_values(line)
        keypos = self.key_positions(line[:start],nkeys)
        if keypos:
            index_pos = [keypos[-1]+5, start]
            longest_line = line
            rowdict = {}
            count,index = 0,-1
            skiplines = []
            lsep = 60
            more = True
            internal_header_skiplines = None
            def count_read(count): return self.readline(), count + 1
            def is_header(line): return all([col in line for col in cols])
            def is_separator(line): return len(line)>lsep and line[1:lsep+1] == line[1]*lsep
            while more:
                keyval = [fix_blockname(line[kp:kp+5]) for kp in keypos]
                if len(keyval) > 1: keyval = tuple(keyval)
                else: keyval = keyval[0]
                indexstr = line[index_pos[0]:index_pos[1]]
                try: index = int(indexstr) - 1
                # To handle overflow (****) in index field: assume indices continue:
                except ValueError: index += 1
                # Use a dictionary to deal with duplicate row indices (TOUGH2_MP):
                rowdict[index] = (count,keyval)
                if len(line.strip()) > len(longest_line): longest_line = line
                pos = self._file.tell()
                last_count = count
                line,count = count_read(count)
                internal_header = False
                if is_header(line): internal_header = True
                elif is_separator(line): # end of table
                    more = False
                    self._file.seek(pos)
                elif not line.strip(): # blank- check next line:
                    pos = self._file.tell()
                    line,count = count_read(count)
                    stripline = line.strip()
                    if is_header(line): internal_header = True
                    elif is_separator(line) or stripline == self.title or not stripline:
                        more = False  # end of table
                        self._file.seek(pos)
                if more and internal_header:
                    if internal_header_skiplines is None:
                        internal_header_skiplines = self.skip_to_results_line(expected_floats)
                        count += internal_header_skiplines
                        line = self.readline()
                    else:
                        for i in range(internal_header_skiplines): line,count = count_read(count)
                skiplines.append(count - last_count - 1)
            indices = sorted(rowdict.keys())
            row_line = [rowdict[index][0] for index in indices]
            rows = [rowdict[index][1] for index in indices]
            numpos = self.parse_table_line(longest_line,start)
            row_format = {'key': keypos, 'index': keypos[-1] + 5, 'values': numpos}
            allow_rev = tablename == 'connection'
            self._table[tablename] = listingtable(cols, rows, row_format, row_line, num_keys = nkeys,
                                                  allow_reverse_keys = allow_rev,
                                                  header_skiplines = header_skiplines,
                                                  skiplines = skiplines)
            self._tablenames.append(tablename)
        else: raise Exception('Error parsing '+tablename+' table keys: table not created.')

    def read_header_AUTOUGH2(self):
        """Reads header info (title and time data) for one set of AUTOUGH2 listing results."""
        self.read_title()
        line = self.readline()
        istart, iend = line.find('AFTER') + 5, line.find('TIME STEPS')
        try: self._step = fortran_int(line[istart:iend])
        except ValueError: self._step = -1 # to handle overflow
        istart = iend + 10
        iend = line.find('SECONDS')
        self._time=fortran_float(line[istart:iend])
        self._file.readline()

    def read_header_TOUGH2(self):
        """Reads header info (time data) for one set of TOUGH2 listing results."""
        strs = self.readline().split()
        self._time, self._step = fortran_float(strs[0]), fortran_int(strs[1])
        marker = ['@@@@@','====='][self.simulator=='TOUGH+']
        self.skipto(marker)
        self.skip_to_nonblank()
        pos = self._file.tell()
        strs = self.readline().split()
        if len(strs) < 4: self.skip_to_nonblank() # to skip over extra lines in EOS7c listings
        else: self._file.seek(pos)

    def read_title_AUTOUGH2(self):
        """Read simulation title for AUTOUGH2 listings, from current position- in all headers."""
        self.title = self.readline().strip()

    def read_title_TOUGH2(self):
        """Reads simulation title for TOUGH2 listings, at top of file."""
        self._file.seek(0)
        line = ' '
        while not ('problem title' in line.lower() and ':' in line) or (line == ''):
            line = self.readline()
        if line == '': self.title = ''
        else:
            colonpos = line.find(':')
            if colonpos >= 0: self.title = line[colonpos + 1:].strip()
            else: self.title = ''

    def read_title_TOUGH2_MP(self):
        """Reads simulation title for TOUGH2_MP listings, at top of file."""
        self._file.seek(0)
        self._file.readline()
        self.title = self.readline().strip()

    def next_tablename(self, tablename):
        """Returns name of table after the specified one, or None if it is the last."""
        if tablename is None: return self._tablenames[0]
        i = self._tablenames.index(tablename)
        if i < len(self._tablenames)-1: return self._tablenames[i+1]
        else: return None

    def read_tables_AUTOUGH2(self):
        tablename = 'element'
        while tablename:
            self.read_header()
            if tablename in self.skip_tables: self.skip_table(tablename)
            else: self.read_table(tablename)
            tablename = self.next_table()

    def read_tables_TOUGH2(self):
        tablename = 'element'
        self.read_header() # only one header at each time
        last_tablename = None
        while tablename:
            if tablename in self.skip_tables: self.skip_table(tablename)
            elif tablename in self._table: self.read_table(tablename)
            else: # tables not present at first time step
                next_tablename = self.next_tablename(last_tablename)
                if next_tablename:
                    self.skip_to_table(next_tablename, last_tablename, 1)
            last_tablename = tablename
            tablename = self.next_table()

    def read_tables_TOUGHplus(self):
        tablename='element'
        self.read_header() # only one header at each time
        nelt_tables = 0
        while tablename:
            if tablename in self.skip_tables: self.skip_table(tablename)
            else: self.read_table(tablename)
            tablename = self.next_table()
            if tablename == 'element':
                nelt_tables += 1
                tablename += str(nelt_tables)

    def read_table_AUTOUGH2(self, tablename):
        fmt = self._table[tablename].row_format
        keyword = tablename[0].upper()*5
        self.skip_to_blank()
        self._file.readline()
        self.skip_to_blank()
        self.skip_to_nonblank()
        line = self.readline()
        row = 0
        while line[1:6] != keyword:
            self._table[tablename][row] = self.read_table_line_AUTOUGH2(line, fmt = fmt)
            row += 1
            line = self.readline()
        self._file.readline()

    def skip_table_AUTOUGH2(self, tablename):
        keyword = tablename[0].upper()*5
        self.skip_to_blank()
        line = self.readline()
        while line[1:6] != keyword: line = self.readline()
        self._file.readline()

    def read_table_line_AUTOUGH2(self, line, num_columns = None, fmt = None):
        start = fmt['values'][0]
        vals = [fortran_float(s) for s in line[start:].strip().split()]
        return vals

    def read_table_line_TOUGH2(self, line, num_columns, fmt):
        """Reads values from a line in a TOUGH2 listing, given the number of columns, and format."""
        nvals = len(fmt['values']) - 1
        return [fortran_float(line[fmt['values'][i]: fmt['values'][i+1]])
                for i in range(nvals)] + [0.0] * (num_columns - nvals)

    def read_table_TOUGH2(self, tablename):
        table = self._table[tablename]
        ncols = table.num_columns
        fmt = table.row_format
        self.skiplines(table.header_skiplines)
        for skip in table.skiplines:
            line = self.readline()
            key = table.key_from_line(line)
            table[key] = self.read_table_line_TOUGH2(line, ncols, fmt)
            self.skiplines(skip)

    def skip_table_TOUGH2(self, tablename):
        if tablename in self._table:
            table = self._table[tablename]
            self.skiplines(table.header_skiplines + table.num_rows + sum(table.skiplines))
        else:
            if self.simulator == 'TOUGH+' and tablename == 'primary': chars = '_____'
            else: chars = '@@@@@'
            self.skipto(chars)

    def history(self, selection, short = True, start_datetime = None):
        """Returns time histories for specified selection of table type, names
           (or indices) and column names.  Table type is specified as
           'e','c','g' or 'p' (upper or lower case) for element table,
           connection table, generation table or primary table
           respectively.  For TOUGH+ results, additional element
           tables may be specified as 'e1' or 'e2'.  If the short
           parameter is True, results from 'short output' (AUTOUGH2
           only) are included in the results. If a start_datetime is
           specified (a Python datetime object) then times will be
           returned as datetimes."""

        # This can obviously be done much more simply using next(), and accessing self._table,
        # but that is too slow for large listing files.  This method reads only the required data lines
        # in each table.

        def tablename_from_specification(tabletype): # expand table specification to table name:
            from string import digits
            namemap = {'e': 'element', 'c': 'connection', 'g': 'generation', 'p': 'primary'}
            type0 = tabletype[0].lower()
            if type0 in namemap:
                name = namemap[type0]
                if tabletype[-1] in digits:
                    name += tabletype[-1] # additional TOUGH+ element tables
                return name
            else: return None

        def ordered_selection(selection, tables, short_types, short_indices):
            """Given the initial history selection, returns a list of tuples of
            table name and table selections.  The tables are in the
            same order as they appear in the listing file.  Each table
            selection is a list of tuples of (table row index, column
            name, reversed, selection index) for each table, ordered
            by table row index.  This ordering means all data can be
            read sequentially to make it more efficient.  There is a
            table selection each for full and short output, to account
            for possible differences in ordering between them."""
            converted_selection = []
            for sel_index,(tspec, key, h) in enumerate(selection):
                # convert keys to indices as necessary, and expand table names:
                tablename = tablename_from_specification(tspec)
                if tablename in tables:
                    if isinstance(key, int): index, reverse = key, False
                    else:
                        index, reverse = None, False
                        if key in tables[tablename].row_name:
                            index = tables[tablename]._row[key]
                        elif len(key) > 1 and tables[tablename].allow_reverse_keys:
                            revkey = key[::-1]
                            if revkey in tables[tablename].row_name:
                                index = tables[tablename]._row[revkey]
                                reverse = True
                    if index is not None:
                        if tables[tablename].row_line:
                            index = tables[tablename].row_line[index] # find line index if needed
                        ishort = None
                        short_keyword = tspec[0].upper() + 'SHORT'
                        if short_keyword in short_types:
                            if index in short_indices[short_keyword]:
                                ishort = short_indices[short_keyword][index]
                        converted_selection.append((tablename, index, ishort, h,
                                                    reverse, sel_index))
            tables = list(set([sel[0] for sel in converted_selection]))
            # need to retain table order as in the file:
            tables = [tname for tname in
                      ['element', 'element1', 'connection',
                       'primary', 'element2', 'generation'] if tname in tables]
            tableselection, short_tableselection = [], []
            for table in tables:
                tselect = [(i, h, rev, sel_index)
                           for (tname, i, ishort, h, rev, sel_index) in
                           converted_selection if tname == table]
                tselect.sort()
                tselect_short = [(ishort,h,rev,sel_index)
                                 for (tname, i, ishort, h, rev, sel_index)
                                 in converted_selection if
                                 tname == table and ishort is not None]
                tselect_short.sort()
                tableselection.append((table, tselect, tselect_short))
            return tableselection

        old_index = self.index
        # If input just one tuple rather than a list of them:
        if isinstance(selection, tuple): selection = [selection] 
        tableselection = ordered_selection(selection, self._table,
                                           self.short_types, self.short_indices)
        if len(tableselection) == 0:
            return None # no valid specifications
        hist = [[] for s in selection]
        self.rewind()

        for ipos, pos in enumerate(self._pos):
            self._file.seek(pos)
            self._index = ipos
            is_short = self._short[ipos]
            if not (is_short and not short):
                last_tname = None
                nelt_tables = -1
                for (tname, tselect, tselect_short) in tableselection:
                    if is_short: tablename = tname[0].upper() + 'SHORT'
                    else: tablename = tname
                    if not (is_short and not (tablename in self.short_types)):
                        self.skip_to_table(tname, last_tname, nelt_tables)
                        if tname.startswith('element'): nelt_tables += 1
                        ncols = self._table[tname].num_columns
                        expected_floats = self.table_expected_floats(tname, ncols)
                        self.skip_to_results_line(expected_floats)
                        fmt = self._table[tname].row_format
                        index = 0
                        line = self.readline()
                        ts = tselect_short if is_short else tselect
                        for (lineindex, colname, reverse, sel_index) in ts:
                            if lineindex is not None:
                                if lineindex > index:
                                    for k in range(lineindex - index - 1):
                                        self._file.readline()
                                    line = self.readline()
                                index = lineindex
                                vals = self.read_table_line(line, ncols, fmt)
                                valindex = self._table[tname]._col[colname]
                                sgn = [1.,-1.][reverse]
                                hist[sel_index].append(sgn*vals[valindex])
                    last_tname = tname

        self._index = old_index
        short_times, all_times = self.times, self.fulltimes
        if start_datetime is not None:
            from datetime import datetime, timedelta
            def datetime_array(t): return np.array([start_datetime + timedelta(0, s) for s in t])
            short_times, all_times = datetime_array(short_times), datetime_array(all_times)
        result = [([short_times, all_times][len(h) == self.num_fulltimes],
                   np.array(h)) for sel_index, h in enumerate(hist)]
        if len(result) == 1: result = result[0]
        return result

    def get_reductions(self):
        """Returns a list of time step indices at which the time step is
        reduced, and the blocks at which the maximum residual occurred
        prior to the reduction."""
        self.rewind()
        line, lastline = '', ''
        keyword = "+++++++++ REDUCE TIME STEP"
        keyend = len(keyword)+1
        rl = []
        finished = False
        while not finished:
            while line[1:keyend] != keyword:
                lastline=line
                line = self.readline()
                if not line:
                    finished = True; break
            if not finished:
                lowerlastline = lastline.lower()
                eltindex = lowerlastline.find('element')
                if eltindex > 0:
                    if lowerlastline.find('eos cannot find parameters') >= 0: space = 9
                    else: space = 8
                    blockname = fix_blockname(lastline[eltindex+space: eltindex+space+5])
                    brackindex, comindex = line.find('('),line.find(',')
                    timestep = fortran_int(line[brackindex + 1: comindex])
                    rl.append((timestep, blockname))
                lastline = line
                line = self.readline()
                if not line: finished = True
        return rl
    reductions=property(get_reductions)

    def get_difference(self, indexa = None, indexb = None):
        """Returns dictionary of maximum differences, and locations of
        difference, of all element table properties between two sets
        of results.  If both indexa and indexb are provided, the
        result is the difference between these two result indices.  If
        only one index is given, the result is the difference between
        the given index and the one before that.  If neither are
        given, the result is the difference between the last and
        penultimate sets of results.
        """
        from copy import deepcopy
        tablename = 'element'
        if indexa is None: self.last()
        else: self.set_index(indexa)
        results2 = deepcopy(self._table[tablename])
        if indexb is None: self.prev()
        else: self.set_index(indexb)
        results1 = self._table[tablename]
        cvg = {}
        for name in results1.column_name:
            iblk = np.argmax(abs(results2[name] - results1[name]))
            blkname = results1.row_name[iblk]
            diff = results2[name][iblk] - results1[name][iblk]
            cvg[name] = (diff,blkname)
        return cvg
    convergence = property(get_difference)

    def get_vtk_data(self, geo, grid = None, flows = False, flux_matrix = None,
                     geo_matches = True, blockmap = {}):
        """Returns dictionary of VTK data arrays from listing file at current
        time.  If flows is True, average flux vectors are also
        calculated from connection data at the block centres.
        """
        from vtk import vtkFloatArray
        natm = geo.num_atmosphere_blocks
        nele = geo.num_underground_blocks
        arrays = {'Block': {}, 'Node': {}}
        elt_tablenames = [key for key in self._table.keys() if key.startswith('element')]
        for tablename in elt_tablenames:
            for name in self._table[tablename].column_name: arrays['Block'][name] = vtkFloatArray()
        flownames = []
        def is_flowname(name):
            name = name.lower()
            return name.startswith('flo') or name.endswith('flo') or \
                name.endswith('flow') or name.endswith('veloc')
        if flows:
            if flux_matrix is None: flux_matrix = grid.flux_matrix(geo, blockmap)
            flownames = [name for name in self.connection.column_name if is_flowname(name)]
            for name in flownames: arrays['Block'][name] = vtkFloatArray()
        array_length = {'Block': nele, 'Node': 0}
        array_data = {'Block': {}, 'Node': {}}
        def mname(blk): return blockmap[blk] if blk in blockmap else blk
        for array_type, array_dict in arrays.items():
            for name, array in array_dict.items():
                if name in flownames:
                    array.SetName(name + '/area')
                    array.SetNumberOfComponents(3)
                    array.SetNumberOfTuples(array_length[array_type])
                    array_data[array_type][name] = flux_matrix * self.connection[name]
                else:
                    array.SetName(name)
                    array.SetNumberOfComponents(1)
                    array.SetNumberOfValues(array_length[array_type])
                    for tablename in elt_tablenames:
                        if geo_matches:
                            array_data[array_type][name] = self._table[tablename][name][natm:] # faster
                        else:  # more flexible
                            array_data[array_type][name] = \
                                np.array([self._table[tablename][mname(blk)][name]
                                          for blk in geo.block_name_list[natm:]])
        for array_type, data_dict in array_data.items():
            for name, data in data_dict.items():
                if name in flownames:
                    for iblk in range(nele):
                        arrays[array_type][name].SetTuple3(iblk, data[3*iblk], data[3*iblk+1], data[3*iblk+2])
                else:
                    for iblk in range(nele):
                        arrays[array_type][name].SetValue(iblk, data[iblk])
        return arrays

    def write_vtk(self, geo, filename, grid = None, indices = None, flows = False,
                  wells = False, start_time = 0.0, time_unit = 's',
                  flux_matrix = None, blockmap = {}, surface_snap = 0.1):
        """Writes VTK files for a vtkUnstructuredGrid object corresponding to
        the grid in 3D with the listing data, with the specified
        filename, for visualisation with VTK.  A t2grid can optionally
        be specified, to include rock type data as well.  A list of
        the required time indices can optionally be specified.  If a
        grid is specified, flows is True, and connection data are
        present in the listing file, approximate average flux vectors
        are also calculated at the block centres from the connection
        data.
        """
        from vtk import vtkXMLUnstructuredGridWriter
        from os.path import splitext
        base, ext = splitext(filename)
        if wells: geo.write_well_vtk()
        geo_matches = geo.block_name_list == self.element.row_name
        doflows = False
        if flows and (self.connection is not None):
            if grid is None:
                raise Exception("t2listing.write_vtk(): if flows == True, " +
                                " a t2grid object must be specified.")
            else:
                if geo_matches or len(blockmap) > 0: doflows = True
                else:
                    raise Exception("t2listing.write_vtk(): if flows == True, " +
                                    "block names in the listing file and geometry must match, or " +
                                    "a block mapping must be specified.")            
        arrays = geo.get_vtk_data(blockmap)
        if grid is not None:
            grid_arrays = grid.get_vtk_data(geo, blockmap)
            for array_type, array_dict in arrays.items():
                array_dict.update(grid_arrays[array_type])
        if doflows and flux_matrix is None: flux_matrix = grid.flux_matrix(geo, blockmap)
        import xml.dom.minidom
        pvd = xml.dom.minidom.Document()
        vtkfile = pvd.createElement('VTKFile')
        vtkfile.setAttribute('type','Collection')
        pvd.appendChild(vtkfile)
        collection = pvd.createElement('Collection')
        initial_index = self.index
        if indices is None: indices = range(self.num_fulltimes)
        timescales = {'s': 1.0, 'h': 3600., 'd': 3600.*24, 'y': 3600. * 24 * 365.25}
        if time_unit in timescales: timescale = timescales[time_unit]
        else: timescale = 1.0
        writer = vtkXMLUnstructuredGridWriter()
        for i in indices:
            self.index = i
            t = start_time + self.time / timescale
            filename_time = base + '_' + str(i) + '.vtu'
            results_arrays = self.get_vtk_data(geo, grid, flows = doflows, flux_matrix = flux_matrix,
                                               geo_matches = geo_matches, blockmap = blockmap)
            for array_type,array_dict in arrays.items():
                array_dict.update(results_arrays[array_type])
            vtu = geo.get_vtk_grid(arrays, surface_snap)
            writer.SetFileName(filename_time)
            if hasattr(writer, 'SetInput'): writer.SetInput(vtu)
            elif hasattr(writer, 'SetInputData'): writer.SetInputData(vtu)
            writer.Write()
            dataset = pvd.createElement('DataSet')
            dataset.setAttribute('timestep',str(t))
            dataset.setAttribute('file',filename_time)
            collection.appendChild(dataset)
        vtkfile.appendChild(collection)
        pvdfile = open(base+'.pvd','w')
        pvdfile.write(pvd.toprettyxml())
        pvdfile.close()
        self.index = initial_index

    def add_side_recharge(self, geo, dat):
        """Adds side recharge generators to a TOUGH2 data object for a production run,
        calculated according to the final results in the listing.  These generators represent side
        inflows due to pressure changes in the blocks on the model's horizontal boundaries.
        Recharge generators are given the names of their blocks- any existing generators with the same
        names will be overwritten."""
        from IAPWS97 import cowat,visc
        from geometry import line_projection
        from t2data import t2generator
        initial_index = self.index
        keyword = {'AUTOUGH2':{'P': 'Pressure', 'T': 'Temperature'},
                   'TOUGH2': {'P': 'P', 'T': 'T'},
                   'TOUGH2_MP': {'P': 'P', 'T': 'T'},
                   'TOUGH+' : {'P': 'Pressure', 'T': 'Temperature'}}
        self.last()
        bdy_nodes = geo.boundary_nodes
        for blk in dat.grid.blocklist[geo.num_atmosphere_blocks:]:
            colname = geo.column_name(blk.name)
            if colname in geo.column:
                col = geo.column[colname]
                if col.num_neighbours < col.num_nodes:
                    k = 0.5 * np.sum(blk.rocktype.permeability[0:2])
                    p0 = self.element[blk.name][keyword[self.simulator]['P']]
                    t0 = self.element[blk.name][keyword[self.simulator]['T']]
                    rho, u = cowat(t0, p0)
                    h = u + p0 / rho
                    xnu = visc(rho,t0) / rho
                    coef = 0.
                    for iface in range(col.num_nodes):
                        facenode = [col.node[i] for i in [iface,(iface + 1)%col.num_nodes]]
                        if all([node in bdy_nodes for node in facenode]):
                            side_length = np.linalg.norm(facenode[1].pos - facenode[0].pos)
                            height = blk.volume / col.area
                            area = side_length * height
                            facepos = line_projection(col.centre, [node.pos for node in facenode])
                            dist = np.linalg.norm(col.centre - facepos)
                            coef += 0.5 * area * k / (xnu * dist) # recharge coefficient
                    gen_name = blk.name
                    dat.add_generator(t2generator(gen_name, blk.name, type = 'RECH',
                                                  gx = coef, ex = h, hg = p0))
        self.index = initial_index

class t2historyfile(object):
    """Class for TOUGH2 FOFT, COFT and GOFT files (history of element,
    connection and generator variables)."""

    def __init__(self, filename = None):
        self.filename = filename
        self.empty()
        if self.filename: self.read(filename)

    def empty(self):
        self.simulator = None
        self.type = None
        self._data = []
        self._row = None
        self.times = []
        self.keys = []
        self.key_name = []
        self.times = []
        self._keyrows = {}
        self.column_name = []
        self.row_name = []

    def get_num_keys(self): return len(self.keys)
    num_keys = property(get_num_keys)
    def get_num_times(self): return len(self.times)
    num_times = property(get_num_times)
    def get_num_rows(self): return len(self._data)
    num_rows = property(get_num_rows)
    def get_num_columns(self): return len(self.column_name)
    num_columns = property(get_num_columns)
    def get_keytype(self):
        return {'FOFT': 'block', 'COFT': 'connection',
                'GOFT': 'generator', None: 'key'}[self.type]
    keytype = property(get_keytype)

    def __repr__(self):
        if self._nkeys > 0: nkeystr = str(self.num_keys)
        else: nkeystr = 'summed'
        return 'History results for ' + nkeystr + ' ' + self.keytype + \
            's at ' + str(self.num_times) + ' times'

    def __getitem__(self, key):
        """Returns results for a given key value, in the form of a dictionary
        with keys corresponding to the column names.  Each value is an
        array of history results for that key and column name- unless
        the key also has the time appended as its third element, in
        which case the dictionary values are just single floating
        point values for each column.
        """
        if self.num_rows > 0:
            if self._nkeys > 0:
                if not isinstance(key, tuple): key = (key,)
                if key in self.keys:
                    keydata = self._data[self._keyrows[key]]
                    return dict([(colname, keydata[:,icol+1]) for
                                 icol, colname in enumerate(self.column_name)])
                elif key in self._row:
                    row = self._data[self._row[key]]
                    return dict([(colname,row[icol+1]) for
                                 icol, colname in enumerate(self.column_name)])
                else: return None
            else: # no keys (e.g. TOUGH+ COFT/GOFT)
                try:
                    icol = self.column_name.index(key)
                    return self._data[:, icol + 1]
                except ValueError: return None
        else: return None

    def read(self, filename):
        """Reads contents of file(s) and stores in appropriate data structures."""
        self._rowindex = 0
        from glob import glob
        files = glob(filename)
        configured = False
        for i, fname in enumerate(files):
            self._file = open(fname,'rU')
            header = self._file.readline()
            if header:
                if not configured:
                    self.detect_simulator(header)
                    self.setup_headers(fname,header)
                self.read_data(configured)
                if self.num_columns>0: configured=True
            self._file.close()
        self.finalize_data()

    def detect_simulator(self, header):
        """Detects simulator (TOUGH2 or TOUGH2_MP) from header line."""
        if 'OFT' in header: self.simulator='TOUGH2_MP'
        elif header.startswith('Time ['): self.simulator='TOUGH+'
        else: self.simulator='TOUGH2'
        internal_fns = ['setup_headers','read_data']
        simname = self.simulator.replace('+', 'plus')
        for fname in internal_fns:
            fname_sim = fname + '_' + simname
            setattr(self, fname, getattr(self, fname_sim))

    def setup_headers_TOUGH2(self, filename, header):
        """Sets up keys and column headings from given filename and header
        line, for TOUGH2 output."""
        if filename.endswith('OFT') and len(filename) >= 4:
            self.type = filename[-4:].strip()
        else: self.type = None
        self.time_index = 1
        self._nkeys = 1
        items = header.strip().split(',')
        if items[-1] == '': del items[-1] # often an extra comma on the end of the lines
        items = items[3:]
        int_index = 0
        for i, item in enumerate(items):
            try:
                int_item = int(item)
                int_index = i
                break
            except: pass
        if int_index == 0: ncols = len(items)
        else: ncols = int_index
        self.column_name=range(ncols)

    def setup_headers_TOUGH2_MP(self, filename, header):
        """Sets up keys and column headings from given filename and header
        line, for TOUGH2_MP output."""
        headers = header.strip().split()
        self.type = headers[0] # FOFT, COFT or GOFT
        time_header = [h for h in headers if h.lower().startswith('time')][0]
        self.time_index = headers.index(time_header)
        if self.type == 'FOFT':
            self.key_index = self.time_index - 1
            self._nkeys = 1
        else: # COFT or GOFT
            self.key_index = self.time_index + 1
            self._nkeys = 2
        self.key_name = headers[self.key_index: self.key_index+self._nkeys]
        prepend_titles, append_titles = ['GAS','GENERATION'], ['flow']
        startcol = self._nkeys + 2
        cols = []
        i = startcol
        while i <= len(headers) - 1:
            title = headers[i]
            if title in prepend_titles:
                cols.append(title + ' ' + headers[i+1])
                i += 1
            elif title in append_titles:
                cols[-1] += ' ' + title
            else: cols.append(title)
            i += 1
        self.column_name = cols
        self.col_start = [header.index(colname) for colname in self.column_name]
        self.key_start = [header.index(key) for key in self.key_name]
        self.time_pos = [header.index(time_header)]
        if self.type == 'FOFT':
            self.key_start.append(self.time_pos[0])
            self.time_pos.append(self.col_start[0])
        else:
            self.key_start.append(self.col_start[0])
            self.time_pos.append(self.key_start[0])

    def setup_headers_TOUGHplus(self, filename, header):
        """Sets up keys and column headings from given filename and header
        line, for TOUGH+ output."""
        headers = header.strip().split('-')
        if filename.endswith('OFT') and len(filename) >= 4:
            self.type = filename[-4:].strip()
        elif '_Time_Series' in filename:
            filetype = {'Elem_Time_Series': 'FOFT', 'Conx_Time_Series': 'COFT',
                        'SS_Time_Series': 'GOFT'}
            for key in filetype.keys():
                if key in filename:
                    self.type = filetype[key]
                    break
        else: self.type = None
        if self.type == 'FOFT':
            self.time_index = 1
            self._nkeys = 1
        else:
            self.time_index = 0
            self._nkeys = 0
        cols = headers[self._nkeys + 1:]
        from re import sub
        cols = [sub('\[.*\]','',col).strip() for col in cols] # remove units
        self.column_name = cols

    def read_data_TOUGH2(self, configured):
        """Reads in the data, for TOUGH2 output."""
        self._file.seek(0)
        lines = self._file.readlines()
        for line in lines:
            items = line.strip().split(',')
            if items[-1] == '': del items[-1]
            time_index = fortran_int(items.pop(0))
            time = float(items.pop(0))
            self.times.append(time)
            nc1 = self.num_columns + 1
            nsets = len(items) // nc1
            for i in range(nsets):
                setvals = items[i*nc1: (i+1)*nc1]
                key = (fortran_int(setvals[0]),)
                vals = [fortran_float(val) for val in setvals[1:]]
                self.row_name.append(key+(time,))
                if not key in self.keys:
                    self._keyrows[key] = []
                    self.keys.append(key)
                self._keyrows[key].append(self._rowindex)
                self._data.append([time]+vals)
                self._rowindex += 1

    def read_data_TOUGH2_MP(self, configured):
        """Reads in the data, for TOUGH2_MP output."""
        def get_key(line):
            return tuple([fix_blockname(line[self.key_start[i]:
                                             self.key_start[i+1]].rstrip()) for
                          i in range(self._nkeys)])
        def get_time(line):
            return fortran_float(line[self.time_pos[0]:self.time_pos[1]])
        def get_vals(line):
            start = self.col_start + [len(line)] # allow for lines of different lengths
            return [fortran_float(line[start[i]: start[i+1]]) for
                    i in range(self.num_columns)]
        lines = self._file.readlines()
        first_key = None
        from copy import copy
        otherfile_keys = copy(self.keys)
        for line in lines:
            if line.strip():
                time = get_time(line)
                key = get_key(line)
                if not first_key: first_key = key
                vals = get_vals(line)
                rowname = key + (time,)
                if not (key in otherfile_keys):
                    if key in self._keyrows:
                        keyrows = self._keyrows[key]
                        keytimes = [self._data[irow][0] for irow in self._keyrows[key]]
                        keyrownames = [self.row_name[irow] for irow in self._keyrows[key]]
                        newtime = not (time in keytimes)
                        newrowname = not (rowname in keyrownames)
                    else: newtime, newrowname = True, True
                    if (not configured) and (key == first_key) and newtime:
                        self.times.append(time)
                    if newrowname:
                        self.row_name.append(rowname)
                        if not key in self.keys:
                            self._keyrows[key] = []
                            self.keys.append(key)
                        self._keyrows[key].append(self._rowindex)
                        self._data.append([time] + vals)
                        self._rowindex += 1

    def read_data_TOUGHplus(self, configured):
        """Reads in the data, for TOUGH+ output."""
        lines = self._file.readlines()
        if self.type == 'FOFT':
            for line in lines:
                items = line.strip().split(',')
                time = fortran_float(items[1])
                self.times.append(time)
                nc1 = self.num_columns + 1
                nsets = (len(items)-2) // nc1
                for i in range(nsets):
                    setvals = items[2+i*nc1: 2+(i+1)*nc1]
                    key = (fortran_int(setvals[0]),)
                    if key[0] > 0:
                        vals = [fortran_float(val) for val in setvals[1:]]
                        self.row_name.append(key + (time,))
                        if not key in self.keys:
                            self._keyrows[key] = []
                            self.keys.append(key)
                        self._keyrows[key].append(self._rowindex)
                        self._data.append([time] + vals)
                        self._rowindex += 1
        else:
            for line in lines:
                vals = [fortran_float(val) for val in line.strip().split()]
                time = vals.pop(0)
                self.times.append(time)
                self._data.append([time] + vals)

    def finalize_data(self):
        self._data = np.array(self._data, float64)
        self.times = np.array(self.times, float64)
        self._row = dict([(r,i) for i,r in enumerate(self.row_name)])


class toughreact_tecplot(object):
    """Class for TOUGHREACT Tecplot output files. These work similarly to
    t2listing objects, but have just a single element table at each
    time. It is possible to navigate through time by using the next()
    and prev() functions to step through, or using the first() and
    last() functions to go to the start or end, or to set the index,
    step (model time step number) or time properties directly.  When
    reading a toughreact_tecplot object from file it is necessary also
    to specify the block names (as these are not stored in the Tecplot
    file).
    """
    def __init__(self, filename, blocks):
        self.filename = filename
        self._file = open(filename, 'rU')
        self.setup_pos()
        self.setup_table(blocks)
        if self.num_times > 0:
            self._index = 0
            self.first()
        else: raise Exception('No results found in TOUGHREACT Tecplot file ' + filename)

    def __repr__(self): return "TOUGHREACT results for " + str(self.element.num_rows) + " blocks"

    def get_index(self): return self._index
    def set_index(self, i):
        self._file.seek(self._pos[i])
        self._index = i
        if self._index < 0: self._index += self.num_times
        self.read_table()
    index = property(get_index, set_index)

    def get_time(self): return self.times[self._index]
    def set_time(self, t):
        if t < self.times[0]: self.index=0
        elif t > self.times[-1]: self.index = -1
        else:
            dt = np.abs(self.times - t)
            self.index = np.argmin(dt)
    time = property(get_time, set_time)

    def get_num_times(self): return len(self.times)
    num_times = property(get_num_times)

    def rewind(self):
        """Rewinds to start (without reading any results)"""
        self._file.seek(0)
        self._index = -1

    def first(self): self.index = 0
    def last(self): self.index = -1
    def next(self):
        """Find and read next set of results; returns false if at end of file"""
        more = self.index < self.num_times - 1
        if more: self.index += 1
        return more
    def prev(self):
        """Find and read previous set of results; returns false if at start of file"""
        more = self.index > 0
        if more: self.index -= 1
        return more

    def skipto(self, keyword, count = False):
        """Advances file to next line starting with specified keyword, returning the line."""
        line = ''
        num_lines = 0
        while not line.startswith(keyword):
            line = self._file.readline()
            num_lines += 1
            if line == '': return None, num_lines
        if count: return line, num_lines
        else: return line

    def find_next_time(self):
        """Advances to set of results at next time, and returns the time value."""
        line, num_lines = self.skipto('ZONE', count = True)
        if line is None: return None, num_lines
        quotepos = line.find('"')
        if quotepos >= 0:
            spacepos = line.find(' ', quotepos)
            if spacepos >= 0: return float(line[quotepos+1:spacepos]), num_lines
            else: return None, num_lines
        else: return None, num_lines

    def setup_pos(self):
        """Sets up _pos list for TOUGHREACT Tecplot files, containing file position at the start
        of each set of results. Also sets up the times array."""
        self._file.seek(0)
        self._pos = []
        t = []
        endfile = False
        num_blocks = None
        while not endfile:
            time, num_lines = self.find_next_time()
            if time is not None:
                self._pos.append(self._file.tell())
                t.append(time)
            else: endfile = True
            if num_lines: num_blocks = num_lines - 1
        self.times = np.array(t)
        self._num_blocks = num_blocks

    def setup_table(self, blocks):
        """Sets up element table structure. Table columns are read from the VARIABLES line in the file.
        Table rows are block names, supplied as a list of strings, or taken from a mulgrid or t2data
        object."""
        import mulgrids as mg
        import t2grids as t2g
        if isinstance(blocks, mg.mulgrid): blocks = blocks.block_name_list
        elif isinstance(blocks, t2g.t2grid): blocks = [blk.name for blk in  blocks.blocklist]
        if len(blocks) != self._num_blocks:
            raise Exception("Specified block name list is the wrong length for " +
                            "TOUGHREACT Tecplot file "+ self.filename)
        self._file.seek(0)
        line = self.skipto('VARIABLES')
        if line is not None:
            eqpos = line.find('=')
            cols = [col.strip() for col in line[eqpos+1:].strip().split(',')[:-1]]
            self.element = listingtable(cols, blocks, num_keys = 1)
        else:
            raise Exception("Could not find variable definitions " +
                            "for TOUGHREACT Tecplot file " + self.filename)

    def read_table_line(self, line):
        """Parses given string and returns an array of float values."""
        return np.fromstring(line, sep = ' ')

    def read_table(self):
        """Reads table data at the current time."""
        for i,blk in enumerate(self.element.row_name):
            line = self._file.readline().strip()
            self.element[i] = self.read_table_line(line)

    def history(self, selection):
        """Returns time histories for specified selection of block names (or
        index) and column names."""

        def ordered_selection(selection):
            osel = []
            for sel_index,(key,h) in enumerate(selection):  # convert keys to indices as necessary
                if isinstance(key, int): index = key
                else:
                    if key in self.element.row_name: index = self.element._row[key]
                    else: index = None
                if index is not None: osel.append((index,h,sel_index))
            osel.sort()
            return osel

        old_index = self.index
        if isinstance(selection, tuple): selection = [selection]
        osel = ordered_selection(selection)
        if len(osel) == 0: return None # no valid specifications
        hist = [[] for s in selection]
        self.rewind()

        for ipos,pos in enumerate(self._pos):
            self._file.seek(pos)
            self._index = ipos
            index = 0
            line = self._file.readline()
            for (lineindex, colname, sel_index) in osel:
                if lineindex is not None:
                    for k in range(lineindex-index): line = self._file.readline()
                    index = lineindex
                    vals = self.read_table_line(line)
                    valindex = self.element._col[colname]
                    hist[sel_index].append(vals[valindex])
        self._index = old_index
        result = [(self.times,np.array(h)) for sel_index,h in enumerate(hist)]
        if len(result) == 1: result = result[0]
        return result

    def get_vtk_data(self, geo, grid = None, geo_matches = True, blockmap = {}):
        """Returns dictionary of VTK data arrays from Tecplot file at current time."""
        from vtk import vtkFloatArray
        natm = geo.num_atmosphere_blocks
        nele = geo.num_underground_blocks
        arrays = {'Block':{}, 'Node':{}}
        for name in self.element.column_name: arrays['Block'][name] = vtkFloatArray()
        array_length = {'Block':nele, 'Node':0}
        array_data = {'Block':{}, 'Node':{}}
        def mname(blk): return blockmap[blk] if blk in blockmap else blk
        for array_type,array_dict in arrays.items():
            for name,array in array_dict.items():
                array.SetName(name)
                array.SetNumberOfComponents(1)
                array.SetNumberOfValues(array_length[array_type])
                if geo_matches: array_data[array_type][name] = self.element[name][natm:] # faster
                else:  # more flexible
                    array_data[array_type][name] = np.array([self.element[mname(blk)][name]
                                                             for blk in geo.block_name_list[natm:]])

        for array_type,data_dict in array_data.items():
            for name,data in data_dict.items():
                for iblk in range(nele):
                    arrays[array_type][name].SetValue(iblk, data[iblk])
        return arrays

    def write_vtk(self, geo, filename, grid = None, indices = None, start_time = 0.0,
                  time_unit = 's', blockmap = {}, surface_snap = 0.1):
        """Writes VTK files for a vtkUnstructuredGrid object corresponding to
        the grid in 3D with the Tecplot data, with the specified
        filename, for visualisation with VTK.  A t2grid can optionally
        be specified, to include rock type data as well.  A list of
        the required time indices can optionally be specified.
        """
        from vtk import vtkXMLUnstructuredGridWriter
        from os.path import splitext
        base, ext = splitext(filename)
        geo_matches = geo.block_name_list == self.element.row_name
        arrays = geo.get_vtk_data(blockmap)
        if grid is not None:
            grid_arrays = grid.get_vtk_data(geo, blockmap)
            for array_type,array_dict in arrays.items():
                array_dict.update(grid_arrays[array_type])
        import xml.dom.minidom
        pvd = xml.dom.minidom.Document()
        vtkfile = pvd.createElement('VTKFile')
        vtkfile.setAttribute('type','Collection')
        pvd.appendChild(vtkfile)
        collection = pvd.createElement('Collection')
        initial_index = self.index
        if indices is None: indices = range(self.num_times)
        yr = 3600. * 24 * 365.25
        timescales = {'s': 1.0, 'h': 3600., 'd': 3600. * 24, 'y': yr}
        if time_unit in timescales:
            timescale = timescales[time_unit] / yr # assumes Tecplot times are in years
        else: timescale = 1.0
        writer = vtkXMLUnstructuredGridWriter()
        for i in indices:
            self.index = i
            t = start_time + self.time / timescale
            filename_time = base + '_' + str(i) + '.vtu'
            results_arrays = self.get_vtk_data(geo, grid, geo_matches = geo_matches,
                                               blockmap = blockmap)
            for array_type,array_dict in arrays.items():
                array_dict.update(results_arrays[array_type])
            vtu = geo.get_vtk_grid(arrays, surface_snap)
            writer.SetFileName(filename_time)
            if hasattr(writer, 'SetInput'): writer.SetInput(vtu)
            elif hasattr(writer, 'SetInputData'): writer.SetInputData(vtu)
            writer.Write()
            dataset = pvd.createElement('DataSet')
            dataset.setAttribute('timestep', str(t))
            dataset.setAttribute('file', filename_time)
            collection.appendChild(dataset)
        vtkfile.appendChild(collection)
        pvdfile = open(base+'.pvd', 'w')
        pvdfile.write(pvd.toprettyxml())
        pvdfile.close()
        self.index = initial_index
