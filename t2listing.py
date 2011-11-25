"""For reading TOUGH2 listing files."""

"""
Copyright 2011 University of Auckland.

This file is part of PyTOUGH.

PyTOUGH is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

PyTOUGH is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with PyTOUGH.  If not, see <http://www.gnu.org/licenses/>."""

import string
try:
    import numpy as np
    from numpy import float64
except ImportError: # try importing Numeric on old installs
    import Numeric as np
    from Numeric import Float64 as float64
from mulgrids import fix_blockname, valid_blockname, fortran_float

class listingtable(object):
    """Class for table in listing file, with values addressable by index (0-based) or row name, and column name:
    e.g. table[i] returns the ith row (as a dictionary), table[rowname] returns the row with the specified name,
    and table[colname] returns the column with the specified name."""
    def __init__(self,cols,rows,row_format=None,row_line=None):
        """The row_format parameter is a dictionary with three keys, 'key','index' and 'values'.  These contain the positions,
        in each row of the table, of the start of the keys, index and data fields.  The row_line parameter is a list containing,
        for each row of the table, the number of lines before it in the listing file, from the start of the table.  This is
        needed for TOUGH2_MP listing files, in which the rows are not in index order and can also be duplicated."""
        self.column_name=cols
        self.row_name=rows
        self.row_format=row_format
        self.row_line=row_line
        self._col=dict([(c,i) for i,c in enumerate(cols)])
        self._row=dict([(r,i) for i,r in enumerate(rows)])
        self._data=np.zeros((len(rows),len(cols)),float64)
    def __repr__(self): return repr(self.column_name)+'\n'+repr(self._data)
    def __getitem__(self,key):
        if isinstance(key,int): return dict(zip(self.column_name,self._data[key,:]))
        else:
            if key in self.column_name: return self._data[:,self._col[key]]
            elif key in self.row_name: return dict(zip(self.column_name,self._data[self._row[key],:]))
            else: return None 
    def __setitem__(self,key,value):
        if isinstance(key,int): self._data[key,:]=value
        else: self._data[self._row[key],:]=value
    def get_num_columns(self):
        return len(self.column_name)
    num_columns=property(get_num_columns)
    def get_num_rows(self):
        return len(self.row_name)
    num_rows=property(get_num_rows)
    def key_from_line(self,line):
        key=[fix_blockname(line[pos:pos+5]) for pos in self.row_format['key']]
        if len(key)==1: return key[0]
        else: return tuple(key)

class t2listing(file):
    """Class for TOUGH2 listing file.  The element, connection and generation tables can be accessed
       via the element, connection and generation fields.  (For example, the pressure in block 'aa100' is
       given by element['aa100']['Pressure'].)  It is possible to navigate through time in the listing by 
       using the next() and prev() functions to step through, or using the first() and last() functions to go to 
       the start or end, or to set the index, step (model time step number) or time properties directly."""
    def __init__(self,filename=None):
        super(t2listing,self).__init__(filename)
        self.detect_simulator()
        if self.simulator==None: print 'Could not detect simulator type.'
        else:
            self.setup_short()
            self.setup_pos()
            if self.num_fulltimes>0:
                self.setup_tables()
                self.first()
            else: print 'No full results found in listing file.'

    def __repr__(self): return self.title

    def get_index(self): return self._index
    def set_index(self,i):
        self.seek(self._fullpos[i])
        self._index=i
        if self._index<0: self._index+=self.num_fulltimes
        self.read_tables()
    index=property(get_index,set_index)

    def get_time(self): return self._time
    def set_time(self,t):
        if t<self.fulltimes[0]: self.index=0
        elif t>self.fulltimes[-1]: self.index=-1
        else:
            i=[j for j,tj in enumerate(self.fulltimes) if tj>=t]
            if len(i)>0: self.index=i[0]
    time=property(get_time,set_time)

    def get_num_times(self): return len(self.times)
    num_times=property(get_num_times)
    def get_num_fulltimes(self): return len(self.fulltimes)
    num_fulltimes=property(get_num_fulltimes)

    def get_step(self): return self._step
    def set_step(self,step):
        if step<self.fullsteps[0]: self.index=0
        elif step>self.fullsteps[-1]: self.index=-1
        else:
            i=[j for j,sj in enumerate(self.fullsteps) if sj>=step]
            if len(i)>0: self.index=i[0]
    step=property(get_step,set_step)

    def rewind(self):
        """Rewinds to start of listing (without reading any results)"""
        self.seek(0)
        self._index=-1
        
    def first(self): self.index=0
    def last(self): self.index=-1
    def next(self):
        """Find and read next set of results; returns false if at end of listing"""
        more=self.index<self.num_fulltimes-1
        if more: self.index+=1
        return more
    def prev(self):
        """Find and read previous set of results; returns false if at start of listing"""
        more=self.index>0
        if more: self.index-=1
        return more

    def skiplines(self,number=1):
        """Skips specified number of lines in listing file"""
        for i in xrange(number):  self.readline()

    def skipto(self,keyword='',start=1):
        """Skips to line starting (after the leading character) with keyword.  keyword can be either a string
        or a list of strings, in which case it skips to a line starting with any of the specified strings.
        Returns the keyword found, or false if it can't find any of them.  The start parameter specifies which
        character in the line is to be considered the first to search."""
        line=''
        if isinstance(keyword,list): keywords=keyword
        else: keywords=[keyword]
        while not any([line[start:].startswith(kw) for kw in keywords]):
            line=self.readline()
            if line=='': return False
        return [kw for kw in keywords if line[start:].startswith(kw)][0]

    def detect_simulator(self):
        """Detects whether the listing has been produced by AUTOUGH2 or TOUGH2/TOUGH2_MP."""
        self.seek(0)
        simulator={'EEEEE':'AUTOUGH2','ESHORT':'AUTOUGH2','@@@@@':'TOUGH2'}
        tableword=self.skipto(simulator.keys())
        if tableword: self.simulator=simulator[tableword]
        else: self.simulator=None
        if self.simulator=='AUTOUGH2':
            self.setup_pos=self.setup_pos_AUTOUGH2
            self.read_tables=self.read_tables_AUTOUGH2
        if self.simulator=='TOUGH2':
            self.setup_pos=self.setup_pos_TOUGH2
            self.read_tables=self.read_tables_TOUGH2
            self.table_keymap={('ELEM.','INDEX','P'):'EEEEE',('ELEM1','ELEM2','INDEX'):'CCCCC',
                               ('ELEMENT','SOURCE','INDEX'):'GGGGG',('ELEM.','INDEX','X1'):'PPPPP'}

    def setup_short(self):
        """Sets up short_types and short_indices, for handling short output."""
        self.short_types=[]
        self.short_indices={}
        self.seek(0)
        if self.simulator=='AUTOUGH2':
            done=False
            while not done:
                shortkw=self.skipto(['ESHORT','CSHORT','GSHORT'])
                if (shortkw in self.short_types) or not shortkw: done=True
                else:
                    self.short_types.append(shortkw)
                    self.short_indices[shortkw]={}
                    self.skipto(shortkw)
                    self.skiplines(2)
                    indexpos=self.readline().index('INDEX')
                    self.skiplines()
                    endtable=False
                    lineindex=0
                    while not endtable:
                        line=self.readline()
                        if line[1:].startswith(shortkw): endtable=True
                        else:
                            index=int(line[indexpos:indexpos+5])-1
                            self.short_indices[shortkw][index]=lineindex
                        lineindex+=1
        # (no SHORT output for TOUGH2)

    def setup_pos_AUTOUGH2(self):
        """Sets up _pos list for AUTOUGH2 listings, containing file position at the start of each set of results.
        Also sets up the times and steps arrays."""
        self.seek(0)
        # set up pos,times, steps and short arrays:
        self._fullpos,self._pos,self._short=[],[],[]
        fullt,fulls,t,s=[],[],[],[]
        keywords=['EEEEE']
        if len(self.short_types)>0: keywords.append(self.short_types[0])
        endfile=False
        while not endfile:
            kwfound=self.skipto(keywords)
            if kwfound:
                self._pos.append(self.tell())
                self.read_header_AUTOUGH2()
                if kwfound=='EEEEE': # full results
                    self._fullpos.append(self._pos[-1])
                    fullt.append(self.time)
                    fulls.append(self.step)
                    self._short.append(False)
                else: self._short.append(True)
                t.append(self.time)
                s.append(self.step)
                self.readline()
                self.skipto(kwfound) # to end of table
            else: endfile=True
        self.times=np.array(t)
        self.steps=np.array(s)
        self.fulltimes=np.array(fullt)
        self.fullsteps=np.array(fulls)

    def setup_pos_TOUGH2(self):
        """Sets up _pos list for TOUGH2 listings, containing file position at the start of each set of results.
        Also sets up the times and steps arrays."""
        self.seek(0)
        # set up pos,times, steps and short arrays:
        self._fullpos,self._pos=[],[]
        t,s=[],[]
        endfile=False
        while not endfile:
            lf_found=self.skipto('\f',0)
            if lf_found:
                pos=self.tell()
                self.skiplines(2)
                line=self.readline()
                if 'OUTPUT DATA AFTER' in line:
                    self._pos.append(pos)
                    self._fullpos.append(pos)
                    self.seek(pos)
                    self.read_header_TOUGH2()
                    t.append(self.time)
                    s.append(self.step)
                    self.skiplines(2)
                    self.skipto('@@@@@')
            else: endfile=True
        self.times=np.array(t)
        self.steps=np.array(s)
        self.fulltimes=np.array(t)
        self.fullsteps=np.array(s)
        self._short=[False for p in self._pos]

    def setup_tables(self):
        self.tablename={'EEEEE':'element','CCCCC':'connection','GGGGG':'generation','PPPPP':'primary'}
        self._table={}
        if self.simulator=='AUTOUGH2':self.setup_tables_AUTOUGH2()
        else: self.setup_tables_TOUGH2()
        if 'EEEEE' in self._table: self.element=self._table['EEEEE']
        else: self.element=None
        if 'CCCCC' in self._table: self.connection=self._table['CCCCC']
        else: self.connection=None
        if 'GGGGG' in self._table: self.generation=self._table['GGGGG']
        else: self.generation=None
        if 'PPPPP' in self._table: self.primary=self._table['PPPPP']
        else: self.primary=None
        
    def setup_tables_AUTOUGH2(self):
        """Sets up configuration of element, connection and generation tables for AUTOUGH2 listings."""
        keyword='EEEEE'
        self.seek(self._fullpos[0])
        while keyword in self.tablename:
            self.read_header_AUTOUGH2()
            self.setup_table_AUTOUGH2(keyword)
            keyword=self.readline()[1:6]

    def setup_tables_TOUGH2(self):
        """Sets up configuration of element, connection and generation tables for TOUGH2 listings."""
        keyword='EEEEE'
        self.seek(self._fullpos[0])
        self.read_header_TOUGH2()
        while keyword in self.tablename:
            self.setup_table_TOUGH2(keyword)
            self.skipto('\f',0)
            pos=self.tell()
            line=self.readline().strip()
            if (line==self.title):
                if (self.num_times==1) or ((self.num_times>1) and (pos<self._fullpos[1])):
                    self.skiplines(3)
                    headerpos=self.tell()
                    headers=tuple(self.readline().strip().split()[0:3])
                    keyword=self.table_keymap[headers]
                    self.seek(headerpos)
                else: keyword=''
            else: keyword=''

    def valid_spaced_blockname(self,name):
        """Tests if a 7-character string is a valid blockname with spaces around it.  Used to detect positions of table keys."""
        return (name[0]==name[6]==' ') and valid_blockname(name[1:6])

    def setup_table_AUTOUGH2(self,keyword):
        """Sets up table from AUTOUGH2 listing file."""
        self.skiplines(3)
        # Read column names (joining lowercase words to previous names):
        headline=self.readline()
        strs=headline.strip().split()
        nkeys=strs.index('INDEX')
        rows,cols=[],[]
        for s in strs[nkeys+1:]:
            if s[0]==s[0].upper(): cols.append(s)
            else: cols[-1]+=' '+s
        self.readline()
        line=self.readline()
        # Double-check number of columns:
        start=headline.index('INDEX')+5
        nvalues=len([s for s in line[start:].strip().split()])
        if (len(cols)==nvalues):
            # work out positions of keys in line:
            keypos=[]
            pos=0
            for k in xrange(nkeys):
                while not self.valid_spaced_blockname(line[pos:pos+7]): pos+=1
                while self.valid_spaced_blockname(line[pos:pos+7]): pos+=1
                keypos.append(pos)
            # determine row names:
            while line[1:6]<>keyword:
                keyval=[fix_blockname(line[kp:kp+5]) for kp in keypos]
                if len(keyval)>1: keyval=tuple(keyval)
                else: keyval=keyval[0]
                rows.append(keyval)
                line=self.readline()
            self._table[keyword]=listingtable(cols,rows)
            self.readline()
        else:
            print 'Error parsing '+self.tablename[keyword]+' table columns: table not created.'

    def setup_table_TOUGH2(self,keyword):
        """Sets up table from TOUGH2 listing file."""
        # Read column names (joining 'RATE' items to previous names):
        headline=self.readline()
        strs=headline.strip().split()
        nkeys=strs.index('INDEX')
        while self.readline().strip(): pass
        rows,cols=[],[]
        for s in strs[nkeys+1:]:
            if s<>'RATE': cols.append(s)
            else: cols[-1]+=' '+s
        line=self.readline()
        # work out positions of keys in line:
        keypos=[]
        pos=0
        for k in xrange(nkeys):
            while not self.valid_spaced_blockname(line[pos:pos+7]): pos+=1
            while self.valid_spaced_blockname(line[pos:pos+7]): pos+=1
            keypos.append(pos)
        # work out position of index:
        index_pos=[keypos[-1]+5]
        pos=line.find('.')
        c=line[pos-2]
        if c in [' ','-']: index_pos.append(pos-2)
        elif c.isdigit(): index_pos.append(pos-1)
        # determine row names:
        longest_line=line
        rowdict={}
        count,index=0,0
        while line.strip():
            keyval=[fix_blockname(line[kp:kp+5]) for kp in keypos]
            if len(keyval)>1: keyval=tuple(keyval)
            else: keyval=keyval[0]
            indexstr=line[index_pos[0]:index_pos[1]]
            try: index=int(indexstr)-1
            except ValueError: index+=1    # to handle overflow (****) in index field: assume indices continue
            rowdict[index]=(count,keyval)  # use a dictionary to deal with duplicate row indices (TOUGH2_MP)
            line=self.readline(); count+=1
            if line.startswith('\f'): # extra headers in the middle of TOUGH2 listings
                while self.readline().strip(): count+=1
                line=self.readline(); count+=2
            if len(line.strip())>len(longest_line): longest_line=line
        # sort rows (needed for TOUGH2_MP):
        indices=rowdict.keys(); indices.sort()
        row_line=[rowdict[index][0] for index in indices]
        rows=[rowdict[index][1] for index in indices]
        # determine row parsing format:
        line=longest_line
        start=keypos[-1]+5
        numpos=[]
        p,done=start,False
        while not done:
            pos=line.find('.',p)
            if pos>2:
                c=line[pos-2]
                if c in [' ','-']: numpos.append(pos-2)
                elif c.isdigit(): numpos.append(pos-1)
                p=pos+1
            else: done=True
        numpos.append(len(line))
        row_format={'key':keypos,'index':keypos[-1]+5,'values':numpos}
        self._table[keyword]=listingtable(cols,rows,row_format,row_line)

    def read_header_AUTOUGH2(self):
        """Reads header info (title and time data) for one set of AUTOUGH2 listing results."""
        self.title=self.readline().strip()
        line=self.readline()
        istart,iend=string.find(line,'AFTER')+5,string.find(line,'TIME STEPS')
        self._step=int(line[istart:iend])
        istart=iend+10
        iend=string.find(line,'SECONDS')
        self._time=fortran_float(line[istart:iend])
        self.readline()

    def read_header_TOUGH2(self):
        """Reads header info (title and time data) for one set of TOUGH2 listing results."""
        self.title=self.readline().strip()
        self.skiplines(6)
        vals=self.readline().split()
        self._time=fortran_float(vals[0])
        self._step=int(vals[1])
        self.skiplines(3)

    def read_tables_AUTOUGH2(self):
        keyword='EEEEE'
        while keyword in self.tablename:
            self.read_header_AUTOUGH2()
            self.read_table_AUTOUGH2(keyword)
            keyword=self.readline()[1:6]

    def read_tables_TOUGH2(self):
        keyword='EEEEE'
        self.read_header_TOUGH2()
        while keyword in self.tablename:
            self.read_table_TOUGH2(keyword)
            self.skipto('\f',0)
            pos=self.tell()
            if (self.num_fulltimes>1) and (self.index<self.num_fulltimes-1):
                if pos>=self._fullpos[self.index+1]: keyword=''
            if keyword<>'':
                line=self.readline().strip()
                if (line==self.title):
                    self.skiplines(3)
                    headers=tuple(self.readline().strip().split()[0:3])
                    keyword=self.table_keymap[headers]
                else: keyword=''

    def read_table_AUTOUGH2(self,keyword):
        self.skiplines(3)
        start=self.readline().index('INDEX')+5
        self.readline()
        line=self.readline()
        row=0
        while line[1:6]<>keyword: 
            vals=[fortran_float(s) for s in line[start:].strip().split()]
            self._table[keyword][row]=vals
            row+=1
            line=self.readline()
        self.readline()

    def read_table_line_TOUGH2(self,keyword,fmt,line):
        """Reads values from a line in a TOUGH2 listing, given the keyword, and format."""
        vals=[fortran_float(line[fmt['values'][i]:fmt['values'][i+1]]) for i in xrange(len(fmt['values'])-1)]
        num_missing=self._table[keyword].num_columns-len(vals)
        for i in xrange(num_missing): vals.append(0.0)
        return vals
        
    def read_table_TOUGH2(self,keyword):
        fmt=self._table[keyword].row_format
        while self.readline().strip(): pass
        line=self.readline()
        while line.strip():
            key=self._table[keyword].key_from_line(line)
            self._table[keyword][key]=self.read_table_line_TOUGH2(keyword,fmt,line)
            line=self.readline()
            if line.startswith('\f'): # extra headers in the middle of TOUGH2 listings
                while self.readline().strip(): pass
                line=self.readline()

    def history(self,selection):
        """Returns time histories for specified selection of table type, names (or indices) and column names.
           Table type is specified as 'e','c','g' or 'p' (upper or lower case) for element table,
           connection table, generation table or primary table respectively."""
        # This can obviously be done much more simply using next(), and accessing self._table,
        # but that is too slow for large listing files.  This method reads only the required data lines
        # in each table.
        if isinstance(selection,tuple): selection=[selection] # if input just one tuple rather than a list of them
        hist=[[] for s in selection]
        sel=[]
        for (tt,key,h) in selection:  # convert keys to indices as necessary
            keyword=tt.upper()*5
            if isinstance(key,int): index=key
            else: index=self._table[keyword]._row[key]
            if self._table[keyword].row_line: index=self._table[keyword].row_line[index] # find line index if needed
            sel.append((tt,index,h))
        tables=list(set([tt for (tt,i,h) in sel]))
        tablelist=[t for t in ['e','c','p','g'] if t in tables]  # need to retain e,c,g,p order
        tagselection=[(tt,i,h,sindex) for sindex,(tt,i,h) in enumerate(sel)]
        tableselection={}
        shortindex={}
        for table in tablelist:
            tableselection[table]=[(i,h,sindex) for (tt,i,h,sindex) in tagselection if tt==table]
            #  work out where (if at all) the items are in short output:
            for (i,h,sindex) in tableselection[table]:
                shortkw=table.upper()+'SHORT'
                if shortkw in self.short_types:
                    if i in self.short_indices[shortkw]: shortindex[sindex]=self.short_indices[shortkw][i]
                    else: shortindex[sindex]=None
                else: shortindex[sindex]=None
            tableselection[table].sort()
        self.rewind()   

        AUTOUGH2=(self.simulator=='AUTOUGH2')
        for ipos,pos in enumerate(self._pos):
            self.seek(pos)
            short=self._short[ipos]
            if short: startkeyword=self.short_types[0]
            else: startkeyword='EEEEE'
            for tt in tablelist:
                longkeyword=tt.upper()*5
                keyword=[longkeyword,tt.upper()+'SHORT'][short]
                if not (short and not (keyword in self.short_types)):
                    # find start of correct table:
                    if AUTOUGH2:
                        if keyword<>startkeyword: self.skipto(keyword)
                        self.skipto(keyword) # start of table
                        self.skiplines(2)
                        start=self.readline().index('INDEX')+5
                    else: # TOUGH2
                        kw=''
                        while kw<>keyword:
                            if keyword==startkeyword:
                                self.skipto('@@@@@'); self.skipto('@@@@@')
                                self.readline()
                            else:
                                self.skipto('\f',0)
                                self.skiplines(4)
                            headers=tuple(self.readline().strip().split()[0:3])
                            kw=self.table_keymap[headers]
                        fmt=self._table[keyword].row_format
                    while self.readline().strip(): pass
                    index=0
                    line=self.readline()
                    for (itemindex,colname,sindex) in tableselection[tt]:
                        lineindex=[itemindex,shortindex[sindex]][short]
                        if lineindex<>None:
                            for k in xrange(lineindex-index): line=self.readline()
                            index=lineindex
                            if AUTOUGH2: vals=[fortran_float(s) for s in line[start:].strip().split()]
                            else: vals=self.read_table_line_TOUGH2(keyword,fmt,line)
                            valindex=self._table[longkeyword]._col[colname]
                            hist[sindex].append(vals[valindex])
                    if AUTOUGH2: self.skipto(keyword)  # end of table
                    else: self.skipto('@@@@@') # TOUGH2

        self.first()
        result=[([self.times,self.fulltimes][shortindex[sindex]==None],np.array(h)) for sindex,h in enumerate(hist)]
        if len(result)==1: result=result[0]
        return result

    def get_reductions(self):
        """Returns a list of time step indices at which the time step is reduced, and the blocks at which the maximum
        residual occurred prior to the reduction."""
        self.rewind()
        line,lastline='',''
        keyword="+++++++++ REDUCE TIME STEP"
        keyend=len(keyword)+1
        rl=[]
        finished=False
        while not finished:
            while line[1:keyend]<>keyword:
                lastline=line
                line=self.readline()
                if not line:
                    finished=True
                    break
            if not finished:
                lowerlastline=lastline.lower()
                eltindex=lowerlastline.find('element')
                if eltindex>0:
                    if lowerlastline.find('eos cannot find parameters')>=0: space=9
                    else: space=8
                    blockname=fix_blockname(lastline[eltindex+space:eltindex+space+5])
                    brackindex,comindex=line.find('('),line.find(',')
                    timestep=int(line[brackindex+1:comindex])
                    rl.append((timestep,blockname))
                lastline=line
                line=self.readline()
                if not line: finished=True
        return rl
    reductions=property(get_reductions)

    def get_difference(self,indexa=None,indexb=None):
        """Returns dictionary of maximum differences, and locations of difference, of all element table properties between two sets of results.
        If both indexa and indexb are provided, the result is the difference between these two result indices.  If only one index is given, the
        result is the difference between the given index and the one before that.  If neither are given, the result is the difference between
        the last and penultimate sets of results."""
        from copy import deepcopy
        keyword='EEEEE'
        if indexa == None: self.last()
        else: self.set_index(indexa)
        results2=deepcopy(self._table[keyword])
        if indexb == None: self.prev()
        else: self.set_index(indexb)
        results1=self._table[keyword]
        cvg={}
        for name in results1.column_name:
            iblk=np.argmax(abs(results2[name]-results1[name]))
            blkname=results1.row_name[iblk]
            diff=results2[name][iblk]-results1[name][iblk]
            cvg[name]=(diff,blkname)
        return cvg
    convergence=property(get_difference)

    def get_vtk_data(self,geo,grid=None,flows=False,flux_matrix=None,geo_matches=True):
        """Returns dictionary of VTK data arrays from listing file at current time.  If flows is True, average flux vectors
        are also calculated from connection data at the block centres."""
        from vtk import vtkFloatArray
        natm=geo.num_atmosphere_blocks
        nele=geo.num_underground_blocks
        arrays={'Block':{},'Node':{}}
        for name in self.element.column_name: arrays['Block'][name]=vtkFloatArray()
        flownames=[]
        if flows:
            if flux_matrix==None: flux_matrix=grid.flux_matrix(geo)
            flownames=[name for name in self.connection.column_name if (name.endswith('flow') or name.startswith('FLO'))]
            for name in flownames: arrays['Block'][name]=vtkFloatArray()
        array_length={'Block':nele,'Node':0}
        array_data={'Block':{},'Node':{}}
        for array_type,array_dict in arrays.items():
            for name,array in array_dict.items():
                if name in flownames:
                    array.SetName(name+'/area')
                    array.SetNumberOfComponents(3)
                    array.SetNumberOfTuples(array_length[array_type])
                    array_data[array_type][name]=flux_matrix*self.connection[name]
                else:
                    array.SetName(name)
                    array.SetNumberOfComponents(1)
                    array.SetNumberOfValues(array_length[array_type])
                    if geo_matches: array_data[array_type][name]=self.element[name][natm:] # faster
                    else: array_data[array_type][name]=np.array([self.element[blk][name] for blk in geo.block_name_list[natm:]]) # more flexible
        for array_type,data_dict in array_data.items():
            for name,data in data_dict.items():
                if name in flownames:
                    for iblk in xrange(nele):
                        arrays[array_type][name].SetTuple3(iblk,data[3*iblk],data[3*iblk+1],data[3*iblk+2])
                else:    
                    for iblk in xrange(nele):
                        arrays[array_type][name].SetValue(iblk,data[iblk])
        return arrays

    def write_vtk(self,geo,filename,grid=None,indices=None,flows=False,wells=False,start_time=0.0,time_unit='s'):
        """Writes VTK files for a vtkUnstructuredGrid object corresponding to the grid in 3D with the listing data,
        with the specified filename, for visualisation with VTK.  A t2grid can optionally be specified, to include rock type
        data as well.  A list of the required time indices can optionally be specified.  If a grid is specified, flows is True,
        and connection data are present in the listing file, approximate average flux vectors are also calculated at the 
        block centres from the connection data."""
        from vtk import vtkXMLUnstructuredGridWriter
        from os.path import splitext
        base,ext=splitext(filename)
        if wells: geo.write_well_vtk()
        geo_matches=geo.block_name_list==self.element.row_name
        doflows=flows and (self.connection<>None) and (grid<>None) and geo_matches
        arrays=geo.vtk_data
        if grid<>None:
            grid_arrays=grid.get_vtk_data(geo)
            for array_type,array_dict in arrays.items():
                array_dict.update(grid_arrays[array_type])
        if doflows: flux_matrix=grid.flux_matrix(geo)
        else: flux_matrix=None
        import xml.dom.minidom
        pvd=xml.dom.minidom.Document()
        vtkfile=pvd.createElement('VTKFile')
        vtkfile.setAttribute('type','Collection')
        pvd.appendChild(vtkfile)
        collection=pvd.createElement('Collection')
        initial_index=self.index
        if indices==None: indices=range(self.num_fulltimes)
        timescales={'s':1.0,'h':3600.,'d':3600.*24,'y':3600.*24*365.25}
        if time_unit in timescales: timescale=timescales[time_unit]
        else: timescale=1.0
        writer=vtkXMLUnstructuredGridWriter()
        for i in indices:
            self.index=i
            t=start_time+self.time/timescale
            filename_time=base+'_'+str(i)+'.vtu'
            results_arrays=self.get_vtk_data(geo,grid,flows=doflows,flux_matrix=flux_matrix,geo_matches=geo_matches)
            for array_type,array_dict in arrays.items():
                array_dict.update(results_arrays[array_type])
            vtu=geo.get_vtk_grid(arrays)
            writer.SetFileName(filename_time)
            writer.SetInput(vtu)
            writer.Write()
            dataset=pvd.createElement('DataSet')
            dataset.setAttribute('timestep',str(t))
            dataset.setAttribute('file',filename_time)
            collection.appendChild(dataset)
        vtkfile.appendChild(collection)
        pvdfile=open(base+'.pvd','w')
        pvdfile.write(pvd.toprettyxml())
        pvdfile.close()
        self.index=initial_index

    def add_side_recharge(self,geo,dat):
        """Adds side recharge generators to a TOUGH2 data object for a production run,
        calculated according to the final results in the listing.  These generators represent side
        inflows due to pressure changes in the blocks on the model's horizontal boundaries.
        Recharge generators are given the names of their blocks- any existing generators with the same
        names will be overwritten."""
        from IAPWS97 import cowat,sat,visc
        from geometry import line_projection
        from t2data import t2generator
        initial_index=self.index
        keyword={'AUTOUGH2':{'P':'Pressure','T':'Temperature'},'TOUGH2':{'P':'P','T':'T'}}
        self.last()
        bdy_nodes=geo.boundary_nodes
        for blk in dat.grid.blocklist[geo.num_atmosphere_blocks:]:
            colname=geo.column_name(blk.name)
            if colname in geo.column:
                col=geo.column[colname]
                if col.num_neighbours<col.num_nodes:
                    k=0.5*np.sum(blk.rocktype.permeability[0:2])
                    p0=self.element[blk.name][keyword[self.simulator]['P']]
                    t0=self.element[blk.name][keyword[self.simulator]['T']]
                    rho,u=cowat(t0,p0)
                    h=u+p0/rho
                    Ps=sat(t0)
                    xnu=visc(rho,t0)/rho
                    coef=0.
                    for iface in xrange(col.num_nodes):
                        facenode=[col.node[i] for i in [iface,(iface+1)%col.num_nodes]]
                        if all([node in bdy_nodes for node in facenode]):
                            side_length=np.linalg.norm(facenode[1].pos-facenode[0].pos)
                            height=blk.volume/col.area
                            area=side_length*height
                            facepos=line_projection(col.centre,[node.pos for node in facenode])
                            dist=np.linalg.norm(col.centre-facepos)
                            coef+=0.5*area*k/(xnu*dist) # recharge coefficient
                    gen_name=blk.name
                    dat.add_generator(t2generator(gen_name,blk.name,type='RECH',gx=coef,ex=h,hg=p0))
        self.index=initial_index
