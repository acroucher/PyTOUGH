"""For reading TOUGH2 listing files."""

"""
Copyright 2012 University of Auckland.

This file is part of PyTOUGH.

PyTOUGH is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

PyTOUGH is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with PyTOUGH.  If not, see <http://www.gnu.org/licenses/>."""

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
    def __init__(self,cols,rows,row_format=None,row_line=None,num_keys=1,allow_reverse_keys=False):
        """The row_format parameter is a dictionary with three keys, 'key','index' and 'values'.  These contain the positions,
        in each row of the table, of the start of the keys, index and data fields.  The row_line parameter is a list containing,
        for each row of the table, the number of lines before it in the listing file, from the start of the table.  This is
        needed for TOUGH2_MP listing files, in which the rows are not in index order and can also be duplicated."""
        self.column_name=cols
        self.row_name=rows
        self.row_format=row_format
        self.row_line=row_line
        self.num_keys=num_keys
        self.allow_reverse_keys=allow_reverse_keys
        self._col=dict([(c,i) for i,c in enumerate(cols)])
        self._row=dict([(r,i) for i,r in enumerate(rows)])
        self._data=np.zeros((len(rows),len(cols)),float64)
    def __repr__(self): return repr(self.column_name)+'\n'+repr(self._data)
    def __getitem__(self,key):
        if isinstance(key,int): return dict(zip(['key']+self.column_name,[self.row_name[key]]+list(self._data[key,:])))
        else:
            if key in self.column_name: return self._data[:,self._col[key]]
            elif key in self.row_name:
                rowindex=self._row[key]
                return dict(zip(['key']+self.column_name,[self.row_name[rowindex]]+list(self._data[rowindex,:])))
            elif len(key)>1 and self.allow_reverse_keys:
                revkey=key[::-1] # try reversed key for multi-key tables
                if revkey in self.row_name:
                    rowindex=self._row[revkey]
                    return dict(zip(['key']+self.column_name,[self.row_name[rowindex][::-1]]+list(-self._data[rowindex,:])))
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
    def rows_matching(self,pattern,index=0,match_any=False):
        """Returns rows in the table with keys matching the specified regular expression pattern
        string.
        For tables with multiple keys, pattern can be a list or tuple of regular expressions.  If
        a single string pattern is given for a multiple-key table, the pattern is matched on the
        index'th key (and any value of the other key- unless match_any is used; see below).
        If match_any is set to True, rows are returned with keys matching any of the specified
        patterns (instead of all of them).  If this option is used in conjunction with a single
        string pattern, the specified pattern is applied to all keys."""
        from re import search
        if self.num_keys==1: return [self[key] for key in self.row_name if search(pattern,key)]
        else: 
            if isinstance(pattern,str): pattern=[pattern]
            else: pattern=list(pattern)
            if len(pattern)<self.num_keys:
                if match_any: default=[pattern[0]]
                else: default=['.*']
                if 0<=index<=self.num_keys:
                    pattern=default*index+pattern+default*(self.num_keys-1-index)
                else: return []
            combine=[all,any][match_any]
            return [self[key] for key in self.row_name if combine([search(p,n) for p,n in zip(pattern,key)])]

class t2listing(file):
    """Class for TOUGH2 listing file.  The element, connection and generation tables can be accessed
       via the element, connection and generation fields.  (For example, the pressure in block 'aa100' is
       given by element['aa100']['Pressure'].)  It is possible to navigate through time in the listing by 
       using the next() and prev() functions to step through, or using the first() and last() functions to go to 
       the start or end, or to set the index, step (model time step number) or time properties directly."""
    def __init__(self,filename=None):
        super(t2listing,self).__init__(filename,'rU')
        self.detect_simulator()
        if self.simulator==None: print 'Could not detect simulator type.'
        else:
            self.setup_short()
            self.setup_pos()
            if self.num_fulltimes>0:
                self._index=0
                self.setup_tables()
                self.set_table_attributes()
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

    def get_table_names(self):
        names=self._table.keys()
        names.sort()
        return names
    table_names=property(get_table_names)

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
        """Skips to line starting  with keyword.  keyword can be either a string or a list of strings, in which case
        it skips to a line starting with any of the specified strings.
        Returns the keyword found, or false if it can't find any of them.  The start parameter specifies which
        character in the line is to be considered the first to search."""
        line=''
        if isinstance(keyword,list): keywords=keyword
        else: keywords=[keyword]
        while not any([line[start:].startswith(kw) for kw in keywords]):
            line=self.readline()
            if line=='': return False
        return [kw for kw in keywords if line[start:].startswith(kw)][0]

    def skip_to_nonblank(self):
        """Skips to start of next non-blank line."""
        pos=self.tell()
        while not self.readline().strip(): pos=self.tell()
        self.seek(pos)

    def skip_to_blank(self):
        """Skips to start of next blank line."""
        pos=self.tell()
        while self.readline().strip(): pos=self.tell()
        self.seek(pos)
        
    def skip_over_next_blank(self):
        """Skips past next blank line."""
        while self.readline().strip(): pass

    def detect_simulator(self):
        """Detects whether the listing has been produced by AUTOUGH2, TOUGH2/TOUGH2_MP or TOUGH+, and sets some internal methods
        according to the simulator type."""
        self.seek(0)
        simulator={'EEEEE':'AUTOUGH2','ESHORT':'AUTOUGH2','BBBBB':'AUTOUGH2','@@@@@':'TOUGH2','=====':'TOUGH+'}
        line=' '
        while not ('output data after' in line or 'output after' in line or line==''): line=self.readline().lower()
        if line=='': self.simulator=None
        else:
            self.readline()
            line=self.readline()
            linechars=line[1:6]
            if linechars in simulator.keys(): self.simulator=simulator[linechars]
            else: self.simulator=None
            if self.simulator:
                # Set internal methods according to simulator type:
                simname=self.simulator.replace('+','plus')
                internal_fns=['setup_pos','table_type','setup_table','setup_tables','read_header','read_table','next_table',
                              'read_tables','skip_to_table','read_table_line']
                for fname in internal_fns:
                    fname_sim=fname+'_'+simname
                    if simname=='TOUGHplus' and not hasattr(self,fname_sim): fname_sim=fname_sim.replace('plus','2')
                    setattr(self,fname,getattr(self,fname_sim))

    def table_type_AUTOUGH2(self,keyword):
        """Returns AUTOUGH2 table name based on the 5-character keyword read at the top of the table."""
        keytable={'EEEEE':'element','CCCCC':'connection','GGGGG':'generation'}
        if keyword in keytable: return keytable[keyword]
        else: return None

    def table_type_TOUGH2(self,headers):
        """Returns TOUGH2 table name based on a tuple of the first three column headings."""
        if headers[0:2]==('ELEM.','INDEX'):
            if headers[2]=='P': return 'element'
            elif headers[2]=='X1': return 'primary'
        else:
            keytable={('ELEM1','ELEM2','INDEX'):'connection',('ELEMENT','SOURCE','INDEX'):'generation'}
            if headers in keytable: return keytable[headers]
        return None

    def table_type_TOUGHplus(self,headers):
        """Returns TOUGH+ table name based on a tuple of the first three column headings."""
        if headers[0:2]==('ELEM','INDEX'):
            if headers[2]=='X1': return 'primary'
            else: return 'element'
        else:
            keytable={('ELEM1','ELEM2','INDEX'):'connection',('ELEMENT','SOURCE','INDEX'):'generation'}
            if headers in keytable: return keytable[headers]
        return None
        
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
        # (no SHORT output for TOUGH2 or TOUGH+)

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

    def setup_pos_TOUGHplus(self):
        """Sets up _pos list for TOUGH+ listings, containing file position at the start of each set of results.
        Also sets up the times and steps arrays."""
        self.seek(0)
        # set up pos,times, steps and short arrays:
        self._fullpos,self._pos=[],[]
        t,s=[],[]
        endfile=False
        while not endfile:
            line=' '
            while not (line.lstrip().startswith('Output data after') or line==''): line=self.readline()
            if line<>'':
                self.skipto('TOTAL TIME',2)
                pos=self.tell()
                self._pos.append(pos)
                self._fullpos.append(pos)
                self.read_header_TOUGHplus()
                t.append(self.time)
                s.append(self.step)
                self.skipto('@@@@@')
            else: endfile=True
        self.times=np.array(t)
        self.steps=np.array(s)
        self.fulltimes=np.array(t)
        self.fullsteps=np.array(s)
        self._short=[False for p in self._pos]

    def set_table_attributes(self):
        """Makes tables in self._table accessible as attributes."""        
        for key,table in self._table.iteritems(): setattr(self,key,table)
        
    def setup_tables_AUTOUGH2(self):
        """Sets up configuration of element, connection and generation tables."""
        self._table={}
        tablename='element'
        self.seek(self._fullpos[0])
        while tablename:
            self.read_header()
            self.setup_table(tablename)
            tablename=self.next_table()

    def setup_tables_TOUGH2(self):
        self._table={}
        tablename='element'
        self.seek(self._fullpos[0])
        self.read_header() # only one header at each time
        while tablename:
            self.setup_table(tablename)
            tablename=self.next_table()

    def setup_tables_TOUGHplus(self):
        self.read_title_TOUGHplus()
        self._table={}
        tablename='element'
        self.seek(self._fullpos[0])
        self.read_header() # only one header at each time
        nelt_tables=0 # can have multiple element tables
        while tablename:
            self.setup_table(tablename)
            tablename=self.next_table()
            if tablename=='element':
                nelt_tables+=1
                tablename+=str(nelt_tables)

    def next_table_AUTOUGH2(self):
        """Goes to start of next table at current time and returns its type- or None if there are no more."""
        keyword=self.readline()[1:6]
        return self.table_type(keyword)

    def next_table_TOUGH2(self):
        if self.skipto('\f',0):
            pos=self.tell()
            if (self.num_fulltimes>1) and (self.index<self.num_fulltimes-1):
                if pos>=self._fullpos[self.index+1]: return None
            line=self.readline().strip()
            if (line==self.title):
                self.skiplines(3)
                headpos=self.tell()
                headers=tuple(self.readline().strip().split()[0:3])
                self.seek(headpos)
                return self.table_type(headers)
            else: return None
        else: return None

    def next_table_TOUGHplus(self):
        if self.skipto('_____',0):
            self.readline()
            pos=self.tell()
            if (self.num_fulltimes>1) and (self.index<self.num_fulltimes-1):
                if pos>=self._fullpos[self.index+1]: return None
            headpos=self.tell()
            headers=tuple(self.readline().strip().split()[0:3])
            self.seek(headpos)
            return self.table_type(headers)
        else: return None

    def skip_to_table_AUTOUGH2(self,tablename,last_tablename,nelt_tables):
        """Skips forwards to headers of table with specified name at the current time."""
        if self._short[self._index]: keyword=tablename[0].upper()+'SHORT'
        else: keyword=tablename[0].upper()*5
        self.skipto(keyword)
        if tablename<>'element': self.skipto(keyword)
        self.skip_to_blank()
        self.skip_to_nonblank()

    def skip_to_table_TOUGH2(self,tablename,last_tablename,nelt_tables):
        if last_tablename==None:
            for i in xrange(2): self.skipto('@@@@@')
            self.skip_to_nonblank()
            tname='element'
        else: tname=last_tablename
        while tname<>tablename:
            self.skipto('@@@@@')
            tname=self.next_table_TOUGH2()

    def skip_to_table_TOUGHplus(self,tablename,last_tablename,nelt_tables):
        if last_tablename==None:
            self.skipto('=====',0)
            self.skip_to_nonblank()
            tname='element'
            nelt_tables=0
        else: tname=last_tablename
        while tname<>tablename:
            if tname=='primary': keyword='_____'
            else: keyword='@@@@@'
            self.skipto(keyword,0)
            tname=self.next_table_TOUGHplus()
            if tname=='element':
                nelt_tables+=1
                tname+=str(nelt_tables)

    def key_positions(self,line,nkeys):
        """Returns detected positions of keys in the start of a table line.  This works on the assumption
        that key values must have a digit present in their last character.  It searches backwards from the 
        end of the position of INDEX in the header line- sometimes keys can overlap into this area."""
        keylength=5
        keypos=[]
        pos=len(line)-1
        indexstr='INDEX'
        indexpos=pos-len(indexstr)
        # skip over index:
        while line[pos]==' ' and pos>indexpos: pos-=1
        while line[pos]<>' ' and pos>indexpos: pos-=1
        for k in xrange(nkeys):
            while not line[pos].isdigit() and pos>=keylength: pos-=1
            pos-=(keylength-1)
            if valid_blockname(line[pos:pos+keylength]):
                keypos.append(pos)
                pos-=1
            else: return None
        keypos.reverse()
        if len(keypos)<>nkeys: return None
        else: return keypos

    def setup_table_AUTOUGH2(self,tablename):
        """Sets up table from AUTOUGH2 listing file."""
        keyword=tablename[0].upper()*5
        self.skiplines(3)
        # Read column names (joining lowercase words to previous names):
        headline=self.readline()
        strs=headline.strip().split()
        indexstr='INDEX'
        nkeys=strs.index(indexstr)
        rows,cols=[],[]
        for s in strs[nkeys+1:]:
            if s[0]==s[0].upper(): cols.append(s)
            else: cols[-1]+=' '+s
        self.readline()
        line=self.readline()
        # Double-check number of columns:
        indexpos=headline.index(indexstr)
        start=indexpos+len(indexstr)
        nvalues=len([s for s in line[start:].strip().split()])
        if (len(cols)==nvalues):
            keypos=self.key_positions(line[:start],nkeys)
            if keypos:
                # determine row names:
                while line[1:6]<>keyword:
                    keyval=[fix_blockname(line[kp:kp+5]) for kp in keypos]
                    if len(keyval)>1: keyval=tuple(keyval)
                    else: keyval=keyval[0]
                    rows.append(keyval)
                    line=self.readline()
                row_format={'values':[start]}
                allow_rev=tablename=='connection'
                self._table[tablename]=listingtable(cols,rows,row_format,num_keys=nkeys,allow_reverse_keys=allow_rev)
                self.readline()
            else: print 'Error parsing '+tablename+' table keys: table not created.'
        else:
            print 'Error parsing '+tablename+' table columns: table not created.'

    def setup_table_TOUGH2(self,tablename):
        """Sets up table from TOUGH2 (or TOUGH+) listing file."""
        # Read column names (joining flow items to previous names):
        if self.simulator=='TOUGH2': flow_headers=['RATE']
        else: flow_headers=['Flow','Veloc']
        headline=self.readline()
        strs=headline.strip().split()
        indexstr='INDEX'
        indexpos=headline.index(indexstr)
        start=indexpos+len(indexstr)
        nkeys=strs.index(indexstr)
        self.skip_over_next_blank()
        rows,cols=[],[]
        for s in strs[nkeys+1:]:
            if s in flow_headers: cols[-1]+=' '+s
            else: cols.append(s)
        line=self.readline()
        keypos=self.key_positions(line[:start],nkeys)
        if keypos:
            # work out position of index:
            index_pos=[keypos[-1]+5]
            pos=line.find('.')
            c=line[pos-2]
            if c in [' ','-']: index_pos.append(pos-2)
            elif c.isdigit(): index_pos.append(pos-1)
            # determine row names:
            longest_line=line
            rowdict={}
            count,index=0,-1
            def count_read(count): return self.readline(),count+1
            while line.strip():
                keyval=[fix_blockname(line[kp:kp+5]) for kp in keypos]
                if len(keyval)>1: keyval=tuple(keyval)
                else: keyval=keyval[0]
                indexstr=line[index_pos[0]:index_pos[1]]
                try: index=int(indexstr)-1
                except ValueError: index+=1    # to handle overflow (****) in index field: assume indices continue
                rowdict[index]=(count,keyval)  # use a dictionary to deal with duplicate row indices (TOUGH2_MP)
                line,count=count_read(count)
                if line.startswith('\f'):
                    line,count=count_read(count)
                    if line.strip()==self.title: break # some TOUGH2_MP output ends with \f
                    else: # extra headers in the middle of TOUGH2 listings
                        while line.strip(): line,count=count_read(count)
                        line,count=count_read(count)
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
            allow_rev=tablename=='connection'
            self._table[tablename]=listingtable(cols,rows,row_format,row_line,num_keys=nkeys,allow_reverse_keys=allow_rev)
        else: print 'Error parsing '+tablename+' table keys: table not created.'

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
        self.skipto('@@@@@')
        self.skip_to_nonblank()

    def read_header_TOUGHplus(self):
        """Reads header info (time data) for one set of TOUGH+ listing results."""
        line=self.readline()
        strs=line.split()
        self._time,self._step=float(strs[0]),int(strs[1])
        self.skipto('=====')
        self.skip_to_nonblank()

    def read_title_TOUGHplus(self):
        """Reads simulation title for TOUGH+ listings, at top of file."""
        self.seek(0)
        line=' '
        while not (line.lstrip().startswith('PROBLEM TITLE:') or line==''): line=self.readline()
        if line=='': self.title=''
        else:
            colonpos=line.find(':')
            if colonpos>=0: self.title=line[colonpos+1:].strip()
            else: self.title=''

    def read_tables_AUTOUGH2(self):
        tablename='element'
        while tablename:
            self.read_header()
            self.read_table(tablename)
            tablename=self.next_table()

    def read_tables_TOUGH2(self):
        tablename='element'
        self.read_header() # only one header at each time
        while tablename:
            self.read_table(tablename)
            tablename=self.next_table()

    def read_tables_TOUGHplus(self):
        tablename='element'
        self.read_header() # only one header at each time
        nelt_tables=0
        while tablename:
            self.read_table(tablename)
            tablename=self.next_table()
            if tablename=='element':
                nelt_tables+=1
                tablename+=str(nelt_tables)

    def read_table_AUTOUGH2(self,tablename):
        fmt=self._table[tablename].row_format
        keyword=tablename[0].upper()*5
        self.skip_to_blank()
        self.readline()
        self.skip_to_blank()
        self.skip_to_nonblank()
        line=self.readline()
        row=0
        while line[1:6]<>keyword:
            self._table[tablename][row]=self.read_table_line_AUTOUGH2(line,fmt=fmt)
            row+=1
            line=self.readline()
        self.readline()

    def read_table_line_AUTOUGH2(self,line,num_columns=None,fmt=None):
        start=fmt['values'][0]
        vals=[fortran_float(s) for s in line[start:].strip().split()]        
        return vals

    def read_table_line_TOUGH2(self,line,num_columns,fmt):
        """Reads values from a line in a TOUGH2 listing, given the number of columns, and format."""
        vals=[fortran_float(line[fmt['values'][i]:fmt['values'][i+1]]) for i in xrange(len(fmt['values'])-1)]
        num_missing=num_columns-len(vals)
        for i in xrange(num_missing): vals.append(0.0)
        return vals
        
    def read_table_TOUGH2(self,tablename):
        ncols=self._table[tablename].num_columns
        fmt=self._table[tablename].row_format
        self.skip_to_blank()
        self.skip_to_nonblank()
        line=self.readline()
        while line.strip():
            key=self._table[tablename].key_from_line(line)
            self._table[tablename][key]=self.read_table_line_TOUGH2(line,ncols,fmt)
            line=self.readline()
            if line.startswith('\f'):
                line=self.readline()
                if line.strip()==self.title: break # some TOUGH2_MP output ends with \f
                else: # extra headers in the middle of TOUGH2 listings
                    self.skip_over_next_blank()
                    line=self.readline()

    def history(self,selection,short=True):
        """Returns time histories for specified selection of table type, names (or indices) and column names.
           Table type is specified as 'e','c','g' or 'p' (upper or lower case) for element table,
           connection table, generation table or primary table respectively.  For TOUGH+ results, additional
           element tables may be specified as 'e1' or 'e2'.  If the short parameter is True, results from 
           'short output' (AUTOUGH2 only) are included in the results."""

        # This can obviously be done much more simply using next(), and accessing self._table,
        # but that is too slow for large listing files.  This method reads only the required data lines
        # in each table.

        def tablename_from_specification(tabletype): # expand table specification to table name:
            from string import digits
            namemap={'e':'element','c':'connection','g':'generation','p':'primary'}
            type0=tabletype[0].lower()
            if type0 in namemap:
                name=namemap[type0]
                if tabletype[-1] in digits: name+=tabletype[-1] # additional TOUGH+ element tables
                return name
            else: return None

        def ordered_selection(selection,tables,short_types,short_indices):
            """Given the initial history selection, returns a list of tuples of table name and table selection.  The tables
            are in the same order as they appear in the listing file.  The table selection is a list of tuples of 
            (row index, column name, selection index, short row index) for each table, ordered by row index.  This ordering
            means all data can be read sequentially to make it more efficient.""" 
            converted_selection=[]
            for (tspec,key,h) in selection:  # convert keys to indices as necessary, and expand table names
                tablename=tablename_from_specification(tspec)
                if isinstance(key,int): index=key
                else:
                    index,reverse=None,False
                    if key in tables[tablename].row_name: index=tables[tablename]._row[key]
                    elif len(key)>1 and tables[tablename].allow_reverse_keys:
                        revkey=key[::-1]
                        if revkey in tables[tablename].row_name:
                            index=tables[tablename]._row[revkey]
                            reverse=True
                if index<>None:
                    if tables[tablename].row_line: index=tables[tablename].row_line[index] # find line index if needed
                    ishort=None
                    short_keyword=tspec[0].upper()+'SHORT'
                    if short_keyword in short_types:
                        if index in short_indices[short_keyword]: ishort=short_indices[short_keyword][index]
                    converted_selection.append((tablename,index,ishort,h,reverse))
            tables=list(set([sel[0] for sel in converted_selection]))
            # need to retain table order as in the file:
            tables=[tname for tname in ['element','element1','connection','primary','element2','generation'] if tname in tables]
            tagselection=[(tname,i,ishort,h,rev,sel_index) for sel_index,(tname,i,ishort,h,rev) in enumerate(converted_selection)]
            tableselection=[]
            shortindex={}
            for table in tables:
                tselect=[(i,ishort,h,rev,sel_index) for (tname,i,ishort,h,rev,sel_index) in tagselection if tname==table]
                tselect.sort()
                tableselection.append((table,tselect))
            return tableselection

        old_index=self.index
        if isinstance(selection,tuple): selection=[selection] # if input just one tuple rather than a list of them
        tableselection=ordered_selection(selection,self._table,self.short_types,self.short_indices)
        hist=[[] for s in selection]
        self.rewind()

        for ipos,pos in enumerate(self._pos):
            self.seek(pos)
            self._index=ipos
            is_short=self._short[ipos]
            if not (is_short and not short):
                last_tname=None
                nelt_tables=-1
                for (tname,tselect) in tableselection:
                    if is_short: tablename=tname[0].upper()+'SHORT'
                    else: tablename=tname
                    if not (is_short and not (tablename in self.short_types)):
                        self.skip_to_table(tname,last_tname,nelt_tables)
                        if tname.startswith('element'): nelt_tables+=1
                        self.skip_to_blank()
                        self.skip_to_nonblank()
                        ncols=self._table[tname].num_columns
                        fmt=self._table[tname].row_format
                        index=0
                        line=self.readline()
                        for (itemindex,ishort,colname,reverse,sel_index) in tselect:
                            lineindex=[itemindex,ishort][is_short]
                            if lineindex<>None:
                                for k in xrange(lineindex-index): line=self.readline()
                                index=lineindex
                                vals=self.read_table_line(line,ncols,fmt)
                                valindex=self._table[tname]._col[colname]
                                sgn=[1.,-1.][reverse]
                                hist[sel_index].append(sgn*vals[valindex])
                    last_tname=tname

        self._index=old_index
        result=[([self.times,self.fulltimes][len(h)==self.num_fulltimes],np.array(h)) for sel_index,h in enumerate(hist)]
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
        tablename='element'
        if indexa == None: self.last()
        else: self.set_index(indexa)
        results2=deepcopy(self._table[tablename])
        if indexb == None: self.prev()
        else: self.set_index(indexb)
        results1=self._table[tablename]
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
        elt_tablenames=[key for key in self._table.keys() if key.startswith('element')]
        for tablename in elt_tablenames:
            for name in self._table[tablename].column_name: arrays['Block'][name]=vtkFloatArray()
        flownames=[]
        def is_flowname(name):
            name=name.lower()
            return name.startswith('flo') or name.endswith('flo') or name.endswith('flow') or name.endswith('veloc')
        if flows:
            if flux_matrix==None: flux_matrix=grid.flux_matrix(geo)
            flownames=[name for name in self.connection.column_name if is_flowname(name)]
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
                    for tablename in elt_tablenames:
                        if geo_matches: array_data[array_type][name]=self._table[tablename][name][natm:] # faster
                        else:  # more flexible
                            array_data[array_type][name]=np.array([self._table[tablename][blk][name] for blk in geo.block_name_list[natm:]])
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
        keyword={'AUTOUGH2':{'P':'Pressure','T':'Temperature'},'TOUGH2':{'P':'P','T':'T'},
                 'TOUGH+':{'P':'Pressure','T':'Temperature'}}
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

class t2historyfile(object):
    """Class for TOUGH2 FOFT, COFT and GOFT files (history of element, connection and generator variables)."""

    def __init__(self,filename=None):
        self.filename=filename
        self.empty()
        if self.filename: self.read(filename)

    def empty(self):
        self.simulator=None
        self.type=None
        self._data=[]
        self._row=None
        self.times=[]
        self.keys=[]
        self.key_name=[]
        self.times=[]
        self._keyrows={}
        self.column_name=[]
        self.row_name=[]

    def get_num_keys(self): return len(self.keys)
    num_keys=property(get_num_keys)
    def get_num_times(self): return len(self.times)
    num_times=property(get_num_times)
    def get_num_rows(self): return len(self._data)
    num_rows=property(get_num_rows)
    def get_num_columns(self): return len(self.column_name)
    num_columns=property(get_num_columns)
    def get_keytype(self): return {'FOFT':'block','COFT':'connection','GOFT':'generator',None:'key'}[self.type]
    keytype=property(get_keytype)

    def __repr__(self):
        if self._nkeys>0: nkeystr=str(self.num_keys)
        else: nkeystr='summed'
        return 'History results for '+nkeystr+' '+self.keytype+'s at '+str(self.num_times)+' times'

    def __getitem__(self,key):
        """Returns results for a given key value, in the form of a dictionary with keys corresponding to the
        column names.  Each value is an array of history results for that key and column name- unless the key also
        has the time appended as its third element, in which case the dictionary values are just single floating
        point values for each column."""
        if self.num_rows>0:
            if self._nkeys>0:
                if not isinstance(key,tuple): key=(key,)
                if key in self.keys:
                    keydata=self._data[self._keyrows[key]]
                    return dict([(colname,keydata[:,icol]) for icol,colname in enumerate(self.column_name)])
                elif key in self._row:
                    row=self._data[self._row[key]]
                    return dict([(colname,row[icol]) for icol,colname in enumerate(self.column_name)])
                else: return None
            else: # no keys (e.g. TOUGH+ COFT/GOFT)
                try:
                    icol=self.column_name.index(key)
                    return self._data[:,icol]
                except ValueError: return None
        else: return None

    def read(self,filename):
        """Reads contents of file(s) and stores in appropriate data structures."""
        self._rowindex=0
        from glob import glob
        files=glob(filename)
        configured=False
        for i,fname in enumerate(files):
            self._file=open(fname,'rU')
            header=self._file.readline()
            if header:
                if not configured:
                    self.detect_simulator(header)
                    self.setup_headers(fname,header)
                self.read_data(configured)
                if self.num_columns>0: configured=True
            self._file.close()
        self.finalize_data()

    def detect_simulator(self,header):
        """Detects simulator (TOUGH2 or TOUGH2_MP) from header line."""
        if 'OFT' in header: self.simulator='TOUGH2_MP'
        elif header.startswith('Time ['): self.simulator='TOUGH+'
        else: self.simulator='TOUGH2'
        internal_fns=['setup_headers','read_data']
        simname=self.simulator.replace('+','plus')
        for fname in internal_fns:
            fname_sim=fname+'_'+simname
            setattr(self,fname,getattr(self,fname_sim))
        
    def setup_headers_TOUGH2(self,filename,header):
        """Sets up keys and column headings from given filename and header line, for TOUGH2 output."""
        if filename.endswith('OFT') and len(filename)>=4: self.type=filename[-4:].strip()
        else: self.type=None
        self.time_index=1
        self._nkeys=1
        items=header.strip().split(',')
        if items[-1]=='': del items[-1] # often an extra comma on the end of the lines
        items=items[3:]
        int_index=0
        for i,item in enumerate(items):
            try:
                int_item=int(item)
                int_index=i
                break
            except: pass
        if int_index==0: ncols=len(items)
        else: ncols=int_index
        self.column_name=range(ncols)

    def setup_headers_TOUGH2_MP(self,filename,header):
        """Sets up keys and column headings from given filename and header line, for TOUGH2_MP output."""
        headers=header.strip().split()
        self.type=headers[0] # FOFT, COFT or GOFT
        time_header=[h for h in headers if h.lower().startswith('time')][0]
        self.time_index=headers.index(time_header)
        if self.type=='FOFT':
            self.key_index=self.time_index-1
            self._nkeys=1
        else: # COFT or GOFT
            self.key_index=self.time_index+1
            self._nkeys=2
        self.key_name=headers[self.key_index:self.key_index+self._nkeys]
        prepend_titles, append_titles=['GAS','GENERATION'],['flow']
        startcol=self._nkeys+2
        cols=[]
        i=startcol
        while i<=len(headers)-1:
            title=headers[i]
            if title in prepend_titles:
                cols.append(title+' '+headers[i+1])
                i+=1
            elif title in append_titles:
                cols[-1]+=' '+title
            else: cols.append(title)
            i+=1
        self.column_name=cols
        self.col_start=[header.index(colname) for colname in self.column_name]
        self.key_start=[header.index(key) for key in self.key_name]
        self.time_pos=[header.index(time_header)]
        if self.type=='FOFT':
            self.key_start.append(self.time_pos[0])
            self.time_pos.append(self.col_start[0])
        else:
            self.key_start.append(self.col_start[0])
            self.time_pos.append(self.key_start[0])

    def setup_headers_TOUGHplus(self,filename,header):
        """Sets up keys and column headings from given filename and header line, for TOUGH+ output."""
        headers=header.strip().split('-')
        if filename.endswith('OFT') and len(filename)>=4: self.type=filename[-4:].strip()
        elif '_Time_Series' in filename:
            filetype={'Elem_Time_Series':'FOFT','Conx_Time_Series':'COFT','SS_Time_Series':'GOFT'}
            for key in filetype.keys():
                if key in filename:
                    self.type=filetype[key]
                    break
        else: self.type=None
        if self.type=='FOFT':
            self.time_index=1
            self._nkeys=1
        else:
            self.time_index=0
            self._nkeys=0
        cols=headers[self._nkeys+1:]
        from re import sub
        cols=[sub('\[.*\]','',col).strip() for col in cols] # remove units
        self.column_name=cols

    def read_data_TOUGH2(self,configured):
        """Reads in the data, for TOUGH2 output."""
        self._file.seek(0)
        lines=self._file.readlines()
        for line in lines:
            items=line.strip().split(',')
            if items[-1]=='': del items[-1]
            time_index=int(items.pop(0))
            time=float(items.pop(0))
            self.times.append(time)
            nc1=self.num_columns+1
            nsets=len(items)/nc1
            for i in xrange(nsets):
                setvals=items[i*nc1:(i+1)*nc1]
                key=(int(setvals[0]),)
                vals=[fortran_float(val) for val in setvals[1:]]
                self.row_name.append(key+(time,))
                if not key in self.keys:
                    self._keyrows[key]=[]
                    self.keys.append(key)
                self._keyrows[key].append(self._rowindex)
                self._data.append(vals)
                self._rowindex+=1
        
    def read_data_TOUGH2_MP(self,configured):
        """Reads in the data, for TOUGH2_MP output."""
        def get_key(line):
            return tuple([fix_blockname(line[self.key_start[i]:self.key_start[i+1]].rstrip()) for i in xrange(self._nkeys)])
        def get_time(line): return fortran_float(line[self.time_pos[0]:self.time_pos[1]])
        def get_vals(line):
            start=self.col_start+[len(line)] # allow for lines of different lengths
            return [fortran_float(line[start[i]:start[i+1]]) for i in xrange(self.num_columns)]
        lines=self._file.readlines()
        last_time=None
        from copy import copy
        otherfile_keys=copy(self.keys)
        for line in lines:
            if line.strip():
                time=get_time(line)
                key=get_key(line)
                vals=get_vals(line)
                if (not configured) and (time<>last_time): self.times.append(time)
                last_time=time
                if not (key in otherfile_keys):
                    self.row_name.append(key+(time,))
                    if not key in self.keys:
                        self._keyrows[key]=[]
                        self.keys.append(key)
                    self._keyrows[key].append(self._rowindex)
                    self._data.append(vals)
                    self._rowindex+=1

    def read_data_TOUGHplus(self,configured):
        """Reads in the data, for TOUGH+ output."""
        lines=self._file.readlines()
        if self.type=='FOFT': 
            for line in lines:
                items=line.strip().split(',')
                time=fortran_float(items[1])
                self.times.append(time)
                nc1=self.num_columns+1
                nsets=(len(items)-2)/nc1
                for i in xrange(nsets):
                    setvals=items[2+i*nc1:2+(i+1)*nc1]
                    key=(int(setvals[0]),)
                    if key[0]>0:
                        vals=[fortran_float(val) for val in setvals[1:]]
                        self.row_name.append(key+(time,))
                        if not key in self.keys:
                            self._keyrows[key]=[]
                            self.keys.append(key)
                        self._keyrows[key].append(self._rowindex)
                        self._data.append(vals)
                        self._rowindex+=1
        else:
            for line in lines:
                vals=[fortran_float(val) for val in line.strip().split()]
                time=vals.pop(0)
                self.times.append(time)
                self._data.append(vals)

    def finalize_data(self):
        self._data=np.array(self._data,float64)
        self.times=np.array(self.times,float64)
        self._row=dict([(r,i) for i,r in enumerate(self.row_name)])

