#!/bin/env python

import os
import sys
import time
from datetime import datetime
from IOVAPI import IOVDB

class dbfolder():

    def __init__(self, host, port, dbname, user_name, folder_name ):

        if folder_name == None:
            print "\nNo folder given. That's bad.  Seek help (-h)"
            sys.exit(0)
            
        print "\nConnecting to DB Folder:"
        print "   host: ", host
	print "   port: ", port
        print "     db: ", dbname
        print " folder: ", folder_name
        print "   user: ", user_name

        self.db = IOVDB(connstr = "dbname="+dbname+" user="+user_name+" host="+host+" port="+port )
        
        self.dbname = dbname
        self.folder_name = folder_name
        #self.io_class = io_class
            
      
    ############################################################################################
    #
    #  Read data from this DB at time, t_query
    #
    ############################################################################################                 
    def read(self, t_query = time.time(), channelid = -1, tag = None):

        print "    tag: ", tag
        if t_query == -1: t_query = time.time()
        t_q = datetime.fromtimestamp(t_query)
        print '\nQuery time: ', t_q
        print '  in seconds: ', t_query

        try:
            folder = self.db.openFolder(self.folder_name)
        except:
            print '\nCould not open folder', self.folder_name, '!'
            sys.exit(2)
            
	columns = folder.Columns
        data = folder.getData(t_query,tag=tag)

        col_width = 8
        for c in columns:
	    if len(str(c)) >= col_width:
	        col_width = len(str(c))+1

        try:

            #print IOV
	    print '\nIOV for requested data: ', data.iov()[0],' - ', data.iov()[1],'\n'
	    
	    #print folder header
	    print 'Channel',
	    for sp in range(0,col_width-7):
	        print '',	    
	    for c in columns: 
	        print c,
		for sp in range(0,col_width-len(str(c))):
		    print '',
		    
		    
            #print folder contents 
	    data = data.asDict()
	    ctr = -1
            for chid, v in data.iteritems():
                ctr = ctr+1
                if ctr==10 and channelid==-1:
		    msg = str(raw_input('\nYou have elected to view all '+str(len(data))+' channels!  Are you sure you want to do this? [y/n]'))
		    if msg != 'y' and msg != 'yes':
		        sys.exit(0)

                if channelid == -1 or chid == channelid:
		    print '\n', chid,
		    for sp in range(0,col_width-len(str(chid))):
		        print '',
                    for c in range(len(v)):
                        print v[c], 
			for sp in range(0,col_width-len(str(v[c]))):
			    print '',
                
        except:
            print '\nDATA NOT FOUND'
		    
            

    ############################################################################################
    #
    #  Write data to this DB
    #
    ############################################################################################
    def write(self, file_name = None, override_future = 'no', t_start = 0):
        if file_name == None:
            print '\nNice try, buddy, you need a file to upload data!'
            sys.exit(2)

        if os.path.exists(file_name) != True:
            print '\nFile', file_name, 'not found!'
            return 1   # file not found

        t0, data, columns = self.parseFile(file_name)

        if len(data) == 0:
            return 2   # no good data in the file

        # if start time not specified, use the time from the file
        if t_start == 0: t_start = t0

        #if still not specified, complain
        if t_start == 0:
            print '\nFile ', file_name, 'does not include a time stamp! ***DATA NOT INSERTED***\n'
            return 3   # IOV time not specified

        #open folder
        try:
            folder = self.db.openFolder(self.folder_name)
        except:
            print '\nFolder ', self.folder_name, 'does not exist!'
            return 4   # Could not open DB folder

        #make sure folder and data agree on number of columns
        folder_columns = folder.Columns
	if len(folder_columns) != len(columns):
	    print '\nFile ', file_name, 'does not have same number of columns as destination folder! ***DATA NOT INSERTED***\n'
	    return 5  #column mismatch
	    

        print '\nStoring data with starting t(%s) = %s' % (t_start, datetime.fromtimestamp(t_start))
        T0 = time.time()
        folder.addData(t_start, data, True, override_future)
        print '\nData stored. Elapsed time = %.1f seconds' % (time.time() - T0,)

        return 0   # no errors, write successful

    ############################################################################################
    #
    #  Create new folder in DB
    #
    ############################################################################################
    def create(self, columnNames, columnTypes):
                   
        columnNamesTypes = []
	for n,t in zip(columnNames,columnTypes):
	   columnNamesTypes.append( (n,t) )
 
        print '\nCreating folder %s with %d columns:\n   %s' % (self.folder_name, len(columnNames), ', '.join(columnNames))

        self.db.createFolder(self.folder_name, columnNamesTypes, False, grants={'uboone_web':'r'})    
    
    ############################################################################################
    #
    #  Tag current state of folder in DB
    #
    ############################################################################################
    def tag(self, tag = None, message = ''):
    
        if tag == None:
	    print '\nNeed to specify tag'
	    return 6
	
	#open folder
        try:
            folder = self.db.openFolder(self.folder_name)
        except:
            print '\nFolder ', self.folder_name, 'does not exist!'
            return 4   # Could not open DB folder
	    
	print '\nApplying tag: '+tag+' to current state of folder: '+self.folder_name
	folder.createTag(tag, message)
    
    ############################################################################################
    #
    #  parse file
    #
    ############################################################################################
    def parseFile(self, fname):

	t0 = 0
	column_names = []
	data = {}
   
        try:
	    content = open(fname,'r').read()
	except:
	    print '\nFailed to open file ',fname,'!'
	    return t0, data, column_names

	content = content.split('\n')
	for line in content:
	    if not line: continue
	    words = line.split()
	    if words[0][0] == '#':
	        if words[1].lower() == 'time':
		    t0 = float(words[2])
		elif words[1].lower() == 'channel':
		    column_names = words[2:]   #expects column names to be last line of header!
		else:
		    print '\nInvalid header in file ',fname,'!'
		    return t0, data, column_names
            else:
	        line = line.strip()
		if line:
		    chid = int(words[0])
		    data[chid] = words[1:]
		      
	
	print '\n%s parsed, %d entries, %d columns, at time %f' % (fname, len(data), len(column_names), t0)
	assert column_names != [] and column_names != None
	
	return t0, data, column_names
