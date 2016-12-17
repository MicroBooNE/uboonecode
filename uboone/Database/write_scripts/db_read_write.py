import sys
import os
import getopt
import __main__
from dbfolder import dbfolder

def help_me():

  print """
  db_read_write.py [OPTIONS]
  
  Author: Brandon Eberly (eberly@slac.stanford.edu)
  (Based on scripts written by Kazu and work done by the MINERvA collaboration)
  
  
  ----------------------------------------------------
  OPTIONS
  ----------------------------------------------------
  -h,--help     : print help

  -a,--action   : do an action [create, write, read, tag] (required)

  -f,--folder   : folder name (required)
  
  -t,--time     : time (seconds) to read/write from database
                  writing: retrieves from file by default
	        	   Can't be used with --is_playlist
	          reading: defaults to NOW
	          
  -i,--input    : writing: single file or playlist of files to load into DB
                  creating folder: file with lines in format: column_name type 
  
  --is_playlist : use in conjunction with -i,--input if file is a playlist of files
  
  -c,--channel  : specify channel when reading from DB (default is all channels)
  
  --override    : causes new data IOV to extend to infinity, even if there
                  are later IOVs in the DB.  Use with caution!  Options are 'hide' 
		  and 'delete'
	          
  --production  : use the production database (default is development)
  
  --tag         : database tag
  
  ----------------------------------------------------
  NOTE ON TIMES
  ----------------------------------------------------
  This routine expects start and end times in integer seconds since the 
  start of the unix epoch in UTC.
  
  
  ----------------------------------------------------
  NOTE ON IOVS AND OVERLAP WHEN WRITING DATA
  ----------------------------------------------------  
  By default, if the timestamp t is later than the start time of all other 
  IOVs in the folder, then the written data is given an IOV that starts at 
  t and has no end time.  The other IOVs are adjusted accordingly.
  
  If the timestamp t is earlier than the start time of at least one IOV, then 
  the written data is given an IOV that starts at t and ends at the begin_time 
  of the next IOV chronologically.  This masks part or all of the IOV that 
  previously preceded the next IOV.  This behavior can be changed using the 
  --override option, so that the written data has an IOV that starts at t and 
  has no end time.  This will mask all other IOVs that have a start time >= t.
  
  
  ----------------------------------------------------
  File Format
  ---------------------------------------------------- 
  Files that are written to the database via the -i(--input) option must be 
  formatted as follows (N columns, M channels, <> indicate substitution required):
  
  # time <t>
  # channel <column name 1> <column name 2> ... <column name N>
  <chan# 1> <number 1,1>    <number 1,2>    ... <number 1,N>
  <chan# 2> <number 2,1>    <number 2,2>    ... <number 2,N>
  ...
  <chan# M> <number M,1>    <number M,2>    ... <number M,N>
  
  Note that <t> is the start of the interval of validity for the file 
  contents in time since the start of the unix epoch.  The column names 
  should match those of the folder that is being written to (specified by 
  -f(--folder)). Also, channel numbers are not required to start at 1 or 
  be consecutive.  If a column is an array, the array contents should be contained 
  in {} and separated by commas without whitespace.
  
  
  ----------------------------------------------------
  Examples:
  ----------------------------------------------------
  
  #creates a folder called detpedestals in the production database
  python db_read_write.py --action create --folder detpedestals --production -i format_file.txt
  
  #tags a folder called detpedestals as 'v1' in the development database
  python db_read_write.py --action tag --folder detpedestals --tag v1
  
  
  #prints the contents of folder pmtgains_data in the development database 
  #valid at the current time  
  python db_read_write.py --action read --folder pmtgains_data

  
  #prints the contents of folder pmtgains_data in the development database 
  #for each interval of validity between the start of the unix epoch and 
  #the current time  
  python db_read_write.py --action read --folder pmtgains_data -t 0
  
  
  #writes the contents of the files in playlist.dat to the larproperties 
  #folder in the development database.  Uses timestamps included in the files
  python db_read_write.py -a write -f larproperties -i playlist.dat --is_playlist
  
  #writes the content of file.txt to the larproperties folder in the 
  #development database.  Start of interval of validity is 1440000000, end is 
  #infinity (unless there is a later IOV in this folder).
  python db_read_write.py -a write -f larproperties -i file.txt -t 1440000000

  """
  
def main():

  try:
    opts, args = getopt.gnu_getopt(sys.argv[1:],
                                   "ha:f:t:i:c:",
				   ["help",
				    "action=",
				    "folder=",
				    "time=",
				    "input=",
				    "tag=",
				    "channel=",
				    "production",
				    "override=",
				    "is_playlist"])
  except getopt.GetoptError, err:
    print str(err)
    help_me()
    sys.exit(2)
    
  #-----------------------------------------------------
  #
  # Default parameter values
  #
  #-----------------------------------------------------
  
  host   = 'ifdb04.fnal.gov'
  port   = '5437'
  dbname = 'microboone_dev'
  user   = os.environ['USER']
  
  folder   = None
  action   = None
  filename = None
  tag      = None
  channel  = -1
  
  override_future = 'no'
  is_playlist     = False
  
  start_sec = 0


  #-----------------------------------------------------
  #
  # Read arguments
  #
  #-----------------------------------------------------
  
  if len(opts) == 0:
    help_me()
    sys.exit(0)
    
  for o, a in opts:
    if o in ("-h", "--help"):
      help_me()
      sys.exit(0)
    elif o in ("-a", "--action"):
      action = a
    elif o in ("-f", "--folder"):
      folder = a
    elif o in ("-i", "--input"):
      filename = a
    elif o in ("-t", "--time"):
      start_sec=float(a)
    elif o in ("--tag"):
      tag = a
    elif o in ("-c", "--channel"):
      channel = float(a)
    elif o in ("--override"):
      override_future = a
    elif o in ("--is_playlist"):
      is_playlist = True
    elif o in ("--production"):
      host = "ifdb05.fnal.gov"
      dbname = "microboone_prod"
    else:
      assert False, "unhandled option "+o
      
  if action == None or folder == None:
    print "Need to specify action and folder\n"
    sys.exit(0)
  
  if action == 'write' and is_playlist and start_sec != 0:
    print "Cannot specify --is_playlist and -t,--time when writing!\n"
    sys.exit(0)
  
  #-----------------------------------------------------
  #
  # Do requested action
  #
  #-----------------------------------------------------   
  
  if   action == 'read':
    dbf = dbfolder(host, port, dbname, user, folder)
    dbf.read(start_sec, channel, tag)
    
  elif action == 'write':  

    if filename == None:
      print "\nWhen writing, must specify a file with -i,--input\n"
      sys.exit(0)

    dbf = dbfolder(host, port, dbname, user, folder)

    # if a single filename was passed in
    if is_playlist == False:
      print override_future
      dbf.write(filename, override_future, start_sec)

    # read in playlist one file at a time
    else:    
      with open(filename,'r') as f:
        for line in f:
	  if not line: continue
	  line2 = line.strip('\n')
	  dbf.write(line2, override_future, 0)
	  
  elif action == 'create':
    if filename == None:
      print "\nWhen creating a folder, must specify a file with -i,--input\n"
      sys.exit(0)
      
    column_names = []
    column_types = []
    with open(filename) as f:
      for line in f:
        words = line.split()
	column_names.append(words[0])
	column_types.append(words[1])
	
    dbf = dbfolder(host, port, dbname, user, folder)
    dbf.create(column_names, column_types)  
  elif action == 'tag':
    dbf = dbfolder(host, port, dbname, user, folder)
    dbf.tag(tag,'')
  else:
    print "Invalid action: "+action+".  Seek help (-h) and try again\n"
    sys.exit(0)

if __name__ == '__main__':
  main()
