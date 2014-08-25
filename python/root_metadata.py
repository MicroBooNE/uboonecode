#!/usr/bin/env python

# Import stuff.

import sys, os, json

# Import ROOT (hide command line arguments).

myargv = sys.argv
sys.argv = myargv[0:1]
sys.argv.append('-n')
os.environ['TERM'] = 'vt100'    # Prevent root from printing garbage on initialization.
import ROOT
ROOT.gErrorIgnoreLevel = ROOT.kError
sys.argv = myargv

# Checksum utilities copied from sam_web_client

def enstoreChecksum(fileobj):
    import zlib
    readblocksize = 1024*1024
    crc = 0
    while 1:
        try:
            s = fileobj.read(readblocksize)
        except (OSError, IOError), ex:
            raise Error(str(ex))
        if not s: break
        crc = zlib.adler32(s,crc)
    crc = long(crc)
    if crc < 0:
        # Return 32 bit unsigned value
        crc  = (crc & 0x7FFFFFFFL) | 0x80000000L
    return { "crc_value" : str(crc), "crc_type" : "adler 32 crc type" }

def fileEnstoreChecksum(path):
    """Calculate enstore compatible CRC value"""
    try:
        f =open(path,'rb')
    except (IOError, OSError), ex:
        raise Error(str(ex))
    try:
        return enstoreChecksum(f)
    finally:
        f.close()

def get_external_metadata(inputfile):

	# define an empty python dictionary
	md = {}

        # Check whether this file exists.
        if not os.path.exists(inputfile):
            return md
            
	# Get the other meta data field parameters						
	md['file_name'] =  os.path.basename(inputfile)
	md['file_size'] =  str(os.path.getsize(inputfile))
	md['crc'] = fileEnstoreChecksum(inputfile)

	# Root checks.

        file = ROOT.TFile.Open(inputfile)
        if file and file.IsOpen() and not file.IsZombie():

            # Root file opened successfully.
            
            obj = file.Get('Events')
            if obj and obj.InheritsFrom('TTree'):

                # This has a TTree names Events.

                nev = obj.GetEntriesFast()
                md['events'] = str(nev)
        else:

            # Root file could not be opened.

            md = {}
	
	return md

if __name__ == "__main__":
	md = get_external_metadata(str(sys.argv[1]))
	#print md	
	mdtext = json.dumps(md, sys.stdout, indent=2, sort_keys=True)
	print mdtext
	sys.exit(0)	
