#!/usr/bin/env python

# Import stuff.

import sys, os, string, subprocess, json, project_utilities

# Import ROOT (hide command line arguments).

myargv = sys.argv
sys.argv = myargv[0:1]
sys.argv.append('-n')
os.environ['TERM'] = 'vt100'    # Prevent root from printing garbage on initialization.
import ROOT
ROOT.gErrorIgnoreLevel = ROOT.kError
sys.argv = myargv

# Convert adler32-1 (used by dcache) to adler32-0 (used by sam).

def convert_1_adler32_to_0_adler32(crc, filesize):
    crc = long(crc)
    filesize = long(filesize)
    size = int(filesize % 65521)
    s1 = (crc & 0xffff)
    s2 = ((crc >> 16) &  0xffff)
    s1 = (s1 + 65521 - 1) % 65521
    s2 = (s2 + 65521 - size) % 65521
    return (s2 << 16) + s1


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

    crc = {}
    srm_url = project_utilities.path_to_srm_url(path)

    if srm_url == path:
        try:
            f = open(path,'rb')
            crc = enstoreChecksum(f)
        except (IOError, OSError), ex:
            raise Error(str(ex))
        finally:
            f.close()
    else:
        try:
            # Following commented commands are old way of calculating checksum by
            # transferring entire file over network.
            # Should work again if uncommented (if srm way breaks).

            #cmd = ['ifdh', 'cp', path, '/dev/fd/1']
            #p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
            #f = p.stdout
            #crc = enstoreChecksum(f)

            # New (clever, efficient, obscure...) way of accessing dCache 
            # stored checksum using srm.
            cmd = ['srmls', '-2', '-l', srm_url]
            srmout = subprocess.check_output(cmd)
            first = True
            crc0 = 0
            for line in string.split(srmout, '\n'):
                if first:
                    size = long(line[2:line.find('/')-1])
                    first = False
                    continue
                if line.find("Checksum value:") > 0:
                    ssum = line[line.find(':') + 2:]
                    crc1 = long( ssum , base = 16 )
                    crc0 = convert_1_adler32_to_0_adler32(crc1, size)
                    break

            crc = {"crc_value": str(crc0), "crc_type": "adler 32 crc type"}
            
        except (IOError, OSError), ex:
            raise Error(str(ex))

        return crc

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

        file = project_utilities.SafeTFile(inputfile)
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
