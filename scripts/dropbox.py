#!/usr/bin/env python
#----------------------------------------------------------------------
#
# Name: dropbox.py
#
# Purpose: Executable script to return file transfer service dropbox
#          directory for a given file.  The same function is available
#          in python via python module uboone_utilities.
#
# Created: 28-Oct-2013  H. Greenlee
#
# Usage:
#
# dropbox.py <filename>
#
#----------------------------------------------------------------------

import sys
import uboone_utilities

if len(sys.argv) <= 1:
    print 'Usage: dropbox.py <filename>'
    sys.exit(0)

filename = sys.argv[1]
dir = uboone_utilities.get_dropbox(filename)
print dir
sys.exit(0)
