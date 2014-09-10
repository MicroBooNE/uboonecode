#!/usr/bin/env python
#----------------------------------------------------------------------
#
# Name: dropbox.py
#
# Purpose: Executable script to return file transfer service dropbox
#          directory for a given file.  The same function is available
#          in python via python module project_utilities.
#
# Created: 28-Oct-2013  H. Greenlee
#
# Usage:
#
# dropbox.py <filename>
#
#----------------------------------------------------------------------

import sys
import project_utilities

if len(sys.argv) <= 1:
    print 'Usage: dropbox.py <filename>'
    sys.exit(0)

filename = sys.argv[1]
dir = project_utilities.get_dropbox(filename)
print dir
sys.exit(0)
