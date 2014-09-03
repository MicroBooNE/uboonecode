#!/usr/bin/env python
#----------------------------------------------------------------------
#
# Name: uboone_utilities.py
#
# Purpose: A python module containing various experiment-specific
#          python utility functions.
#
# Created: 28-Oct-2013  H. Greenlee
#
#----------------------------------------------------------------------

import sys, os, stat, time
import subprocess
import ROOT

proxy_ok = False


# Don't fail if samweb is not available.

try:
    import samweb_cli
except ImportError:
    pass

# Function to return the current sam experiment.
# The following places for obtaining this information are
# tried (in order):
#
# 1.  Environment variable $SAM_EXPERIMENT.
# 2.  Environment variable $EXPERIMENT.
# 3.  Environment variable $GROUP.
#
# Raise an exception if none of the above environment
# variables is defined.
#

def get_experiment():

    exp = ''
    for ev in ('SAM_EXPERIMENT', 'EXPERIMENT', 'GROUP'):
        if os.environ.has_key(ev):
            exp = os.environ[ev]
            break

    if not exp:
        raise RuntimeError, 'Unable to determine experiment.'

    return exp

# Function to return the fictitious disk server node
# name to use for bluearc disks.

def get_bluearc_server():
    return get_experiment() + 'data:'

# Function to return the fictitious disk server node
# name to use for dCache disks.

def get_dcache_server():
    return 'fnal-dcache:'

# Function to determine dropbox directory based on sam metadata.
# Raise an exception if the specified file doesn't have metadata.

def get_dropbox(filename):

    # Get metadata.

    md = {}
    exp = get_experiment()
    samweb = samweb_cli.SAMWebClient(experiment=exp)
    try:
        md = samweb.getMetadata(filenameorid=filename)
    except:
        pass

    # Extract the metadata fields that we need.
    
    file_type = ''
    group = ''
    data_tier = ''

    if md.has_key('file_type'):
        file_type = md['file_type']
    if md.has_key('group'):
        group = md['group']
    if md.has_key('data_tier'):
        data_tier = md['data_tier']

    if not file_type or not group or not data_tier:
        raise RuntimeError, 'Missing or invalid metadata for file %s.' % filename

    # Construct dropbox path.

    path = '/uboone/data/uboonepro/dropbox/%s/%s/%s' % (file_type, group, data_tier)
    return path

# Function to optionally convert a filesystem path into an xrootd url.
# Only affects paths in /pnfs space.

def path_to_url(path):
    url = path
    if path[0:6] == '/pnfs/':
        url = 'root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/' + path[6:]
    return url

# Function to optionally convert a filesystem path into an srm url.
# Only affects paths in /pnfs space.

def path_to_srm_url(path):
    srm_url = path
    if path[0:6] == '/pnfs/':
        srm_url = 'srm://fndca1.fnal.gov:8443/srm/managerv2?SFN=/pnfs/fnal.gov/usr/' + path[6:]
    return srm_url

# dCache-safe method to test whether path exists without opening file.

def safeexist(path):
    try:
        os.stat(path)
        return True
    except:
        return False

# Test whether user has a valid grid proxy.  Exit if no.

def test_proxy():
    global proxy_ok
    if not proxy_ok:
        try:
            subprocess.check_call(['voms-proxy-info', '-exists'], stdout=-1)
            proxy_ok = True
        except:
            print 'Please get a grid proxy.'
            os._exit(1)
    return proxy_ok

# dCache-safe method to return contents (list of lines) of file.

def saferead(path):
    lines = []
    if path[0:6] == '/pnfs/':
        test_proxy()
        proc = subprocess.Popen(['ifdh', 'cp', path, '/dev/fd/1'], stdout=subprocess.PIPE)
        lines = proc.stdout.readlines()
    else:
        lines = open(path).readlines()
    return lines

# Like os.path.isdir, but faster by avoiding unnecessary i/o.

def fast_isdir(path):
    result = False
    if path[-5:] != '.list' and \
            path[-5:] != '.root' and \
            path[-4:] != '.txt' and \
            path[-4:] != '.fcl' and \
            path[-4:] != '.out' and \
            path[-4:] != '.err' and \
            path[-3:] != '.sh' and \
            path[-5:] != '.stat' and \
            os.path.isdir(path):
        result = True
    return result

# Wait for file to appear on local filesystem.

def wait_for_stat(path):

    ntry = 60
    while ntry > 0:
        if os.access(path, os.R_OK):
            return 0
        #print 'Waiting ...'
        time.sleep(1)
        ntry = ntry - 1

    # Timed out.

    return 1

# DCache-safe TFile-like class for opening files for reading.
#
# Class SafeTFile acts as follows.
#
# 1.  When initialized with a pnfs path (starts with "/pnfs/"), SafeTFile uses
#     one of the following methods to open the file.
#
#     a) Open as a regular file (posix open).
#
#     b) Convert the path to an xrootd url (xrootd open).
#
#     c) Copy the file to a local temp disk using ifdh (copy to $TMPDIR or
#        local directory) using ifdh, and open the local copy.
#
# 2.  When initialized with anything except a pnfs path (including regular
#     filesystem paths and urls), SafeTFile acts exactly like TFile.
#
# Implementation notes.
#
# This class has the following two data members.
#
# root_tfile - a ROOT.TFile object.
# local      - Name of local file (for ifdh access).
#
# This class aggregates rather than inherits from ROOT.TFile because the owned
# TFile can itself be polymorphic.
#
#

class SafeTFile:

    # Default constructor.

    def __init__(self):
        self.root_tfile = None
        self.local = ''

    # Initializing constructor.

    def __init__(self, path):
        self.root_tfile = None
        self.local = ''
        self.Open(path)

    # Destructor.

    def __del__(self):
        self.Close()

    # Unbound (static) Open method.

    def Open(path):
        return SafeTFile(path)

    # Bound Open method.

    def Open(self, path):

        global proxy_ok
        self.root_tfile = None
        self.local = ''

        # Open file, with special handling for pnfs paths.

        url = path_to_url(path)
        if url != path:
            if not proxy_ok:
                proxy_ok = test_proxy()

            # Decide whether to use xrootd or ifdh.
            # See if there is a writable directory $TMPDIR.  If there is,
            # construct a path for a local copy of the pnfs file.

            if os.environ.has_key('TMPDIR'):
                tmpdir = os.environ['TMPDIR']
                mode = os.stat(tmpdir).st_mode
                if stat.S_ISDIR(mode) and os.access(tmpdir, os.W_OK):
                    self.local = os.path.join(tmpdir, os.path.basename(path))

                    # Make sure local path doesn't already exist (ifdh cp may fail).

                    if os.path.exists(self.local):
                        os.remove(self.local)


            if self.local != '':

                # Use ifdh to make local copy of file.

                #print 'Copying %s to %s.' % (path, self.local)
                rc = subprocess.call(['ifdh', 'cp', path, self.local])
                if rc == 0:
                    rc = wait_for_stat(self.local)
                    if rc == 0:
                        self.root_tfile = ROOT.TFile.Open(self.local)

                        # Now that the local copy is open, we can safely delete it already.

                        os.remove(self.local)
                        self.local = ''

                if rc != 0:

                    # Ifdh copy did not succeed.

                    self.root_tfile = None
                    self.local = ''

            else:

                # No local path, use xrootd url.

                #print 'Using url %s' % url
                self.root_tfile = ROOT.TFile.Open(url)

        else:

            # Not a pnfs-path.  Use regular TFile open.

            self.root_tfile = ROOT.TFile.Open(path)
            self.local = ''

    # Close method.

    def Close(self):

        # Close file and delete temporary file (if any and not already deleted).

        if self.root_tfile != None and self.root_tfile.IsOpen():
            self.root_tfile.Close()
            self.root_tfile = None
        if self.local != '':
            os.remove(self.local)
            self.local = ''

    # Copies of regular TFile methods used in project.py.

    def IsOpen(self):
        return self.root_tfile.IsOpen()

    def IsZombie(self):
        return self.root_tfile.IsZombie()

    def GetListOfKeys(self):
        return self.root_tfile.GetListOfKeys()

    def Get(self, objname):
        return self.root_tfile.Get(objname)
