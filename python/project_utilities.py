#!/usr/bin/env python
#----------------------------------------------------------------------
#
# Name: project_utilities.py
#
# Purpose: A python module containing various python utility functions
#          and classes used by project.py and other python scripts.
#
# Created: 28-Oct-2013  H. Greenlee
#
#----------------------------------------------------------------------

import sys, os, stat, time, types
import socket
import subprocess
import getpass
import ROOT

proxy_ok = False


# Function to return the current experiment.
# The following places for obtaining this information are
# tried (in order):
#
# 1.  Environment variable $EXPERIMENT.
# 2.  Environment variable $SAM_EXPERIMENT.
# 3.  Hostname (up to "gpvm").
#
# Raise an exception if none of the above methods works.
#

def get_experiment():

    exp = ''
    for ev in ('EXPERIMENT', 'SAM_EXPERIMENT'):
        if os.environ.has_key(ev):
            exp = os.environ[ev]
            break

    if not exp:
        hostname = socket.gethostname()
        n = hostname.find('gpvm')
        if n > 0:
            exp = hostname[:n]

    if not exp:
        raise RuntimeError, 'Unable to determine experiment.'

    return exp

# Function to return the fictitious disk server node
# name used by sam for bluearc disks.

def get_bluearc_server():
    return get_experiment() + 'data:'

# Function to return the fictitious disk server node
# name used by sam for dCache disks.

def get_dcache_server():
    return 'fnal-dcache:'

# Function to determine dropbox directory based on sam metadata.
# Raise an exception if the specified file doesn't have metadata.
# This function should be overridden in <experiment>-utilities module.

def get_dropbox(filename):
    raise RuntimeError, 'Function get_dropbox not implemented.'

# Get role (normally 'Analysis' or 'Production').

def get_role():

    # If environment variable ROLE is defined, use that.  Otherwise, make
    # an educated guess based on user name.

    result = 'Analysis'   # Default role.

    # Check environment variable $ROLE.

    if os.environ.has_key('ROLE'):
        result = os.environ['ROLE']

    # Otherwise, check user.

    else:
        prouser = get_experiment() + 'pro'
        user = getpass.getuser()
        if user == prouser:
            result = 'Production'

    return result

# Get proxy path.

def get_proxy():
    exp = get_experiment()
    role = get_role()
    user = getpass.getuser()
    if role == 'Production':
        result = '/scratch/%s/grid/%s.%s.Production.proxy' % (user, user, exp)
    else:
        result = '/scratch/%s/grid/%s.%s.proxy' % (user, user, exp)
    return result

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
        print 'Waiting ...'

        # Reading the parent directory seems to make files be visible faster.

        os.listdir(os.path.dirname(path))
        time.sleep(1)
        ntry = ntry - 1

    # Timed out.

    return 1

# Method to optionally make a copy of a pnfs file.

def path_to_local(path):

    # Depending on the input path and the environment, this method
    # will do one of the following things.
    #
    # 1.  If the input path is a pnfs path (starts with "/pnfs/"), and
    #     if $TMPDIR is defined and is accessible, the pnfs file will
    #     be copied to $TMPDIR using ifdh, and the path of the local
    #     copy will be returned.
    #
    # 2.  If the input path is a pnfs path, and if $TMPDIR is not
    #     defined, is not accessible, or if the ifdh copy fails,
    #     this method will return the empty string ("").
    #
    # 3.  If the input path is anything except a pnfs path, this
    #     method will not do any copy and will return the input path.
    #     

    global proxy_ok
    result = ''

    # Is this a pnfs path?

    if path[0:6] == '/pnfs/':

        # Is there a temp directory?

        local = ''
        if os.environ.has_key('TMPDIR'):
            tmpdir = os.environ['TMPDIR']
            mode = os.stat(tmpdir).st_mode
            if stat.S_ISDIR(mode) and os.access(tmpdir, os.W_OK):
                local = os.path.join(tmpdir, os.path.basename(path))

        if local != '':

            # Do local copy.

            if not proxy_ok:
                proxy_ok = test_proxy()

            # Make sure local path doesn't already exist (ifdh cp may fail).

            if os.path.exists(local):
                os.remove(local)

            # Use ifdh to make local copy of file.

            #print 'Copying %s to %s.' % (path, local)
            rc = subprocess.call(['ifdh', 'cp', path, local])
            if rc == 0:
                rc = wait_for_stat(local)
                if rc == 0:

                    # Copy succeeded.

                    result = local

    else:

        # Not a pnfs path.

        result = path

    return result


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
# This class has the following data member.
#
# root_tfile - a ROOT.TFile object.
#
# This class aggregates rather than inherits from ROOT.TFile because the owned
# TFile can itself be polymorphic.
#
#

class SafeTFile:

    # Default constructor.

    def __init__(self):
        self.root_tfile = None

    # Initializing constructor.

    def __init__(self, path):
        self.Open(path)

    # Destructor.

    def __del__(self):
        self.Close()

    # Unbound (static) Open method.

    def Open(path):
        return SafeTFile(path)

    # Bound Open method.

    def Open(self, path):

        self.root_tfile = None

        # Open file, with special handling for pnfs paths.

        local = path_to_local(path)
        if local != '':

            # Open file or local copy.

            self.root_tfile = ROOT.TFile.Open(local)

            # Now that the local copy is open, we can safely delete it already.

            if local != path:
                os.remove(local)

        else:

            # Input path is pnfs, but we could not get a local copy.
            # Get xrootd url instead.a

            global proxy_ok
            if not proxy_ok:
                proxy_ok = test_proxy()
            url = path_to_url(path)
            #print 'Using url %s' % url
            self.root_tfile = ROOT.TFile.Open(url)

    # Close method.

    def Close(self):

        # Close file and delete temporary file (if any and not already deleted).

        if self.root_tfile != None and self.root_tfile.IsOpen():
            self.root_tfile.Close()
            self.root_tfile = None

    # Copies of regular TFile methods used in project.py.

    def IsOpen(self):
        return self.root_tfile.IsOpen()

    def IsZombie(self):
        return self.root_tfile.IsZombie()

    def GetListOfKeys(self):
        return self.root_tfile.GetListOfKeys()

    def Get(self, objname):
        return self.root_tfile.Get(objname)

# Import experiment-specific utilities.  In this imported module, one can 
# override any function or symbol defined above, or add new ones.

from experiment_utilities import *
