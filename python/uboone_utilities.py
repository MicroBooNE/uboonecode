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

import sys, os

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
