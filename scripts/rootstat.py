#! /usr/bin/env python
######################################################################
#
# Name: stat.py
#
# Purpose: Analyze art root file and dump object statistics.
#
# Created: 27-Nov-2012  Herbert Greenlee
#
# Usage:
#
# stat.py <options> [@filelist] [file1 file2 ...]
#
# Options:
#
# [-h|--help] - Print help message.
# --level n   - Branch level (default 1).  Use --level 1 to see top
#               branches only.  Use --level 2 to also see subbranches.
# --nfile n   - Number of files to analyze (default all).
# --all       - Print analysis of each file (default is only summary).
#
# Arguments:
#
# @filelist       - File list containing one input file per line.
# file1 file2 ... - Input files.
#
######################################################################

import sys, os, string
import project_utilities

# Import ROOT module.
# Globally turn off root warnings.
# Don't let root see our command line options.

myargv = sys.argv
sys.argv = myargv[0:1]
import ROOT
ROOT.gErrorIgnoreLevel = ROOT.kError
sys.argv = myargv

# Print help.

def help():

    filename = sys.argv[0]
    file = open(filename)

    doprint=0
    
    for line in file.readlines():
        if line[2:9] == 'stat.py':
            doprint = 1
        elif line[0:6] == '######' and doprint:
            doprint = 0
        if doprint:
            if len(line) > 2:
                print line[2:],
            else:
                print

# Analyze root file.

def analyze(root, level, gtrees, gbranches, doprint):

    trees = {}
    events = None
    keys = root.GetListOfKeys()
    for key in keys:
        objname = key.GetName()
        if not trees.has_key(objname):
            obj = root.Get(objname)
            if obj and obj.InheritsFrom('TTree'):
                trees[objname] = obj
                if objname == 'Events':
                    events = obj

    # Print summary of trees.

    if doprint:
        print '\nTrees:\n'
    for key in sorted(trees.keys()):
        tree = trees[key]
        nentry = tree.GetEntriesFast()
        if doprint:
            print '%s has %d entries.' % (key, nentry)

        # Remember information about trees.

        if gtrees.has_key(key):
            gtrees[key] = gtrees[key] + nentry
        else:
            gtrees[key] = nentry

    # Print summary of branches in Events tree.

    if doprint:
        print '\nBranches of Events tree:\n'

    # If level is zero, we are done (don't analyze branches).

    if level == 0:
        return

    if events:

        if doprint:
            print '   Total bytes  Zipped bytes   Comp.  Branch name'
            print '   -----------  ------------   -----  -----------'
            
        branches = events.GetListOfBranches()
        ntotall = 0
        nzipall = 0

        # Loop over branche of Events tree.

        for branch in sorted(branches):
            branch_class = branch.GetClass().GetName()

            # Only look at data products (class art::Wrapper<T>).
            
            if branch_class[0: 13] == 'art::Wrapper<':

                # Loop over subbranches.
                
                subbranches = branch.GetListOfBranches()
                for subbranch in sorted(subbranches):
                    name = subbranch.GetName()

                    # Only look at '.obj' subbranch (wrapped object).
                    
                    if name[-4:] == '.obj':
                        ntot = subbranch.GetTotBytes("*")
                        nzip = subbranch.GetZipBytes("*")
                        ntotall = ntotall + ntot
                        nzipall = nzipall + nzip
                        if doprint:
                            comp = float(ntot) / float(nzip)
                            print '%14d%14d%8.2f  %s' % (ntot, nzip, comp, name)

                        # Remember information about branches.
                        
                        if gbranches.has_key(name):
                            gbranches[name][0] = gbranches[name][0] + ntot
                            gbranches[name][1] = gbranches[name][1] + nzip
                        else:
                            gbranches[name] = [ntot, nzip]

                        # Loop over subsubbranches (attributes of wrapped object).
                        
                        if level > 1:
                            subsubbranches = subbranch.GetListOfBranches()
                            for subsubbranch in sorted(subsubbranches):
                                name = subsubbranch.GetName()
                                ntot = subsubbranch.GetTotBytes("*")
                                nzip = subsubbranch.GetZipBytes("*")
                                if doprint:
                                    comp = float(ntot) / float(nzip)
                                    print '%14d%14d%8.2f  %s' % (ntot, nzip, comp,
                                                                 subsubbranch.GetName())

                                # Remember information about branches.
                        
                                if gbranches.has_key(name):
                                    gbranches[name][0] = gbranches[name][0] + ntot
                                    gbranches[name][1] = gbranches[name][1] + nzip
                                else:
                                    gbranches[name] = [ntot, nzip]

        # Do summary of all branches.

        name = 'All branches'
        if doprint:
            comp = float(ntotall) / float(nzipall)
            print '%14d%14d%8.2f  %s' % (ntotall, nzipall, comp, name)

            # Print average event size.

            nev = events.GetEntriesFast()
            nevtot = 1.e-6 * float(ntotall) / float(nev)
            nevzip = 1.e-6 * float(nzipall) / float(nev)
            print
            print '%10d events.' % nev
            print '%7.2f Mb average size per event.' % nevtot
            print '%7.2f Mb average zipped size per event.' % nevzip

        if gbranches.has_key(name):
            gbranches[name][0] = gbranches[name][0] + ntotall
            gbranches[name][1] = gbranches[name][1] + nzipall
        else:
            gbranches[name] = [ntotall, nzipall]


    # Done.                     
    
    return
                
# Main program.

def main(argv):

    # Parse arguments.

    input_files = []
    level = 1
    nfilemax = 0
    all = 0

    args = argv[1:]
    while len(args) > 0:
        if args[0] == '-h' or args[0] == '--help':

            # Help.
            
            help()
            return 0

        elif args[0] == '--level' and len(args) > 1:

            # Analyze level.

            level = int(args[1])
            del args[0:2]
            
        elif args[0] == '--nfile' and len(args) > 1:

            # Number of files.

            nfilemax = int(args[1])
            del args[0:2]
            
        elif args[0] == '--all':

            # All files flag.

            all = 1
            del args[0]
            
        elif args[0][0] == '-':

            # Unknown option.

            print 'Unknown option %s' % args[0]
            return 1
            
        elif args[0][0] == '@':

            # Read in file list to input files.
            
            filelistname = args[0][1:]
            if project_utilities.safeexist(filelistname):
                for filename in project_utilities.saferead(filelistname):
                    input_files.append(string.strip(filename))
            else:
                print 'File list %s does not exist.' % filelistname
                return 1
            del args[0]
        else:

            # Add single file to input files.
            
            input_files.append(args[0])
            del args[0]

    # Loop over input files.

    gtrees = {}
    gbranches = {}
    nfile = 0

    for input_file in input_files:

        if nfilemax > 0 and nfile >= nfilemax:
            break
        nfile = nfile + 1

        if not project_utilities.safeexist(input_file):
            print 'Input file %s does not exist.' % input_file
            return 1

        print '\nOpening %s' % input_file
        root = project_utilities.SafeTFile(input_file)
        if not root.IsOpen() or root.IsZombie():
            print 'Failed to open %s' % input_file
            return 1

        # Analyze this file.
        
        analyze(root, level, gtrees, gbranches, all)

    print '\n%d files analyzed.' % nfile
                    
    # Print summary of trees.

    print '\nTrees from all files:\n'
    for key in sorted(gtrees.keys()):
        nentry = gtrees[key]
        print '%s has %d total entries.' % (key, nentry)

    # Print summary of branches.

    if level > 0:
        print '\nBranches of Events tree from all files:\n'
        print '   Total bytes  Zipped bytes   Comp.  Branch name'
        print '   -----------  ------------   -----  -----------'
    allname = 'All branches'
    ntot = 0
    nzip = 0
    for key in sorted(gbranches.keys()):
        if key != allname:
            ntot = gbranches[key][0]
            nzip = gbranches[key][1]
            comp = float(ntot) / float(nzip)
            print '%14d%14d%8.2f  %s' % (ntot, nzip, comp, key)
    if gbranches.has_key(allname):
        ntot = gbranches[allname][0]
        nzip = gbranches[allname][1]
        comp = float(ntot) / float(nzip)
        print '%14d%14d%8.2f  %s' % (ntot, nzip, comp, allname)

    # Print average event size.

    if gtrees.has_key('Events'):
        nev = gtrees['Events']
        nevtot = 1.e-6 * float(ntot) / float(nev)
        nevzip = 1.e-6 * float(nzip) / float(nev)
        print
        print '%10d events.' % nev
        if level > 0:
            print '%7.2f Mb average size per event.' % nevtot
            print '%7.2f Mb average zipped size per event.' % nevzip
    

    # Done.

    return 0

# Invoke main program.

if __name__ == '__main__':
    sys.exit(main(sys.argv))
