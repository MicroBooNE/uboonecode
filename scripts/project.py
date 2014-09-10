#! /usr/bin/env python
######################################################################
#
# Name: project.py
#
# Purpose: Production project script.
#
# Created: 11-Sep-2012  Herbert Greenlee
#
# Usage:
#
# project.py <options>
#
# Project options:
#
# --xml <-|file|url>  - Xml file containing project description.
# --project <project> - Project name (required if xml file contains
#                       more than one project description).
# --stage <stage>     - Project stage (required if project contains
#                       more than one stage).
# --histmerge <program> - Override default histogram merging program
#                         (no effect on histogram merging specified
#                         at stage level).
# --tmpdir <tempdir>  - Override TMPDIR internally.  If TMPDIR is set
#                       use ifdh cp instead of xrootd for accessing
#                       content of root files in dCache.
#
# Actions (specify one):
#
# [-h|--help]  - Print help (this message).
# [-xh|--xmlhelp] - Print xml help.
# --submit     - Submit all jobs for specified stage.
# --check      - Check results for specified stage and print message.
# --checkana   - Check analysis results for specified stage and print message.
# --status     - Print status of each stage.
# --makeup     - Submit makeup jobs for specified stage.
# --clean      - Delete output from specified project and stage.
# --declare    - Declare files to sam.
# --add_locations    - Check sam disk locations and add missing ones.
# --clean_locations  - Check sam disk locations and remove non-existent ones.
# --remove_locations - Remove all sam disk locations, whether or not file exists.
# --upload     - Upload files to enstore.
# --define     - Make sam dataset definition.
# --undefine   - Delete sam dataset definition.
#
# Information only actions:
#
# --outdir     - Print the name of the output directory for stage.
# --fcl        - Print the fcl file name and version for stage.
# --defname    - Print sam dataset definition name for stage.
# --input_files        - Print all input files.
# --check_declarations - Check whether data files are declared to sam.
# --test_declarations  - Print a summary of files returned by sam query.
# --check_locations    - Check sam locations and report the following:
#                        a) Files that lack any location.
#                        b) Disk locations that can be added.
#                        c) Incorrect disk locations that should be removed.
# --check_tape         - Check sam tape locations.
#                        Reports any files that lack tape (enstore) locations.
# --check_definition   - Reports whether the sam dataset definition associated
#                        with this project/stage exists, or needs to be created.
# --test_definition    - Print a summary of files returned by dataset definition.
#
######################################################################
#
# XML file structure
# ------------------
#
# The xml file must contain one or more elements with tag "project."
#
# The project element must have attribute "name."
#
# The following element tags withing the project element are recognized.
#
# <group>   - Group (default $GROUP).
#
# <numevents> - Total number of events (required).
# <numjobs> - Number of worker jobs (default 1).  This value can be
#             overridden for individual stages by <stage><numjobs>.
# <os>      - Specify batch OS (comma-separated list: SL5,SL6).
#             Default let jobsub decide.
#
# <script>  - Name of batch worker script (default condor_lar.sh).
#             The batch script must be on the execution path.
#
# <larsoft> - Information about larsoft release.
# <larsoft><tag> - Frozen release tag (default "development").
# <larsoft><qual> - Build qualifier (default "debug", or "prof").
# <larsoft><local> - Local test release directory or tarball (default none).
#
# <ubfcl>   - Ubfcl version (default none).
# <ubxml>   - Ubxml version (default none).
# <filetype> - Sam file type ("data" or "mc", default none).
# <runtype>  - Sam run type (normally "physics", default none).
# <histmerge> - Default histogram merging program (default "hadd -T", 
#               can be overridden at each stage).
#
# <stage name="stagename"> - Information about project stage.  There can
#             be multiple instances of this tag with different name
#             attributes.  The name attribute is optional if there is
#             only one project stage.
#
# <stage><fcl> - Name of fcl file (required).  Specify just the filename,
#             not the full path.
# <stage><outdir> - Output directory (required).  A subdirectory with the
#             project name is created underneath this directory.  Individual
#             workers create an additional subdirectory under that with
#             names like <cluster>_<process>.
# <stage><workdir> - Specify work directory (required).  This directory acts
#             as the submission directory for the batch job.  Fcl file, batch
#             script, and input file list are copied here.  A subdirectory with
#             the name of the project and "/work" are appended to this path.
#             This directory should be grid-accessible and located on an
#             executable filesystem (use /expt/app rather than /expt/data).
# <stage><inputfile> - Specify a single input file (full path).  The number
#             of batch jobs must be one.
# <stage><inputlist> - Specify input file list (a file containing a list
#             of input files, one per line, full path).
# <stage><inputdef>  - Specify input sam dataset definition.
#
#             It is optional to specify an input file or input list (Monte
#             Carlo generaiton doesn't need it, obviously).  It is also
#             optional for later production stages.  If no input is specified,
#             the list of files produced by the previous production stage
#             (if any) will be used as input to the current production stage
#             (must have been checked using option --check).
# <stage><numjobs> - Number of worker jobs (default 1).
# <stage><targetsize> - Specify target size for output files.  If specified,
#                       this attribute may override <numjobs> in the downward
#                       direction (i.e. <numjobs> is the maximum number of jobs).
# <stage><defname> - Sam output dataset defition name (default none).
# <stage><datatier> - Sam data tier (default none).
# <stage><initscript> - Worker initialization script (condor_lar.sh --init-script).
# <stage><initsource> - Worker initialization bash source script (condor_lar.sh --init-source).
# <stage><endscript>  - Worker end-of-job script (condor_lar.sh --end-script).
#                       Initialization/end-of-job scripts can be specified using an
#                       absolute or relative path relative to the current directory.
# <stage><histmerge>  - Name of histogram merging program or script (default specified
#                       at project level, or global default "hadd -T").
#
#
#
# <fcldir>  - Directory in which to search for fcl files (optional, repeatable).
#             Fcl files are searched for in the following directories, in order.
#             1.  Fcl directories specified using <fcldir>.
#             2.  The specified ubfcl version directory.
#             3.  $FHICL_FILE_PATH.
#             Regardless of where an fcl file is found, a copy is placed
#             in the work directory before job submission.  It is an
#             error of the fcl file isn't found.
#
######################################################################

import sys, os, stat, string, subprocess, shutil, urllib, getpass, json
from xml.dom.minidom import parse
import project_utilities, root_metadata

# Initialize the global symbol ROOT to be None.
# We will update this to point to the ROOT module later if we need to.

ROOT = None

# Do the same for samweb_cli module and global SAMWebClient object.

samweb_cli = None
samweb = None           # SAMWebClient object
extractor_dict = None   # Metadata extractor
proxy_ok = False

# XML exception class.

class XMLError(Exception):

    def __init__(self, s):
        self.value = s
        return

    def __str__(self):
        return self.value

# File-like class for writing files in pnfs space.
# Temporary file is opened in current directory and copied to final destination
# on close using ifdh.

class SafeFile:

    # Constructor.

    def __init__(self, destination=''):
        self.destination = ''
        self.filename = ''
        self.file = None
        if destination != '':
            self.open(destination)
        return

    # Open method.
    
    def open(self, destination):
        global proxy_ok
        if not proxy_ok:
            proxy_ok = project_utilities.test_proxy()
        self.destination = destination
        if project_utilities.safeexist(self.destination):
            subprocess.call(['ifdh', 'rm', self.destination])
        self.filename = os.path.basename(destination)
        if os.path.exists(self.filename):
            os.remove(self.filename)
        self.file = open(self.filename, 'w')
        return self.file

    # Write method.

    def write(self, line):
        self.file.write(line)
        return

    # Close method.

    def close(self):
        if self.file is not None and not self.file.closed:
            self.file.close()
        subprocess.call(['ifdh', 'cp', self.filename, self.destination])
        os.remove(self.filename)
        self.destination = ''
        self.filename = ''
        self.file = None
        return

# Function for opening files for writing using either python's built-in
# file object or SafeFile for dCache/pnfs files, as appropriate.

def safeopen(destination):
    if destination[0:6] == '/pnfs/':
        file = SafeFile(destination)
    else:
        file = open(destination, 'w')
    return file

# Stage definition class.

class StageDef:

    # Constructor.

    def __init__(self, stage_element, default_input_list, default_num_jobs, default_histmerge):

        # Assign default values.
        
        self.name = ''         # Stage name.
        self.fclname = ''      # Fcl name (just name, not path).
        self.outdir = ''       # Output directory.
        self.workdir = ''      # Work directory.
        self.inputfile = ''    # Single input file.
        self.inputlist = ''    # Input file list.
        self.inputdef = ''     # Input sam dataset definition.
        self.num_jobs = default_num_jobs # Number of jobs.
        self.target_size = 0   # Target size for output files.
        self.defname = ''      # Sam dataset definition name.
        self.data_tier = ''    # Sam data tier.
        self.init_script = ''  # Worker initialization script.
        self.init_source = ''  # Worker initialization bash source script.
        self.end_script = ''   # Worker end-of-job script.
        self.histmerge = default_histmerge # Histogram merging program

        # Extract values from xml.

        # Stage name (attribute).

        if stage_element.attributes.has_key('name'):
            self.name = stage_element.attributes['name'].firstChild.data
        if self.name == '':
            raise XMLError, "Stage name not specified."

        # Fcl file name (subelement).

        fclname_elements = stage_element.getElementsByTagName('fcl')
        if fclname_elements:
            self.fclname = fclname_elements[0].firstChild.data
        if self.fclname == '':
            raise XMLError, 'Fcl name not specified for stage %s.' % self.name

        # Output directory (subelement).

        outdir_elements = stage_element.getElementsByTagName('outdir')
        if outdir_elements:
            self.outdir = outdir_elements[0].firstChild.data
        if self.outdir == '':
            raise XMLError, 'Output directory not specified for stage %s.' % self.name

        # Work directory (subelement).

        workdir_elements = stage_element.getElementsByTagName('workdir')
        if workdir_elements:
            self.workdir = workdir_elements[0].firstChild.data
        if self.workdir == '':
            raise XMLError, 'Work directory not specified for stage %s.' % self.name

        # Single input file (subelement).

        inputfile_elements = stage_element.getElementsByTagName('inputfile')
        if inputfile_elements:
            self.inputfile = inputfile_elements[0].firstChild.data

        # Input file list (subelement).

        inputlist_elements = stage_element.getElementsByTagName('inputlist')
        if inputlist_elements:
            self.inputlist = inputlist_elements[0].firstChild.data

        # Input sam dataset dfeinition (subelement).

        inputdef_elements = stage_element.getElementsByTagName('inputdef')
        if inputdef_elements:
            self.inputdef = inputdef_elements[0].firstChild.data

        # It is an error to specify both input file and input list.

        if self.inputfile != '' and self.inputlist != '':
            raise XMLError, 'Input file and input list both specified for stage %s.' % self.name

        # It is an error to specify either input file or input list together
        # with a sam input dataset.

        if self.inputdef != '' and (self.inputfile != '' or self.inputlist != ''):
            raise XMLError, 'Input dataset and input files specified for stage %s.' % self.name

        # If none of input definition, input file, nor input list were specified, set
        # the input list to the dafault input list.

        if self.inputfile == '' and self.inputlist == '':
            self.inputlist = default_input_list

        # Number of jobs (subelement).

        num_jobs_elements = stage_element.getElementsByTagName('numjobs')
        if num_jobs_elements:
            self.num_jobs = int(num_jobs_elements[0].firstChild.data)

        # Target size for output files (subelement).

        target_size_elements = stage_element.getElementsByTagName('targetsize')
        if target_size_elements:
            self.target_size = int(target_size_elements[0].firstChild.data)

        # Sam dataset definition name (subelement).

        defname_elements = stage_element.getElementsByTagName('defname')
        if defname_elements:
            self.defname = defname_elements[0].firstChild.data

        # Sam data tier (subelement).

        data_tier_elements = stage_element.getElementsByTagName('datatier')
        if data_tier_elements:
            self.data_tier = data_tier_elements[0].firstChild.data

        # Worker initialization script (subelement).

        init_script_elements = stage_element.getElementsByTagName('initscript')
        if init_script_elements:
            self.init_script = init_script_elements[0].firstChild.data

        # Worker initialization source script (subelement).

        init_source_elements = stage_element.getElementsByTagName('initsource')
        if init_source_elements:
            self.init_source = init_source_elements[0].firstChild.data

        # Worker end-of-job script (subelement).

        end_script_elements = stage_element.getElementsByTagName('endscript')
        if end_script_elements:
            self.end_script = end_script_elements[0].firstChild.data

        # Histogram merging program.

        histmerge_elements = stage_element.getElementsByTagName('histmerge')
        if histmerge_elements:
            self.histmerge = histmerge_elements[0].firstChild.data

        # Done.

        return
        
    # String conversion.

    def __str__(self):
        result = 'Stage name = %s\n' % self.name
        result += 'Fcl filename = %s\n' % self.fclname
        result += 'Output directory = %s\n' % self.outdir
        result += 'Work directory = %s\n' % self.workdir
        result += 'Input file = %s\n' % self.inputfile
        result += 'Input list = %s\n' % self.inputlist
        result += 'Input sam dataset = %s\n' % self.inputdef
        result += 'Number of jobs = %d\n' % self.num_jobs
        result += 'Output file target size = %d\n' % self.target_size
        result += 'Dataset definition name = %s\n' % self.defname
        result += 'Data tier = %s' % self.defname
        result += 'Worker initialization script = %s\n' % self.init_script
        result += 'Worker initialization source script = %s\n' % self.init_source
        result += 'Worker end-of-job script = %s\n' % self.end_script
        result += 'Histogram merging program = %s\n' % self.histmerge
        return result

    # Raise an exception if any specified input file/list doesn't exist.
    # (We don't currently check sam input datasets).

    def checkinput(self):
        if self.inputfile != '' and not project_utilities.safeexist(self.inputfile):
            raise IOError, 'Input file %s does not exist.' % self.inputfile
        if self.inputlist != '' and not project_utilities.safeexist(self.inputlist):
            raise IOError, 'Input list %s does not exist.' % self.inputlist

        # If target size is nonzero, and input is from a file list, calculate
        # the ideal number of output jobs and override the current number 
        # of jobs.

        if self.target_size != 0 and self.inputlist != '':
            input_filenames = project_utilities.saferead(self.inputlist)
            size_tot = 0
            for line in input_filenames:
                filename = string.split(line)[0]
                filesize = os.stat(filename).st_size
                size_tot = size_tot + filesize
            new_num_jobs = size_tot / self.target_size
            if new_num_jobs < 1:
                new_num_jobs = 1
            if new_num_jobs > self.num_jobs:
                new_num_jobs = self.num_jobs
            print "Ideal number of jobs based on target file size is %d." % new_num_jobs
            if new_num_jobs != self.num_jobs:
                print "Updating number of jobs from %d to %d." % (self.num_jobs, new_num_jobs)
                self.num_jobs = new_num_jobs


    # Raise an exception if output directory or work directory doesn't exist.

    def checkdirs(self):
        if not os.path.exists(self.outdir):
            raise IOError, 'Output directory %s does not exist.' % self.outdir
        if not os.path.exists(self.workdir):
            raise IOError, 'Work directory %s does not exist.' % self.workdir
        return
    
    # Make output and work directory, if they don't exist.

    def makedirs(self):
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)
        self.checkdirs()
        return
    
# Project definition class contains data parsed from project defition xml file.

class ProjectDef:

    # Constructor.
    # project_element argument can be an xml element or None.

    def __init__(self, project_element, override_histmerge):

        # Assign default values.
        
        self.name= ''                     # Project name.
        self.group = ''                   # Experiment or group.
        if os.environ.has_key('GROUP'):
            self.group = os.environ['GROUP']
        self.num_events = 0               # Total events (all jobs).
        self.num_jobs = 1                 # Number of jobs.
        self.os = ''                      # Batch OS.
        self.histmerge = 'hadd -T'        # Default histogram merging program.
        self.release_tag = ''             # Larsoft release tag.
        self.release_qual = 'debug'       # Larsoft release qualifier.
        self.local_release_dir = ''       # Larsoft local release directory.
        self.local_release_tar = ''       # Larsoft local release tarball.
        self.ubfcl = ''                   # Ubfcl version.
        self.ubxml = ''                   # Ubxml version.
        self.file_type = ''               # Sam file type.
        self.run_type = ''                # Sam run type.
        self.script = 'condor_lar.sh'     # Batch script.
        self.start_script = 'condor_start_project.sh'  # Sam start project script.
        self.stop_script = 'condor_stop_project.sh'    # Sam stop project script.
        self.fclpath = []                 # Fcl search path.
        self.stages = []                  # List of stages (StageDef objects).
            
        # Extract values from xml.

        # Project name (attribute)

        if project_element.attributes.has_key('name'):
            self.name = project_element.attributes['name'].firstChild.data
        if self.name == '':
            raise XMLError, 'Project name not specified.'

        # Group (subelement).

        group_elements = project_element.getElementsByTagName('group')
        if group_elements:
            self.group = group_elements[0].firstChild.data
        if self.group == '':
            raise XMLError, 'Group not specified.'

        # Total events (subelement).

        num_events_elements = project_element.getElementsByTagName('numevents')
        if num_events_elements:
            self.num_events = int(num_events_elements[0].firstChild.data)
        if self.num_events == 0:
            raise XMLError, 'Number of events not specified.'

        # Number of jobs (subelement).

        num_jobs_elements = project_element.getElementsByTagName('numjobs')
        if num_jobs_elements:
            self.num_jobs = int(num_jobs_elements[0].firstChild.data)

        # OS (subelement).

        os_elements = project_element.getElementsByTagName('os')
        if os_elements:
            self.os = os_elements[0].firstChild.data

        # Histmerge (subelement).

        histmerge_elements = project_element.getElementsByTagName('histmerge')
        if histmerge_elements:
            if histmerge_elements[0].firstChild:
                self.histmerge = histmerge_elements[0].firstChild.data
            else:
                self.histmerge = ''
        if override_histmerge != '':
            self.histmerge = override_histmerge

        # Larsoft (subelement).

        larsoft_elements = project_element.getElementsByTagName('larsoft')
        if larsoft_elements:

            # Release tag (subelement).

            tag_elements = larsoft_elements[0].getElementsByTagName('tag')
            if tag_elements and tag_elements[0].firstChild != None:
                self.release_tag = tag_elements[0].firstChild.data

            # Release qualifier (subelement).

            qual_elements = larsoft_elements[0].getElementsByTagName('qual')
            if qual_elements:
                self.release_qual = qual_elements[0].firstChild.data

            # Local release directory or tarball (subelement).
            # 

            local_elements = larsoft_elements[0].getElementsByTagName('local')
            if local_elements:
                local = local_elements[0].firstChild.data
                if os.path.isdir(local):
                    self.local_release_dir = local
                else:
                    self.local_release_tar = local

        # Make sure local test release directory/tarball exists, if specified.
        # Existence of non-null local_release_dir has already been tested.

        if self.local_release_tar != '' and not os.path.exists(self.local_release_tar):
            raise IOError, "Local release directory/tarball %s does not exist." % self.local_release_tar
            
        # Ubfcl version (subelement).

        ubfcl_elements = project_element.getElementsByTagName('ubfcl')
        if ubfcl_elements:
            self.ubfcl = ubfcl_elements[0].firstChild.data

        # Ubxml version (subelement).

        ubxml_elements = project_element.getElementsByTagName('ubxml')
        if ubxml_elements:
            self.ubxml = ubxml_elements[0].firstChild.data

        # Sam file type (subelement).

        file_type_elements = project_element.getElementsByTagName('filetype')
        if file_type_elements:
            self.file_type = file_type_elements[0].firstChild.data

        # Sam run type (subelement).

        run_type_elements = project_element.getElementsByTagName('runtype')
        if run_type_elements:
            self.run_type = run_type_elements[0].firstChild.data

        # Batch script (subelement).

        script_elements = project_element.getElementsByTagName('script')
        if script_elements:
            self.script = script_elements[0].firstChild.data

        # Make sure batch script exists, and convert into a full path.

        script_path = ''
        try:
            proc = subprocess.Popen(['which', self.script], stdout=subprocess.PIPE)
            script_path = proc.stdout.readlines()[0].strip()
            proc.wait()
        except:
            pass
        if script_path == '' or not os.access(script_path, os.X_OK):
            raise IOError, 'Script %s not found.' % self.script
        self.script = script_path

        # Also convert start and stop project scripts into full path.
        # It is not treated as an error here if these aren't found.

        # Start project script.
        
        script_path = ''
        try:
            proc = subprocess.Popen(['which', self.start_script], stdout=subprocess.PIPE)
            script_path = proc.stdout.readlines()[0].strip()
            proc.wait()
        except:
            pass
        self.start_script = script_path

        # Stop project script.
        
        script_path = ''
        try:
            proc = subprocess.Popen(['which', self.stop_script], stdout=subprocess.PIPE)
            script_path = proc.stdout.readlines()[0].strip()
            proc.wait()
        except:
            pass
        self.stop_script = script_path

        # Fcl search path (repeatable subelement).

        fclpath_elements = project_element.getElementsByTagName('fcldir')
        for fclpath_element in fclpath_elements:
            self.fclpath.append(fclpath_element.firstChild.data)

        # Add ubfcl version directory.

        if self.ubfcl != '':

            # Loop over products directories.

            proc = subprocess.Popen(['ups', 'list', '-K@PROD_DIR', 'ubfcl', self.ubfcl],
                                    stdout=subprocess.PIPE)
            lines = proc.stdout.readlines()
            proc.wait()
            for line in lines:
                dir = line.replace('"','').strip()
                if os.path.exists(dir):
                    self.fclpath.append(dir)

        # Add $FHICL_FILE_PATH.

        if os.environ.has_key('FHICL_FILE_PATH'):
            for fcldir in string.split(os.environ['FHICL_FILE_PATH'], ':'):
                if os.path.exists(fcldir):
                    self.fclpath.append(fcldir)

        # Make sure all directories of fcl search path exist.

        for fcldir in self.fclpath:
            if not os.path.exists(fcldir):
                raise IOError, "Fcl search directory %s does not exist." % fcldir

        # Project stages (repeatable subelement).

        stage_elements = project_element.getElementsByTagName('stage')
        default_input_list = ''
        for stage_element in stage_elements:
            self.stages.append(StageDef(stage_element, 
                                        default_input_list, 
                                        self.num_jobs, 
                                        self.histmerge))
            default_input_list = os.path.join(self.stages[-1].outdir, 'files.list')

        # Done.
                
        return

    # String conversion.

    def __str__(self):
        result = 'Project name = %s\n' % self.name
        result += 'Group = %s\n' % self.group
        result += 'Total events = %d\n' % self.num_events
        result += 'Number of jobs = %d\n' % self.num_jobs
        result += 'OS = %s\n' % self.os
        result += 'Histogram merging program = %s\n' % self.histmerge
        result += 'Larsoft release tag = %s\n' % self.release_tag
        result += 'Larsoft release qualifier = %s\n' % self.release_qual
        result += 'Local test release directory = %s\n' % self.local_release_dir
        result += 'Local test release tarball = %s\n' % self.local_release_tar
        result += 'Ubfcl verion = %s\n' % self.ubfcl
        result += 'Ubxml verion = %s\n' % self.ubxml
        result += 'File type = %s\n' % self.file_type
        result += 'Run type = %s\n' % self.run_type
        result += 'Batch script = %s\n' % self.script
        result += 'Start sam project script = %s\n' % self.start_script
        result += 'Stop sam project script = %s\n' % self.stop_script
        result += 'Fcl search path:\n'
        for fcldir in self.fclpath:
            result += '    %s\n' % fcldir
        result += '\nStages:'
        for stage in self.stages:
            result += '\n\n' + str(stage)
        return result

    # Get the specified stage.

    def get_stage(self, stagename):

        if len(self.stages) == 0:
            raise LookupError, "Project does not have any stages."

        elif stagename == '' and len(self.stages) == 1:
            return self.stages[0]

        else:
            for stage in self.stages:
                if stagename == stage.name:
                    return stage

        # If we fell through to here, we didn't find an appropriate stage.

        raise RuntimeError, 'No stage %s.' % stagename

    # Find fcl file on fcl search path.

    def get_fcl(self, fclname):
        fcl = ''
        for fcldir in self.fclpath:
            fcl = os.path.join(fcldir, fclname)
            if os.path.exists(fcl):
                break
        if fcl == '' or not os.path.exists(fcl):
            raise IOError, 'Could not find fcl file %s.' % fclname
        return fcl


# Function to make sure ROOT module is imported.
# This function should be called before accessing global ROOT object.

def import_root():

    # Import ROOT module.
    # Globally turn off root warnings.
    # Don't let root see our command line options.

    global ROOT
    if ROOT is None:
        myargv = sys.argv
        sys.argv = myargv[0:1]
        import ROOT
        ROOT.gErrorIgnoreLevel = ROOT.kError
        sys.argv = myargv


# Function to make sure samweb_cli module is imported.
# Also initializes global SAMWebClient object.
# This function should be called before using samweb.

def import_samweb():

    # Import samweb_cli module, if not already done.

    global samweb_cli
    global samweb
    global extractor_dict

    exp = project_utilities.get_experiment()
    
    if samweb_cli == None:
        import samweb_cli
        samweb = samweb_cli.SAMWebClient(experiment=exp)
        import extractor_dict


# Clean function.

def doclean(project, stagename):

    # Loop over stages.
    # Clean all stages beginning with the specified stage.
    # For empty stagename, clean all stages.
    #
    # For safety, only clean directories if the uid of the
    # directory owner matches the current uid or effective uid.
    # Do this even if the delete operation is allowed by filesystem
    # permissions (directories may be group- or public-write
    # because of batch system).

    match = 0
    uid = os.getuid()
    euid = os.geteuid()
    for stage in project.stages:
        if stagename == '' or stage.name == stagename:
            match = 1
        if match:

            # Clean this stage outdir.

            if os.path.exists(stage.outdir):
                dir_uid = os.stat(stage.outdir).st_uid
                if dir_uid == uid or dir_uid == euid:
                    print 'Clean directory %s.' % stage.outdir
                    shutil.rmtree(stage.outdir)
                else:
                    print 'Owner mismatch, delete %s manually.' % stage.outdir
                    sys.exit(1)

            # Clean this stage workdir.

            if os.path.exists(stage.workdir):
                dir_uid = os.stat(stage.workdir).st_uid
                if dir_uid == uid or dir_uid == euid:
                    print 'Clean directory %s.' % stage.workdir
                    shutil.rmtree(stage.workdir)
                else:
                    print 'Owner mismatch, delete %s manually.' % stage.workdir
                    sys.exit(1)

    # Done.

    return

# Stage status fuction.

def dostatus(project):

    print 'Project %s:' % project.name

    # Loop over stages.

    for stage in project.stages:

        stagename = stage.name
        outdir = stage.outdir

        # Test whether output directory exists.

        if project_utilities.safeexist(outdir):

            # Output directory exists.

            nfile = 0
            nev = 0
            nmiss = 0

            # Count good files and events.

            eventsfile = os.path.join(outdir, 'events.list')
            if project_utilities.safeexist(eventsfile):
                lines = project_utilities.saferead(eventsfile)
                for line in lines:
                    words = string.split(line)
                    if len(words) >= 2:
                        nfile = nfile + 1
                        nev = nev + int(words[1])

            # Count bad files.

            missingfile = os.path.join(outdir, 'missing.txt')
            if project_utilities.safeexist(missingfile):
                lines = project_utilities.saferead(missingfile)
                nmiss = nmiss + len(lines)

            print 'Stage %s: %d good output files, %d events, %d bad or missing output files.' % (
                stagename, nfile, nev, nmiss)

        else:

            # Output directory does not exist.

            print 'Stage %s output directory does not exist.' % stagename

    return

# This is a recursive function that looks for project
# subelements within the input element.  It returns the
# first matching element that it finds, or None.

def find_project(element, projectname):

    # First check if the input element is a project.
    # If it is, and if the name matches, return that.

    if element.nodeName == 'project' and \
       (projectname == '' or \
        (element.attributes.has_key('name') and
         projectname == element.attributes['name'].firstChild.data)):
        return element

    # Loop over subelements.

    subelements = element.getElementsByTagName('*')
    for subelement in subelements:
        project = find_project(subelement, projectname)
        if project is not None:
            return project

    # If we fell out of the loop, return None.

    return None
    

# Extract the specified project element from xml file.
# Project elements can exist at any depth inside xml file.
# Return project Element or None.

def get_project(xmlfile, projectname):

    # Parse xml (returns xml document).

    if xmlfile == '-':
        xml = sys.stdin
    else:
        xml = urllib.urlopen(xmlfile)
    doc = parse(xml)

    # Extract root element.

    root = doc.documentElement

    # Find project element.
    
    project = find_project(root, projectname)
    if project is None:
        return None

    # Check that project element has "name" attribute.
    
    if not project.attributes.has_key('name'):
        print 'Project element does not have "name" attribute.'
        return None

    # Success.

    return project

# Parse directory name of type <cluster>_<process> and return
# a 2-tuple of integers or None.

def parsedir(dirname):

    # This method separates a directory name of the form
    # <cluster>_<process> into a 2-tuple of integers (cluster, process).
    # If the directory is not of the proper form, return None.

    result = None

    if string.count(dirname, '_') == 1:
        n = string.index(dirname, '_')
        try:
            cluster = int(dirname[0:n])
            process = int(dirname[n+1:])
            result = (cluster, process)
        except:
            pass
    
    return result

# Check a single root file by reading precalculted metadata information.

def check_root_json(json_path):

    # Default result

    nevroot = -1

    # Read contents of json file and convert to a single string.

    lines = project_utilities.saferead(json_path)
    s = ''
    for line in lines:
        s = s + line

    # Convert json string to python dictionary.

    md = json.loads(s)

    # Extract number of events from metadata.

    if len(md.keys()) > 0:
        nevroot = 0
        if md.has_key('events'):
            nevroot = int(md['events'])

    return nevroot


# Check a single root file.
# Returns an int containing following information.
# 1.  Number of event (>0) in TTree named "Events."
# 2.  Zero if root file does not contain an Events TTree, but is otherwise valid (openable).
# 3.  -1 for error (root file is not openable).

def check_root_file(path):

    global proxy_ok
    nevroot = -1

    # See if we have precalculated metadata for this root file.

    json_path = path + '.json'
    if project_utilities.safeexist(json_path):

        # Get number of events from precalculated metadata.

        nevroot = check_root_json(json_path)

    else:

        # Make root metadata.

        url = project_utilities.path_to_url(path)
        if url != path and not proxy_ok:
            proxy_ok = project_utilities.test_proxy()
        print 'Generating root metadata for file %s.' % os.path.basename(path)
        md = root_metadata.get_external_metadata(path)
        if md.has_key('events'):

            # Art root file if dictionary has events key.

            nevroot = int(md['events'])

        elif len(md.keys()) > 0:

            # No events key, but non-empty dictionary, so histo/ntuple root file.

            nevroot = 0
        else:

            # Empty dictionary is invalid root file.

            nevroot = -1

        # Save root metadata in .json file.

        mdtext = json.dumps(md, sys.stdout, indent=2, sort_keys=True)
        json_file = safeopen(json_path)
        json_file.write(mdtext)
        json_file.close()

    return nevroot


# Check root files (*.root) in the specified directory.

def check_root(dir):

    # This method looks for files with names of the form *.root.
    # If such files are found, it also checks for the existence of
    # an Events TTree.
    #
    # Returns a 3-tuple containing the following information.
    # 1.  Total number of events in art root files.
    # 2.  A list of 2-tuples with an entry for each art root file.
    #     The 2-tuple contains the following information.
    #     a) Number of events
    #     b) Filename.
    # 3.  A list of histogram root files.

    nev = -1
    roots = []
    hists = []

    print 'Checking root files in directory %s.' % dir
    filenames = os.listdir(dir)
    for filename in filenames:
        if filename[-5:] == '.root':
            path = os.path.join(dir, filename)
            nevroot = check_root_file(path)
            if nevroot > 0:
                if nev < 0:
                    nev = 0
                nev = nev + nevroot
                roots.append((os.path.join(dir, filename), nevroot))

            elif nevroot == 0:

                # Valid root (histo/ntuple) file, not an art root file.

                hists.append(os.path.join(dir, filename))

            else:

                # Found a .root file that is not openable.
                # Print a warning, but don't trigger any other error.

                print 'Warning: File %s in directory %s is not a valid root file.' % (filename, dir)

    # Done.

    return (nev, roots, hists)
    

# Get the list of input files for a project stage.

def get_input_files(stage):

    # In case of single file or file list input, files are returned exactly 
    # as specified, which would normallly be as the full path.
    # In case of sam input, only the file names are returned (guaranteed unique).

    result = []
    if stage.inputfile != '':
        result.append(stage.inputfile)

    elif stage.inputlist != '' and project_utilities.safeexist(stage.inputlist):
        try:
            input_filenames = project_utilities.saferead(stage.inputlist)
            for line in input_filenames:
                words = string.split(line)
                result.append(words[0])
        except:
            pass

    elif stage.inputdef != '':
        import_samweb()
        result = samweb.listFiles(defname=stage.inputdef)

    # Done.

    return result

# Check project results in the specified directory.

def docheck(dir, num_events, num_jobs, has_input_files, input_def, ana, has_metadata, histmerge):

    # This method expects to find several subdirectories named as
    # <cluster>_<process>, where the <process> part of the subdirectory
    # name is a small integer.
    #
    # In case of file (list) input, or no input, we expect to have unique
    # process numbers in the range [0, num_jobs).  In the case of sam input,
    # process numbers are not necessarily unique, and they have no fixed
    # range.  However the pair (cluster, process) should still be unique.
    #
    # Return 0 if all checks are OK, otherwise return number of missing files.
    #
    # The following checks are performed.
    #
    # 1.  Make sure subdirectories are as expected (complain about missing
    #     subdirectories).
    #
    # 2.  Look for at least one art root file in each subdirectory (*.root)
    #     containing a valid Events TTree.  Complain about any
    #     that do not contain such a root file.
    #
    # 3.  Check that the number of events in the Events tree are as expected.
    #
    # 4.  Complain about any duplicated art root file names (if sam metadata is defined).
    #
    # 5.  Check job exit status (saved in lar.stat).
    #
    # 6.  For sam input, make sure that file sam_project.txt and cpid.txt is present.
    #
    # 7.  Check that any non-art root files are openable.
    #
    # In analysis mode (if argumment ana != 0), skip checks 2-4, but still do
    # checks 1 and 5-7.
    #
    # This function also creates the following files in the specified directory.
    #
    # 1.  files.list  - List of good root files.
    # 2.  events.list - List of good root files and number of events in each file.
    # 3.  missing.txt - List of missing (cluster, process).
    # 4.  sam_projects.list - List of successful sam projects.
    # 5.  cpids.list        - list of successful consumer process ids.
    # 6.  hists.list  - List of non-art histogram root files (histograms and/or ntuples).
    # 7.  hist.root   - Merged histogram file (excludes TTrees).
    #                   Made using "hadd -T hist.root @hlists.list" (or equivalent).

    import_samweb()

    print 'Checking directory %s' % dir

    # Default result is success.

    result = 0

    # Count total number of events and root files.

    nev_tot = 0
    nroot_tot = 0

    # Loop over subdirectories (ignore files and directories named *_start and *_stop).

    cluster = 0    # Most recent cluster.
    process = 0    # Most recent process.
    procmap = {}   # procmap[process][cluster] = <list of art root files and event counts>
    hists = []     # List of non-art root files.
    sam_projects = []
    cpids = []

    subdirs = os.listdir(dir)
    for subdir in subdirs:
        subpath = os.path.join(dir, subdir)
        dirok = project_utilities.fast_isdir(subpath)

        # Update list of sam projects from start job.

        if dirok and subpath[-6:] == '_start':
            filename = os.path.join(subpath, 'sam_project.txt')
            if project_utilities.safeexist(filename):
                sam_project = project_utilities.saferead(filename)[0].strip()
                if sam_project != '' and not sam_project in sam_projects:
                    sam_projects.append(sam_project)

        # Regular worker jobs checked here.

        if dirok and not subpath[-6:] == '_start' and not subpath[-5:] == '_stop':

            # Found a subdirectory.

            # First check for naming problems.
            
            bad = 0
            clus_proc = parsedir(subdir)
            if clus_proc is None:
                bad = 1
            else:
                cluster = clus_proc[0]
                process = clus_proc[1]

                # Check process.

                if input_def == '':
                    if process < 0 or process >= num_jobs:
                        print 'Process number %d out of range.' % process
                        bad = 1
                    if procmap.has_key(process):
                        print 'Duplicate process number %d.' % process
                        bad = 1

            # Check lar exit status (if any).

            if not bad:
                stat_filename = os.path.join(subpath, 'lar.stat')
                if project_utilities.safeexist(stat_filename):
                    status = 0
                    try:
                        status = int(project_utilities.saferead(stat_filename)[0].strip())
                        if status != 0:
                            print 'Job in subdirectory %s ended with non-zero exit status %d.' % (
                                subdir, status)
                            bad = 1
                    except:
                        print 'Bad file lar.stat in subdirectory %s.' % subdir
                        bad = 1

            # Check existence of sam_project.txt.
            # Update sam_projects.

            if not bad and input_def != '':
                filename = os.path.join(subpath, 'sam_project.txt')
                if not project_utilities.safeexist(filename):
                    bad = 1
                if not bad:
                    sam_project = project_utilities.saferead(filename)[0].strip()
                    if not sam_project in sam_projects:
                        sam_projects.append(sam_project)

            # Check existence of cpid.txt
            # Update cpids.

            if not bad and input_def != '':
                filename = os.path.join(subpath, 'cpid.txt')
                if not project_utilities.safeexist(filename):
                    bad = 1
                if not bad:
                    cpid = project_utilities.saferead(filename)[0].strip()
                    if not cpid in cpids:
                        cpids.append(cpid)

            # Now check root files in this subdirectory.

            if not bad:
                nev = 0
                roots = []
                nev, roots, subhists = check_root(subpath)
                if not ana:
                    if len(roots) == 0 or nev < 0:
                        print 'Problem with root file(s) in subdirectory %s.' % subdir
                        bad = 1

                if not bad:

                    # Check if any root file names match previously found root file names.
                    # If so, count both this process as bad.

                    if has_metadata:
                        for root in roots:
                            rootname = os.path.basename(root[0])
                            for p in procmap.keys():
                                for c in procmap[p].keys():
                                    oldroots = procmap[p][c]
                                    for oldroot in oldroots:
                                        oldrootname = os.path.basename(oldroot[0])
                                        if rootname == oldrootname:
                                            print 'Duplicate filename %s in subdirectory %s' % (rootname,
                                                                                                subdir)
                                            olddir = os.path.basename(os.path.dirname(oldroot[0]))
                                            print 'Previous subdirectory %s' % olddir
                                            bad = 1

                    if not bad:
                    
                        # Save good root files.
                    
                        if not procmap.has_key(process):
                            procmap[process] = {}
                        procmap[process][cluster] = roots

                        # Save good histogram files.

                        hists.extend(subhists)

                        # Count good events and root files.

                        nev_tot = nev_tot + nev
                        nroot_tot = nroot_tot + len(roots)

                        # Save good transferred files.

            # Print/save result of checks for one subdirectory.

            if bad:
                print 'Bad subdirectory %s.' % subdir

    # Done looping over subdirectoryes.
    # Dictionary procmap now contains a list of good processes
    # and root files.

    # Open files.

    filelistname = os.path.join(dir, 'files.list')
    filelist = safeopen(filelistname)

    eventslistname = os.path.join(dir, 'events.list')
    eventslist = safeopen(eventslistname)

    missingname = os.path.join(dir, 'missing.txt')
    missing = safeopen(missingname)

    histlistname = os.path.join(dir, 'hists.list')
    histlist = safeopen(histlistname)

    histurlsname_temp = 'histurls.list'
    if os.path.exists(histurlsname_temp):
        os.remove(histurlsname_temp)
    histurls = open(histurlsname_temp, 'w')

    # See if there are any missing processes and generate file "missing.txt."
    # Skip this step for sam input.

    nmiss = 0
    nerror = 0
    if input_def == '':
        for process in xrange(num_jobs):
            if not procmap.has_key(process):
                nmiss = nmiss + 1
                nerror = nerror + 1
                missing.write('%d %d\n' % (cluster, process))
                print 'Missing process %d.' % process

    # Generate "files.list" and "events.list."

    nproc = 0
    for p in procmap.keys():
        for c in procmap[p].keys():
            nproc = nproc + 1
            for root in procmap[p][c]:
                filelist.write('%s\n' % root[0])
                eventslist.write('%s %d\n' % root)

    # Generate "hists.list."

    for hist in hists:
        histlist.write('%s\n' % hist)
        histurls.write('%s\n' % project_utilities.path_to_url(hist))

    # Print summary.

    if ana:
        print "%d processes completed successfully." % nproc
        print "%d total good histogram files." % len(hists)
    else:
        print "%d total good events." % nev_tot
        print "%d total good root files." % nroot_tot
        print "%d total good histogram files." % len(hists)

    # Close files.

    filelist.close()
    eventslist.close()
    missing.close()
    histlist.close()
    histurls.close()

    # Make merged histogram file using histmerge.

    if len(hists) > 0 and histmerge != '':
        print "Merging %d histogram files using %s." % (len(hists), histmerge)

        histname = os.path.join(dir, 'hist.root')
        if histname[0:6] == '/pnfs/':
            histname_temp = 'hist.root'
        else:
            histname_temp = histname
        if os.path.exists(histname_temp):
            os.remove(histname_temp)
        comlist = histmerge.split()
        comlist.extend(["-v", "0", "-f", "-k", histname_temp, '@' + histurlsname_temp])
        rc = subprocess.call(comlist)
        if rc != 0:
            print "%s exit status %d" % (histmerge, rc)
        if histname != histname_temp:
            if project_utilities.safeexist(histname):
                os.remove(histname)
            if os.path.exists(histname_temp):
                subprocess.call(['ifdh', 'cp', histname_temp, histname])
                os.remove(histname_temp)
    os.remove(histurlsname_temp)

    # Make sam files.

    if input_def != '':

        # List of successful sam projects.
        
        sam_projects_filename = os.path.join(dir, 'sam_projects.list')
        sam_projects_file = safeopen(sam_projects_filename)
        for sam_project in sam_projects:
            sam_projects_file.write('%s\n' % sam_project)
        sam_projects_file.close()

        # List of successfull consumer process ids.

        cpids_filename = os.path.join(dir, 'cpids.list')
        cpids_file = safeopen(cpids_filename)
        for cpid in cpids:
            cpids_file.write('%s\n' % cpid)
        cpids_file.close()

        # Get number of consumed files.

        cpids_list = ''
        sep = ''
        for cpid in cpids:
            cpids_list = cpids_list + '%s%s' % (sep, cpid)
            sep = ','
        if cpids_list != '':
            dim = 'consumer_process_id %s and consumed_status consumed' % cpids_list
            import_samweb()
            nconsumed = samweb.countFiles(dimensions=dim)
        else:
            nconsumed = 0

        # Get number of unconsumed files.

        if cpids_list != '':
            udim = '(defname: %s) minus (%s)' % (input_def, dim)
        else:
            udim = 'defname: %s' % input_def
        nunconsumed = samweb.countFiles(dimensions=udim)
        nerror = nerror + nunconsumed
        
        # Sam summary.

        print '%d sam projects.' % len(sam_projects)
        print '%d successful consumer process ids.' % len(cpids)
        print '%d files consumed.' % nconsumed
        print '%d files not consumed.' % nunconsumed

        # Check project statuses.
        
        for sam_project in sam_projects:
            print '\nChecking sam project %s' % sam_project
            import_samweb()
            url = samweb.findProject(sam_project, project_utilities.get_experiment())
            if url != '':
                result = samweb.projectSummary(url)
                nd = 0
                nc = 0
                nf = 0
                nproc = 0
                nact = 0
                if result.has_key('consumers'):
                    consumers = result['consumers']
                    for consumer in consumers:
                        if consumer.has_key('processes'):
                            processes = consumer['processes']
                            for process in processes:
                                nproc = nproc + 1
                                if process.has_key('status'):
                                    if process['status'] == 'active':
                                        nact = nact + 1
                                if process.has_key('counts'):
                                    counts = process['counts']
                                    if counts.has_key('delivered'):
                                        nd = nd + counts['delivered']
                                    if counts.has_key('consumed'):
                                        nc = nc + counts['consumed']
                                    if counts.has_key('failed'):
                                        nf = nf + counts['failed']
                print 'Status: %s' % result['project_status']
                print '%d total processes' % nproc
                print '%d active processes' % nact
                print '%d files in snapshot' % result['files_in_snapshot']
                print '%d files delivered' % (nd + nc)
                print '%d files consumed' % nc
                print '%d files failed' % nf
                print

    # Done

    print '%d processes with errors.' % nerror
    return nmiss

# Construct dimension string for project, stage.

def dimensions(project, stage):

    dim = 'file_type %s' % project.file_type
    dim = dim + ' and data_tier %s' % stage.data_tier
    dim = dim + ' and ub_project.name %s' % project.name
    dim = dim + ' and ub_project.stage %s' % stage.name
    if project.ubxml != '': 
        dim = dim + ' and ub_project.version %s' % project.ubxml
    else:
        dim = dim + ' and ub_project.version %s' % project.release_tag
    dim = dim + ' and availability: anylocation'
    return dim

# Check sam declarations.

def docheck_declarations(outdir, declare):

    # Initialize samweb.

    import_samweb()

    # Loop over root files listed in files.list.

    roots = []
    fnlist = os.path.join(outdir, 'files.list')
    if project_utilities.safeexist(fnlist):
        roots = project_utilities.saferead(fnlist)
    else:
        print 'No files.list file found, run project.py --check'
        sys.exit(1)

    for root in roots:
        path = string.strip(root)
        fn = os.path.basename(path)

        # Check metadata

        has_metadata = False
        try:
            md = samweb.getMetadata(filenameorid=fn)
            has_metadata = True
        except samweb_cli.exceptions.FileNotFound:
            pass

        # Report or declare file.

        if has_metadata:
            print 'Metadata OK: %s' % fn
        else:
            if declare:
                print 'Declaring: %s' % fn
                md = extractor_dict.getmetadata(path)
                samweb.declareFile(md=md)
            else:
                print 'Not declared: %s' % fn

    return 0

# Print summary of files returned by sam query.

def dotest_declarations(dim):

    # Initialize samweb.

    import_samweb()

    # Do query

    result = samweb.listFilesSummary(dimensions=dim)
    for key in result.keys():
        print '%s: %s' % (key, result[key])

    return 0

# Check sam dataset definition.

def docheck_definition(defname, dim, define):

    if defname == '':
        return 1

    # Initialize samweb.

    import_samweb()

    # See if this definition already exists.

    def_exists = False
    try:
        desc = samweb.descDefinition(defname=defname)
        def_exists = True
    except samweb_cli.exceptions.DefinitionNotFound:
        pass

    # Make report and maybe make definition.

    if def_exists:
        print 'Definition already exists: %s' % defname
    else:
        if define:
            print 'Creating definition %s.' % defname
            samweb.createDefinition(defname=defname, dims=dim)
        else:
            print 'Definition should be created: %s' % defname

    return 0

# Print summary of files returned by dataset definition.

def dotest_definition(defname):

    # Initialize samweb.

    import_samweb()

    # Do query

    result = samweb.listFilesSummary(defname=defname)
    for key in result.keys():
        print '%s: %s' % (key, result[key])

    return 0

# Delete sam dataset definition.

def doundefine(defname):

    if defname == '':
        return 1

    # Initialize samweb.

    import_samweb()

    # See if this definition already exists.

    def_exists = False
    try:
        desc = samweb.descDefinition(defname=defname)
        def_exists = True
    except samweb_cli.exceptions.DefinitionNotFound:
        pass

    # Make report and maybe make definition.

    if def_exists:
        print 'Deleting definition: %s' % defname
        samweb.deleteDefinition(defname=defname)
    else:
        print 'No such definition: %s' % defname

    return 0

# Check disk locations.  Maybe add or remove locations.

def docheck_locations(dim, outdir, add, clean, remove, upload):

    # Initialize samweb.

    import_samweb()

    # Loop over files queried by dimension string.

    filelist = samweb.listFiles(dimensions=dim, stream=True)
    while 1:
        try:
            filename = filelist.next()
        except StopIteration:
            break

        # Got a filename.

        # Look for locations on disk.
        # Look in first level subdirectories of outdir.

        disk_locs = []
        for subdir in os.listdir(outdir):
            subpath = os.path.join(outdir, subdir)
            if project_utilities.fast_isdir(subpath):
                for fn in os.listdir(subpath):
                    if fn == filename:
                        filepath = os.path.join(subpath, fn)
                        disk_locs.append(os.path.dirname(filepath))

        # Also get sam locations.

        sam_locs = samweb.locateFile(filenameorid=filename)
        if len(sam_locs) == 0 and not upload:
            print 'No location: %s' % filename

        # Make a double loop over disk and sam locations, in order
        # to identify locations that should added.
        # Note that we ignore the node part of the sam location.

        locs_to_add = []
        for disk_loc in disk_locs:
            should_add = True
            for sam_loc in sam_locs:
                if sam_loc['location_type'] == 'disk':
                    if disk_loc == sam_loc['location'].split(':')[-1]:
                        should_add = False
                        break
            if should_add:
                locs_to_add.append(disk_loc)

        # Loop over sam locations, in order to identify locations
        # that should be removed.  Note that for this step, we don't
        # necessarily assume that we found the disk location
        # in the directory search above, rather check the existence
        # of the file directly.

        locs_to_remove = []
        for sam_loc in sam_locs:
            if sam_loc['location_type'] == 'disk':

                # If remove is specified, uncondiontally remove this location.

                if remove:
                    locs_to_remove.append(sam_loc['location'])

                # Otherwise, check if file exists.

                else:

                    # Split off the node, if any, from the location.

                    local_path = os.path.join(sam_loc['location'].split(':')[-1], filename)
                    if not project_utilities.safeexist(local_path):
                        locs_to_remove.append(sam_loc['location'])

        # Loop over sam locations and identify files that can be uploaded.
        # If this file has no disk locations, don't do anything (not an error).
        # In case we decide to upload this file, always upload from the first
        # disk location.

        locs_to_upload = {}    # locs_to_upload[disk-location] = dropbox-directory
        should_upload = False
        if upload and len(disk_locs) > 0:
            should_upload = True
            for sam_loc in sam_locs:
                if sam_loc['location_type'] == 'tape':
                    should_upload = False
                    break
            if should_upload:
                dropbox = project_utilities.get_dropbox(filename)
                locs_to_upload[disk_locs[0]] = dropbox
        
        # Report results and do the actual adding/removing/uploading.

        for loc in locs_to_add:
            node = project_utilities.get_bluearc_server()
            if loc[0:6] == '/pnfs/':
                node = project_utilities.get_dcache_server()
            loc = node + loc.split(':')[-1]
            if add:
                print 'Adding location: %s.' % loc
                samweb.addFileLocation(filenameorid=filename, location=loc)
            elif not upload:
                print 'Can add location: %s.' % loc

        for loc in locs_to_remove:
            if clean or remove:
                print 'Removing location: %s.' % loc
                samweb.removeFileLocation(filenameorid=filename, location=loc)
            elif not upload:
                print 'Should remove location: %s.' % loc

        for loc in locs_to_upload.keys():
            dropbox = locs_to_upload[loc]

            # Make sure dropbox directory exists.
            
            if not os.path.isdir(dropbox):
                print 'Dropbox directory %s does not exist.' % dropbox
            else:

                # Test whether this file has already been copied to dropbox directory.
                
                dropbox_filename = os.path.join(dropbox, filename)
                if project_utilities.safeexist(dropbox_filename):
                    print 'File %s already exists in dropbox %s.' % (filename, dropbox)
                else:

                    # Copy file to dropbox.

                    loc_filename = os.path.join(loc, filename)
                    print 'Copying %s to dropbox directory %s.' % (filename, dropbox)
                    subprocess.call(['ifdh', 'cp', loc_filename, dropbox_filename])

    return 0

# Check disk locations.  Maybe add or remove locations.

def docheck_tape(dim):

    # Initialize samweb.

    import_samweb()

    # Loop over files queried by dimension string.

    nbad = 0
    ntot = 0
    filelist = samweb.listFiles(dimensions=dim, stream=True)
    while 1:
        try:
            filename = filelist.next()
        except StopIteration:
            break

        # Got a filename.

        ntot = ntot + 1

        # Look for sam tape locations.

        is_on_tape = False
        sam_locs = samweb.locateFile(filenameorid=filename)
        for sam_loc in sam_locs:
            if sam_loc['location_type'] == 'tape':
                is_on_tape = True
                break

        if is_on_tape:
            print 'On tape: %s' % filename
        else:
            nbad = nbad + 1
            print 'Not on tape: %s' % filename

    print '%d files.' % ntot
    print '%d files need to be store on tape.' % nbad

    return 0

# Print help.

def help():

    filename = sys.argv[0]
    file = open(filename, 'r')

    doprint=0
    
    for line in file.readlines():
        if line[2:12] == 'project.py':
            doprint = 1
        elif line[0:6] == '######' and doprint:
            doprint = 0
        if doprint:
            if len(line) > 2:
                print line[2:],
            else:
                print


# Print xml help.

def xmlhelp():

    filename = sys.argv[0]
    file = open(filename, 'r')

    doprint=0
    
    for line in file.readlines():
        if line[2:20] == 'XML file structure':
            doprint = 1
        elif line[0:6] == '######' and doprint:
            doprint = 0
        if doprint:
            if len(line) > 2:
                print line[2:],
            else:
                print


# Main program.

def main(argv):

    # Parse arguments.

    xmlfile = ''
    projectname = ''
    stagename = ''
    override_histmerge = ''
    submit = 0
    check = 0
    checkana = 0
    stage_status = 0
    makeup = 0
    clean = 0
    outdir = 0
    fcl = 0
    defname = 0
    do_input_files = 0
    declare = 0
    define = 0
    undefine = 0
    check_declarations = 0
    test_declarations = 0
    check_definition = 0
    test_definition = 0
    add_locations = 0
    check_locations = 0
    upload = 0
    check_tape = 0
    clean_locations = 0
    remove_locations = 0

    args = argv[1:]
    while len(args) > 0:
        if args[0] == '-h' or args[0] == '--help' :
            help()
            return 0
        elif args[0] == '-xh' or args[0] == '--xmlhelp' :
            xmlhelp()
            return 0
        elif args[0] == '--xml' and len(args) > 1:
            xmlfile = args[1]
            del args[0:2]
        elif args[0] == '--project' and len(args) > 1:
            projectname = args[1]
            del args[0:2]
        elif args[0] == '--stage' and len(args) > 1:
            stagename = args[1]
            del args[0:2]
        elif args[0] == '--histmerge' and len(args) > 1:
            override_histmerge = args[1]
            del args[0:2]
        elif args[0] == '--tmpdir' and len(args) > 1:
            os.environ['TMPDIR'] = args[1]
            del args[0:2]
        elif args[0] == '--submit':
            submit = 1
            del args[0]
        elif args[0] == '--check':
            check = 1
            del args[0]
        elif args[0] == '--checkana':
            checkana = 1
            del args[0]
        elif args[0] == '--status':
            stage_status = 1
            del args[0]
        elif args[0] == '--makeup':
            makeup = 1
            del args[0]
        elif args[0] == '--clean':
            clean = 1
            del args[0]
        elif args[0] == '--outdir':
            outdir = 1
            del args[0]
        elif args[0] == '--fcl':
            fcl = 1
            del args[0]
        elif args[0] == '--defname':
            defname = 1
            del args[0]
        elif args[0] == '--input_files':
            do_input_files = 1
            del args[0]
        elif args[0] == '--declare':
            declare = 1
            del args[0]
        elif args[0] == '--define':
            define = 1
            del args[0]
        elif args[0] == '--undefine':
            undefine = 1
            del args[0]
        elif args[0] == '--check_declarations':
            check_declarations = 1
            del args[0]
        elif args[0] == '--test_declarations':
            test_declarations = 1
            del args[0]
        elif args[0] == '--check_definition':
            check_definition = 1
            del args[0]
        elif args[0] == '--test_definition':
            test_definition = 1
            del args[0]
        elif args[0] == '--add_locations':
            add_locations = 1
            del args[0]
        elif args[0] == '--check_locations':
            check_locations = 1
            del args[0]
        elif args[0] == '--upload':
            upload = 1
            del args[0]
        elif args[0] == '--check_tape':
            check_tape = 1
            del args[0]
        elif args[0] == '--clean_locations':
            clean_locations = 1
            del args[0]
        elif args[0] == '--remove_locations':
            remove_locations = 1
            del args[0]
        else:
            print 'Unknown option %s' % args[0]
            return 1

    # Make sure xmlfile was specified.

    if xmlfile == '':
        print 'No xml file specified.  Type "project.py -h" for help.'
        return 1
    
    # Make sure that no more than one action was specified (except clean and info options).

    num_action = submit + check + checkana + stage_status + makeup + define + undefine + declare
    if num_action > 1:
        print 'More than one action was specified.'
        return 1

    # Get the project element.

    project_element = get_project(xmlfile, projectname)
    if project_element is None:
        if projectname == '':
            print 'Could not find unique project in xml file %s.' % xmlfile
        else:
            print 'Could not find project %s in xml file %s.' % (projectname, xmlfile)
        return 1

    # Convert the project element into a ProjectDef object.

    project = ProjectDef(project_element, override_histmerge)

    # Make sure that we have a kerberos ticket if we might need one to submit jobs.

    if submit or makeup:
        ok = subprocess.call(['klist', '-s'])
        if ok != 0:
            print 'Please get a kerberos ticket.'
            return 1

    # Do clean action now.  Cleaning can be combined with submission.

    if clean:
        doclean(project, stagename)

    # Do stage_status now.

    if stage_status:
        dostatus(project)
        return 0

    # Get the current stage definition.

    stage = project.get_stage(stagename)

    # Do outdir action now.

    if outdir:
        print stage.outdir

    # Do fcl action now.

    if fcl:
        print stage.fclname, project.ubfcl

    # Do defname action now.

    if defname:
        if stage.defname != '':
            print stage.defname

    # Do input_names action now.

    if do_input_files:
        input_files = get_input_files(stage)
        for input_file in input_files:
            print input_file

    # Maybe make output and work directories.

    if submit:
        stage.makedirs()

    # Check input file/list, output directory, and work directory.

    if submit or check or checkana or makeup:
        stage.checkinput()
        stage.checkdirs()

        # If output is on dcache, make output directory group-writable.

        if stage.outdir[0:6] == '/pnfs/':
            mode = os.stat(stage.outdir).st_mode
            if not mode & stat.S_IWGRP:
                mode = mode | stat.S_IWGRP
                os.chmod(stage.outdir, mode)

        # For now, also make output directory world-writable.

        #if not mode & stat.S_IWOTH:
        #    mode = mode | stat.S_IWOTH
        #    os.chmod(stage.outdir, mode)

    # For first submission, make sure output directory is empty.

    if submit and len(os.listdir(stage.outdir)) != 0:
        raise RuntimeError, 'Output directory %s is not empty.' % stage.outdir

    # Do the following sections only if there is the possibility
    # of submtting jobs (submit or makeup action).

    if submit or makeup:

        # If there is an input list, copy it to the work directory.

        input_list_name = ''
        if stage.inputlist != '':
            input_list_name = os.path.basename(stage.inputlist)
            work_list_name = os.path.join(stage.workdir, input_list_name)
            if stage.inputlist != work_list_name:
                input_files = project_utilities.saferead(stage.inputlist)
                work_list = open(work_list_name, 'w')
                for input_file in input_files:
                    work_list.write('%s\n' % input_file.strip())
                work_list.close()

        # Now locate the fcl file on the fcl search path.

        fcl = project.get_fcl(stage.fclname)

        # Copy the fcl file to the work directory.

        workfcl = os.path.join(stage.workdir, os.path.basename(stage.fclname))
        if fcl != workfcl:
            shutil.copy(fcl, workfcl)

        # Construct a wrapper fcl file (called "wrapper.fcl") that will include
        # the original fcl, plus any overrides that are dynamically generated
        # in this script.

        wrapper_fcl = open(os.path.join(stage.workdir, 'wrapper.fcl'), 'w')
        wrapper_fcl.write('#include "%s"\n' % os.path.basename(stage.fclname))
        wrapper_fcl.write('\n')

        # Generate overrides for sam metadata fcl parameters.
        # Only do this if our xml file appears to contain sam metadata.

        xml_has_metadata = project.file_type != '' or \
                           project.run_type != '' or \
                           project.ubxml != ''
        if xml_has_metadata:

            # Add overrides for FileCatalogMetadata.

            if project.release_tag != '':
                wrapper_fcl.write('services.FileCatalogMetadata.applicationVersion: "%s"\n' % \
                                  project.release_tag)
            else:
                wrapper_fcl.write('services.FileCatalogMetadata.applicationVersion: "test"\n')
            if project.group:
                wrapper_fcl.write('services.FileCatalogMetadata.group: "%s"\n' % \
                                  project.group)
            if project.file_type:
                wrapper_fcl.write('services.FileCatalogMetadata.fileType: "%s"\n' % \
                                  project.file_type)

            # Add overrides for FileCatalogMetadataExtras.
        
            sep = 'services.user.FileCatalogMetadataExtras.Metadata:\n  ['

            # Add overrides for project (name, stage, version).
            # These are the main metadata attributes for identifying mc samples.
            
            wrapper_fcl.write('%s "ubProjectName", "%s"\n  ' % (sep, project.name))
            sep = ','
            if stage.name:
                wrapper_fcl.write('%s "ubProjectStage", "%s"\n  ' % (sep, stage.name))
            if project.ubxml:
                wrapper_fcl.write('%s "ubProjectVersion", "%s"\n  ' % (sep, project.ubxml))
            elif project.release_tag:
                wrapper_fcl.write('%s "ubProjectVersion", "%s"\n  ' % (sep, project.release_tag))

            # Run type.

            if project.run_type:
                wrapper_fcl.write('%s "runType", "%s"\n  ' % (sep, project.run_type))
                sep = ','

            # Other xml attributes may or may not represent sam metadata.
            # We include them only if we have added other overrides for
            # FileCatalogMetadataExtras.
            
            if sep == ',':
                wrapper_fcl.write('%s "fileFormat", "root"\n  ' % sep)
                if project.group:
                    wrapper_fcl.write('%s "group", "%s"\n  ' % (sep, project.group))
                wrapper_fcl.write('%s "fclName", "%s"\n  ' % (sep,
                                                              os.path.basename(stage.fclname)))
                wrapper_fcl.write('%s "fclVersion", "%s"\n  ' % (sep, project.ubfcl))
            wrapper_fcl.write(']\n')

        wrapper_fcl.close()

        # Copy and rename batch script to the work directory.

        workname = '%s-%s' % (stage.name, project.name)
        workname = workname + os.path.splitext(project.script)[1]
        workscript = os.path.join(stage.workdir, workname)
        if project.script != workscript:
            shutil.copy(project.script, workscript)

        # Copy and rename sam start project script to work directory.

        workstartscript = ''
        workstartname = ''
        if project.start_script != '':
            workstartname = 'start-%s' % workname
            workstartscript = os.path.join(stage.workdir, workstartname)
            if project.start_script != workstartscript:
                shutil.copy(project.start_script, workstartscript)

        # Copy and rename sam stop project script to work directory.

        workstopscript = ''
        workstopname = ''
        if project.stop_script != '':
            workstopname = 'stop-%s' % workname
            workstopscript = os.path.join(stage.workdir, workstopname)
            if project.stop_script != workstopscript:
                shutil.copy(project.stop_script, workstopscript)

        # Copy worker initialization script to work directory.

        if stage.init_script != '':
            if not os.path.exists(stage.init_script):
                print 'Worker initialization script %s does not exist.\n' % stage.init_script
                return 1
            work_init_script = os.path.join(stage.workdir, os.path.basename(stage.init_script))
            if stage.init_script != work_init_script:
                shutil.copy(stage.init_script, work_init_script)

        # Copy worker initialization source script to work directory.

        if stage.init_source != '':
            if not os.path.exists(stage.init_source):
                print 'Worker initialization source script %s does not exist.\n' % stage.init_source
                return 1
            work_init_source = os.path.join(stage.workdir, os.path.basename(stage.init_source))
            if stage.init_source != work_init_source:
                shutil.copy(stage.init_source, work_init_source)

        # Copy worker end-of-job script to work directory.

        if stage.end_script != '':
            if not os.path.exists(stage.end_script):
                print 'Worker end-of-job script %s does not exist.\n' % stage.end_script
                return 1
            work_end_script = os.path.join(stage.workdir, os.path.basename(stage.end_script))
            if stage.end_script != work_end_script:
                shutil.copy(stage.end_script, work_end_script)

        # If this is a makeup action, find list of missing
        # (cluster, process) pairs.
        # If sam information is present (cpids.list), create a makeup dataset.

        if makeup:
            
            missing_pairs = []
            missing_filename = os.path.join(stage.outdir, 'missing.txt')
            try:
                missing_files = project_utilities.saferead(missing_filename)
            except:
                print 'Cound not open %s.' % missing_filename
                return 1
            for line in missing_files:
                words = string.split(line)
                missing_pair = (int(words[0]), int(words[1]))
                missing_pairs.append(missing_pair)

                # Also delete the bad directory.

                missing_dir = os.path.join(stage.outdir, '%d_%d' % missing_pair)
                if os.path.exists(missing_dir):
                    shutil.rmtree(missing_dir)

            # Prepare sam-related makeup information.

            import_samweb()

            # Get list of successful consumer process ids.

            cpids = []
            cpids_filename = os.path.join(stage.outdir, 'cpids.list')
            try:
                cpids_files = project_utilities.saferead(cpids_filename)
                for line in cpids_files:
                    cpids.append(line.strip())
            except:
                pass

            # Create makeup dataset definition.

            makeup_defname = ''
            makeup_count = 0
            if len(cpids) > 0:
                makeup_defname = samweb.makeProjectName(stage.inputdef) + '_makeup'

                # Construct comma-separated list of consumer process ids.

                cpids_list = ''
                sep = ''
                for cpid in cpids:
                    cpids_list = cpids_list + '%s%s' % (sep, cpid)
                    sep = ','

                # Construct makeup dimension.
                
                dim = '(defname: %s) minus (consumer_process_id %s and consumed_status consumed)' % (stage.inputdef, cpids_list)

                # Create makeup dataset definition.

                print 'Creating makeup sam dataset definition %s' % makeup_defname
                samweb.createDefinition(defname=makeup_defname, dims=dim)
                makeup_count = samweb.countFiles(defname=makeup_defname)
                print 'Makeup dataset contains %d files.' % makeup_count

    # Do actions.

    rc = 0
    if check or checkana:

        # Check results from specified project stage.
        
        xml_has_metadata = project.file_type != '' or \
                           project.run_type != '' or \
                           project.ubxml != ''
        #input_files = get_input_files(stage)
        #print input_files
        has_input_files = stage.inputfile != '' or stage.inputlist != ''
        rc = docheck(stage.outdir, project.num_events, stage.num_jobs,
                     has_input_files, stage.inputdef,
                     checkana, xml_has_metadata, stage.histmerge)

    if check_definition or define:

        # Make sam dataset definition.

        if stage.defname == '':
            print 'No sam dataset definition name specified for this stage.'
            return 1
        dim = dimensions(project, stage)
        rc = docheck_definition(stage.defname, dim, define)

    if test_definition:

        # Print summary of files returned by dataset definition.

        if stage.defname == '':
            print 'No sam dataset definition name specified for this stage.'
            return 1
        rc = dotest_definition(stage.defname)

    if undefine:

        # Delete sam dataset definition.

        if stage.defname == '':
            print 'No sam dataset definition name specified for this stage.'
            return 1
        rc = doundefine(stage.defname)

    if check_declarations or declare:

        # Check sam declarations.

        rc = docheck_declarations(stage.outdir, declare)

    if test_declarations:

        # Print summary of declared files.

        dim = dimensions(project, stage)
        rc = dotest_declarations(dim)

    if check_locations or add_locations or clean_locations or remove_locations or upload:

        # Check sam disk locations.

        dim = dimensions(project, stage)
        rc = docheck_locations(dim, stage.outdir,
                               add_locations, clean_locations, remove_locations,
                               upload)

    if check_tape:

        # Check sam tape locations.

        dim = dimensions(project, stage)
        rc = docheck_tape(dim)

    if submit or makeup:

        # Sam stuff.

        # Get input sam dataset definition name.
        # Can be from xml or a makeup dataset that we just created.

        inputdef = stage.inputdef
        if makeup and makeup_defname != '':
            inputdef = makeup_defname

        # Sam project name.

        prjname = ''
        if inputdef != '':
            import_samweb()
            prjname = samweb.makeProjectName(inputdef)

        # Get proxy.

        proxy = project_utilities.get_proxy()

        # Construct jobsub command line for workers.

        command = ['jobsub']
        command_njobs = 1

        # Jobsub options.
        
        command.append('-q')       # Mail on error (only).
        command.append('--grid')
        command.append('--opportunistic')
        command.append('--group=%s' % project.group)
        if proxy != '':
            command.append('-x %s' % proxy)
        if project.os != '':
            command.append('--OS=%s' % project.os)
        if submit:
            command_njobs = stage.num_jobs
            command.extend(['-N', '%d' % command_njobs])
        elif makeup:
            if inputdef != '':
                command_njobs = min(makeup_count, stage.num_jobs)
                command.extend(['-N', '%d' % command_njobs])
                
            else:
                command_njobs = len(missing_pairs)
                command.extend(['-N', '%d' % command_njobs])

        # Batch script.

        command.append(workname)

        # Larsoft options.

        command.extend(['--group', project.group])
        command.extend(['-g'])
        #command.extend(['-c', os.path.basename(stage.fclname)])
        command.extend(['-c', 'wrapper.fcl'])
        if project.release_tag != '':
            command.extend(['-r', project.release_tag])
        command.extend(['-b', project.release_qual])
        if project.local_release_dir != '':
            command.extend(['--localdir', project.local_release_dir])
        if project.local_release_tar != '':
            command.extend(['--localtar', project.local_release_tar])
        if project.ubfcl != '':
            command.extend(['--ubfcl', project.ubfcl])
        command.extend(['--workdir', stage.workdir])
        command.extend(['--outdir', stage.outdir])
        if stage.inputfile != '':
            command.extend(['-s', stage.inputfile])
        elif input_list_name != '':
            command.extend(['-S', input_list_name])
        elif inputdef != '':
            command.extend(['--sam_defname', inputdef,
                            '--sam_project', prjname])
        command.extend(['-n', '%d' % project.num_events])
        command.extend(['--njobs', '%d' % stage.num_jobs ])
        if stage.init_script != '':
            command.extend(['--init-script',
                            os.path.join('.', os.path.basename(stage.init_script))])
        if stage.init_source != '':
            command.extend(['--init-source',
                            os.path.join('.', os.path.basename(stage.init_source))])
        if stage.end_script != '':
            command.extend(['--end-script',
                            os.path.join('.', os.path.basename(stage.end_script))])

        # If input is from sam, also construct a dag file.

        if prjname != '':

            # At this point, it is an error of the start and stop project
            # scripts were not found.

            if workstartname == '' or workstopname == '':
                print 'Sam start or stop project script not found.'
                sys.exit(1)

            # Start project jobsub command.
                
            start_command = ['jobsub']

            # General options.
            
            start_command.append('-q')       # Mail on error (only).
            start_command.append('--grid')
            start_command.append('--opportunistic')
            start_command.append('--group=%s' % project.group)

            # Start project script.

            start_command.append(workstartname)

            # Sam options.

            start_command.extend(['--sam_defname', inputdef,
                                  '--sam_project', prjname,
                                  '-g'])

            # Output directory.

            start_command.extend(['--outdir', stage.outdir])

            # Stop project jobsub command.
                
            stop_command = ['jobsub']

            # General options.
            
            stop_command.append('-q')       # Mail on error (only).
            stop_command.append('--grid')
            stop_command.append('--opportunistic')
            stop_command.append('--group=%s' % project.group)

            # Stop project script.

            stop_command.append(workstopname)

            # Sam options.

            stop_command.extend(['--sam_project', prjname,
                                 '-g'])

            # Output directory.

            stop_command.extend(['--outdir', stage.outdir])

            # Create dagNabbit.py configuration script in the work directory.

            dagfilepath = os.path.join(stage.workdir, 'submit.dag')
            dag = open(dagfilepath, 'w')

            # Write start section.

            dag.write('<serial>\n')
            first = True
            for word in start_command:
                if not first:
                    dag.write(' ')
                dag.write(word)
                if word == 'jobsub':
                    dag.write(' -n')
                first = False
            dag.write('\n</serial>\n')

            # Write main section.

            dag.write('<parallel>\n')
            for process in range(command_njobs):
                first = True
                skip = False
                for word in command:
                    if skip:
                        skip = False
                    else:
                        if word == '-N':
                            skip = True
                        else:
                            if not first:
                                dag.write(' ')
                            dag.write(word)
                            if word == 'jobsub':
                                dag.write(' -n')
                            first = False
                dag.write(' --process %d\n' % process)
            dag.write('</parallel>\n')

            # Write stop section.

            dag.write('<serial>\n')
            first = True
            for word in stop_command:
                if not first:
                    dag.write(' ')
                dag.write(word)
                if word == 'jobsub':
                    dag.write(' -n')
                first = False
            dag.write('\n</serial>\n')
            dag.close()

            # Update the main submission command to use dagNabbit.py instead of jobsub.

            command = ['dagNabbit.py', '-i', dagfilepath, '-s']

        os.chdir(stage.workdir)

        if submit:

            # For submit action, invoke the job submission command (jobsub or dagNabbit.py).

            subprocess.call(command)

        elif makeup:

            if inputdef == '':
                if len(missing_pairs) > 0:

                    # Construct missing process map file in workdir.

                    procmapname = 'procmap.txt'
                    procmapfile = open(procmapname, 'w')
                    for missing_pair in missing_pairs:
                        procmapfile.write('%d\n' % missing_pair[1])
                    procmapfile.close()

                    # Update command.
                
                    new_command = command
                    #new_command.extend(['--cluster', '%d' % missing_pairs[0][0]])
                    new_command.extend(['--procmap', procmapname])
                    subprocess.call(new_command)

            else:
                if makeup_count > 0:
                    subprocess.call(command)
                    
        
    # Done.

    return rc

# Invoke main program.

if __name__ == '__main__':
    sys.exit(main(sys.argv))
