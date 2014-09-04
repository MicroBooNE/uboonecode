#!/usr/bin/env python
import sys, getopt
import os
import subprocess
from subprocess import Popen, PIPE
import time
import samweb_cli
#from samweb_client.utility import fileEnstoreChecksum
import ast
import project_utilities, root_metadata

def getmetadata(inputfile):
	# Set up the experiment name for samweb Python API	
	samweb = samweb_cli.SAMWebClient(experiment=project_utilities.get_experiment())

	# Extract metadata into a pipe.
	local = project_utilities.path_to_local(inputfile)
	if local != '':
		proc = subprocess.Popen(["sam_metadata_dumper", local], stdout=subprocess.PIPE)
	else:
		url = project_utilities.path_to_url(inputfile)
		proc = subprocess.Popen(["sam_metadata_dumper", url], stdout=subprocess.PIPE)
	lines = proc.stdout.readlines()
	if local != '' and local != inputfile:
		os.remove(local)

	# Count the number of lines in the file (for later use!)
	num_lines = len(lines)

	# define an empty python dictionary
	md = {}	

	#Read tbe columns from the file and fill the dictionary
	c = 0
	p = 0
	parents = []
	PName = False
	gen = False
	for line in lines:
		c = c+1
		columns = line.split(" ")
		columns = [col.strip() for col in columns]	
		if c>=4 and c<=num_lines-2:
			if columns[1] == 'dataTier':
				md['data_tier'] = columns[-1]
				if columns[-1] == 'generated':
					gen = True
		       	elif columns[1] == 'endTime':
				E  = time.localtime(int(columns[-1]))		
				md['end_time'] = str(E[0])+'-'+str(E[1])+'-'+str(E[2])+'T'+str(E[3])+':'+str(E[4])+':'+str(E[5])
			elif columns[1] == 'startTime':
				S  = time.localtime(int(columns[-1]))
				md['start_time'] = str(S[0])+'-'+str(S[1])+'-'+str(S[2])+'T'+str(S[3])+':'+str(S[4])+':'+str(S[5])
			elif columns[1] == 'group':
				md['group']  = columns[-1]
			elif columns[1] == 'eventCount':
				md['event_count']  = columns[-1]
			elif columns[1] == 'fclName':
				md['fcl.name']  = columns[-1]
			elif columns[1] == 'fclVersion':
				md['fcl.version']  = columns[-1]
			elif columns[1] == 'fileFormat':
				md['file_format']  = columns[-1]
			elif columns[1] == 'ubProjectStage':
				md['ub_project.stage']  = columns[-1]
			elif columns[1] == 'ubProjectVersion':
				md['ub_project.version']  = columns[-1]
			elif columns[1] == 'lastEvent':
				md['last_event']  = columns[-1]
			elif columns[1] == 'firstEvent':
				md['first_event']  = columns[-1]
			elif columns[1] == 'fileType':
				md['file_type']  = columns[-1]
			elif columns[1] == 'group':
				md['group']  = columns[-1]
			elif columns[1] == 'group':
				md['group']  = columns[-1]
			elif columns[1] == 'run':
				run = columns[-1]
			elif columns[1] == 'runType':
				run_type = columns[-1]
			elif columns[1] == 'applicationFamily':
				app_family = columns[-1]
			elif columns[1] == 'applicationVersion':
				app_version = columns[-1]
			elif columns[1] == 'process_name':
				app_name = columns[-1]
			elif columns[1] == 'ubProjectName':
				PName = True
				md['ub_project.name'] = columns[-1]
			elif columns[1] == 'parent':
				parents.append({'file_name': columns[-1]})

	# Get the other meta data field parameters						
	md['file_name'] =  inputfile.split("/")[-1]
	md['file_size'] =  os.path.getsize(inputfile)
	# For now, skip the checksum for dCache files.
	md['crc'] = root_metadata.fileEnstoreChecksum(inputfile)
	md['runs']      =  [[run, run_type]]
	md['application'] = {'family': app_family, 'name': app_name, 'version': app_version}
	md['parents'] = parents

	# If ub_project.name is not in the internal metadata,
	# for generator files, get the ub_project.name from the fcl_filename (without the '.fcl' part) for gen files.
	# for all other stages, get this from the parents
	if gen == True:
		md['parents'] = []
		if PName == False:
			md['ub_project.name'] = md['fcl.name'].split(".fcl")[0]
	else:
		if PName == False:
			if md.has_key('parents'):
				parent = md['parents'][0]['file_name']
				mdparent = samweb.getMetadata(parent)
				if mdparent.has_key('ub_project.name'):
					md['ub_project.name'] = mdparent['ub_project.name']

	return md

if __name__ == "__main__":
	md = getmetadata(str(sys.argv[1]))
	#print md	
	mdtext = samweb_cli.json.dumps(md, sys.stdout, indent=2, sort_keys=True)
	print mdtext
	sys.exit(0)	
