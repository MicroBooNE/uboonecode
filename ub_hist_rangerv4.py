#Evan Peters 20150920
#v4: Adapted to take Gabe's new rangefile format: 'leafname subleafnum min max'; 0 = leaf
#v3: Generalizing OptionParser for N files
#v2: Finalizing execute(); methods now populate 'master' dictionaries
#A module for dealing reranging histograms from uB
#Adapted from Gabe Nowak's uB_v2.py
#Haha, okay, really it was stolen from Gabe
import sys, os 
from optparse import OptionParser
import array
import time #for timing
#print sys.version

#plotting things
import numpy as np
import array
import matplotlib.pyplot as plt
import matplotlib as mpl

#load/set ROOT options if necessary
import ROOT as R
from ROOT import gDirectory
R.gROOT.Reset()



class datadatadata:

	def __init__(self, script=[], quiet=False):	
		
		self.global_quiet = quiet
		self.script = script
		
		#Debug variables
		self.debug = True
		self.debug_length = 100000
		self.test_dct = {'a':{'aa':{'aaa':0, 'aab':0, 'aac': [0, 0, 1]}, 'ab':{'aba':0, 'abb':{'x':'?'}, 'abc': [0, 0, 1]}}, 'b':{'ba':{'baa':0, 'bab':0, 'bac': [0, 0, 1]}, 'bb':{'bba':0, 'bbb':{'x':'?'}, 'bbc': [0, 0, 1]}}, 'c':{'cc':{'ccc':{'cccc':'5tier'}}}}
		
		#Parser Parameters
		if not self.debug:
			self.nfiles = int(raw_input("Enter the number of file sets you gave as input:\n"))
		else:
			self.nfiles = 2
			
		self.data_list = []
		self.hout = None
		
		#Data management
		self.tree_dct = {'tree number': {'size': None, 'handle': None} }
		self.tree_dct = {}
		self.tree_count = 0

		
		#Parameters and dictionary for range_finder
		self.rfile = None
		self.rfile_finished = False
		self.bigleaf_master_dct = {}
		self.ranges_found = False
		
		
		#length of mini-arrays for calculating maxima of leaves
		self.opt_len = 10000
		
		#Dictionary of branch dct properties (format in file_organizer())
		self.br_hist_master_dct = {}
		self.master_ranges_found = False
		
		#Dictionaries of the most extreme ranges for all subleafs (leaf_master_ranges), and leafs (br_master_ranges)
		self.leaf_master_ranges = {}
		self.br_master_ranges = {}
		#Read in files
		self.read_in()

		
		#dct_printer vars
		self.spacing = 0
		
		#file_explorer vars
		self.explorer_dct = {}
		
		#Data Description variables
		self.NUM_CHANNELS = 8256 #The expected number of channels when constructing profile histograms
		
		
		# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
		#EXECUTE REGION	
		if "debug" in self.script:
			self.debug = True
			self.execute()
			#self.dct_printer(self.br_hist_master_dct)
			#self.dct_printer(self.br_hist_dct)
			#self.spacing = 0
		elif "test" in self.script:
			self.dct_printer(self.file_explorer(self.data_list[0], self.explorer_dct, get_leafs=False))
			
			
			#self.dct_printer(self.file_explorer(self.data1, self.explorer_dct, get_leafs=False))
			#self.spacing = 0

		# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
		
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	def read_in(self,  ):
	#Creates option parser, imports ROOT objects for queries
	
		parser = OptionParser()
		#nargs turns this option into a tuple of input files
		parser.add_option("--input",dest="input",help="Ordered input files, .root only",type='string', nargs=self.nfiles)
		parser.add_option("--range",dest = "rfile", help = "Range file input. Format:'branch min max")
		parser.add_option("--profile",dest = "profile", help = "Produce profile histogram for self.NUM_CHANNELS-channel Leafs (Default True)", default = True, action = "store_false")
		parser.add_option("--weight",dest = "wgt",help = "Use weights when drawing from the corresponding Rfile. Entry must be a list of bools of equal length to list of files. The TTree must have a weight branch.", default = [False]*self.nfiles, nargs=self.nfiles)
		#Assign options and args variables
		(self.options,self.args)=parser.parse_args()
		
		#Open input files into self.data_list
		for file_name in self.options.input:
			self.data_list.append(R.TFile.Open(file_name,'r'))
		
		#rfile management
		if self.options.rfile and not os.path.isfile(self.options.rfile):
			self.safe_print("Error at readin(): The range file %s does not exist." % options.rfile)
		self.rfile = self.options.rfile
	
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
	def __file_checker(self, index, quiet=False):
	#Browses the input file to check for proper formatting
	
		#Inputs must exist
		try:
			main = self.data_list[index]
		except: 
			"Error at file_organizer(%s): Could not read in File %s" % (index, index)	
			sys.exit()
			return
		#Inputs must not have subdirectories
		try:
			bad = main.GetListOfKeys()[1]
		except:
			pass
		else:
			print "Error: More than one directory in input file"
			sys.exit()
			return	
		
		#If weights are enabled, the current input tree needs a weight branch
		if self.options.wgt[index]:
			for key in main.GetListOfKeys():
				folder = main.Get(key.GetName())
			for key in folder.GetListOfKeys():
				tree = folder.Get(key.GetName())
				#print str(type(tree))
				#Weight branches are either named 'wgt' or weight
				if str(type(tree)) == "<class 'ROOT.TTree'>":
					for name in ["wgt", "weight"]:
						if name not in tree.GetListOfBranches():
							self.safe_print("Error at file_organizer(%s): TTree %s Does not have a wgt branch" % \
							(index, tree.GetName()), quiet)	
							sys.exit()
							return
							
		#All tests passed, continue on
		return main
		
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -			
	def __name_creator(self, index):
	#Given a file set, creates an appropriate .root filename with output histograms
		
		file_name = self.options.input[index]
		file_name.replace(".root", "")		
		file_name += "_histograms"
		
		if self.options.wgt[index]:
			file_name += "_weighted"
		if self.rfile:
			rfile_name = self.rfile.replace(".txt", "")
			file_name += "(%s)" % rfile_name
		file_name += ".root"
		return file_name
			
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -	-	
	#!!ADD AUTO PARAMETER?
	def file_organizer(self, index, quiet=False ):
	#Populates br_hist_dct - a dictionary of branches (leaves) and histograms
		self.safe_print("Running file_organizer(%s)" % index, quiet)
		
		main = self.__file_checker(index, quiet=quiet)
		br_hist_dct = {'item': {"handle": None, "range": [None, None], "entries": 0} }
		br_hist_dct = {}
		
		subdir = main.GetListOfKeys()[0]
		subdir_name = subdir.GetName()
		subdir_obj = main.Get(subdir_name)
		for obj in subdir_obj.GetListOfKeys():
			obj_type = obj.GetClassName()
			obj_name = obj.GetName()
			#Initialize subdictionary for the tree corresponding to the given file
			
			if str(obj_type) == "TTree":
				tree = subdir_obj.Get(obj_name)
				#tree_dct contains lists at the bottom level - accounts for files with many trees
				self.tree_dct[index] = {"handle": [], "size": []}
				self.tree_dct[index]["handle"].append(tree)
				self.tree_dct[index]["size"].append(tree.GetEntries())
				obj_iter = subdir_obj.Get(obj_name).GetListOfBranches()
				for branch in obj_iter:
					br_name = branch.GetName()
					#List index corresponds to the tree the branch belongs to
					br_hist_dct[br_name] = {"handle": []}
					#Send INSTANCE to br_hist_dct (does not retain methods of pointer(?)
					br_hist_dct[br_name]["handle"].append(tree.GetBranch(br_name))
			elif str(obj_type) == "TH1":
				br_hist_dct[obj_name] = {}
				br_hist_dct[obj_name]["handle"] = subdir_obj.Get(obj_name)
			#Check for subdirectory
			elif str(obj_type) == "TDirectoryFile":
				print "Error: Subdirectory found in main directory"
				return
			
		#Populate br_hist_dct with ranges, write out to master dct
		self.br_hist_master_dct[index] = br_hist_dct
		#Check that ANY sort of range was given
		if self.rfile:
			self.safe_print("File set %s organized" % index, quiet)
		else:
			print "Fatal error at file_organizer(%i): No rfile provided" % index
			sys.exit()
		
		#Handling output file
		hout_name = self.__name_creator(index)
		if self.debug:
			self.hout = R.TFile("test%s.root" % index,"RECREATE")
		else:
			self.hout = R.TFile(hout_name, "RECREATE")
		self.safe_print("Writing to %s" % self.hout, quiet)
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -	
	def dct_initializer(self, index):
	#For future reference, this would be a good idea ...
		return
		
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -		
	def safe_print(self, string, quiet):
		if not quiet and not self.global_quiet:
			print string
		
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -	
	#PENDING MORE INFO ON FUTURE INPUT FILE FORMATS
	def range_reader(self, index, quiet=False):
	#Read in text file and populate br_hist_master_dct
			
		self.safe_print("Running range_reader(%s)" % index, quiet)
		#NOTE - ONLY SUPPORTS A SINGLE RFILE, AND COPIES THOSE RANGES TO ALL INPUT FILE BRANCHES
		if self.rfile_finished:
			self.br_hist_master_dct[index] = self.br_hist_master_dct[0].copy()
			return
					
		try:
			rfile = open(self.rfile, 'r')
		except:
			print "Error at range_reader(%s): bad rfile provided for file set" % index
			return
		
		#Initialize br_hist_master_dct using leafs counted in leaf_range_finder(index)
		ran = False
		for i_tree, tree in enumerate(self.tree_dct[index]["handle"]):
			if not ran:
				for br_name in self.bigleaf_master_dct[index]:
					self.br_hist_master_dct[index][br_name] = {}
					for i in xrange(self.bigleaf_master_dct[index][br_name][i_tree]):
						#An entry for every subleaf (a single entry at 0 for no subleafs)
						self.br_hist_master_dct[index][br_name][i] = []
				ran = True
			elif ran:
				#Check for consistency with all other trees in file
				#Currently throws errors if not all trees contain the same branches
				for branch in tree.GetListOfBranches():
					br_name = branch.GetName()
					if br_name not in self.br_hist_master_dct[index].keys():
						print "range_reader(%s) Error: input(%s) trees contain different branches" % (index, index)
						print "Ranges cannot be read in this format"
						sys.exit()
						return

		#NOTE - CURRENTLY USES THE SAME BRANCH RANGES FOR EVERY TREE IN THE FILESET		 
		#Append a range for the branch
		
		for i_tree, tree in enumerate(self.tree_dct[index]["handle"]):
			rfile.seek(0)
			for line_i, line in enumerate(rfile):
				#We need to make sure the rfile is in the correct format
				try:
					line_list = line.split()
					input_br_name = line_list[0]
					leaf_num = int(line_list[1])
				except:
					print "Fatal Error at range_reader(%s): Improper format of rfile"
					sys.exit()
					return
				if input_br_name not in self.br_hist_master_dct[index]:
					print "Warning at range_reader(%s): branch name given in rfile not found, skipping line %i" % (index, line_i)
					continue
				self.br_hist_master_dct[index][input_br_name][leaf_num].append([line_list[2], line_list[3]])

		self.rfile_finished = True	
		self.safe_print("File set %s branch ranges populated from input file" % index, quiet)
			
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -		
	def leaf_range_finder(self, index, fast=False, quiet=False):
	#Finds ranges of BIG leafs with subleafs via histograming method
	#Option fast determines whether every leaf will be checked for having a consistent number of subleafs over all events
	
		self.safe_print("Running leaf_range_finder(%s)" % index, quiet)
		#Initialize dictionary of branchname: list of leafsize
		self.bigleaf_master_dct[index] = {}
		
		#Initialize bigleaf_master_dct: dictionary of N_subleafs for all files, trees
		#Branches were already checked for consistency in range_reader()
		temp_tree = self.tree_dct[index]["handle"][0]
		for branch in temp_tree.GetListOfBranches():
			br_name = branch.GetName()
			self.bigleaf_master_dct[index][br_name] = []
				
		#Iterate over all trees read from the given file
		for i, tree in enumerate(self.tree_dct[index]["handle"]):
			tree = self.tree_dct[index]["handle"][i]
			#Populate big_leaf_dct with dimensions of each leaf
			for branch in tree.GetListOfBranches():
				br_name = branch.GetName()
				N_sl = self.__N_subleaf_counter(tree, br_name, fast=fast, quiet=True)
				#The bottom layer of bigleaf_master_dct is a list w index corresponding to tree
				self.bigleaf_master_dct[index][br_name].append(N_sl)

		self.safe_print("Leaf ranges found for file set %s" % index, quiet)
		return self.bigleaf_master_dct
		
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -		
	def __N_subleaf_counter(self, tree, br_name, fast=False, quiet=False):
	#tree is the tree containing the branch with GetName() of br_name, containing the desired leaf
	#Throws errors for leafs that have a variable number of subleafs; otherwise returns the number of subleafs in a leaf
	#Option fast determines whether every leaf will be checked for having a consistent number of subleafs over all events
		self.safe_print("Running N_subleaf_counter(fast=%s)" % fast	, quiet)	
		#Give the number of subleafs for the first event, in this branch
		if fast:
			tree.GetEntry(0)
			return tree.GetLeaf(br_name).GetLen()
		
		tree_size = tree.GetEntries()
		N_subleaf_list = [0]*tree_size
		#Find out the subleaf population for every event	
		for event in xrange(tree_size):
			tree.GetEntry(event)
			leaf = tree.GetLeaf(br_name)
			N_subleaf_list[event] = leaf.GetLen()
			
		#Check that the number of subleafs in the given branch doesn't vary over events
		unq_N =  list(set(N_subleaf_list))
		variation = len(unq_N)
		if variation == 1:
			return unq_N[0]
		elif variation > 1:
			print "N_subleaf_counter(%s, %s) Error: %s.Len() varies over events" (tree, br_name, br_name)
			if variation < 20:
				print "Number of subleafs varies among: ", unq_N
			else:
				print "Over 20 different values found for the number of subleafs"
			return None
		
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	def dct_printer(self, obj, started=False, trim=True):
	#Prints a dictionary's contents, with some formatting for readability
	#Do not alter option started
		if not started:
			self.spacing = 0
		if type(obj) == dict:
			self.spacing += 2
			for count, key in enumerate(obj):
				try:
					obj[key].iteritems()
				#Print out lines for bottom layer of dct
				except:
					if count > 20 and trim:
						break
					print "".join([" "]*self.spacing),"%s: %s" % (key, obj[key])

				else:
					print "".join([" "]*self.spacing), key
					self.dct_printer(obj[key], started=True, trim=trim)
					self.spacing -= 2
		else:
			return
			
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -	
	#currently broken...
	def file_explorer(self, obj_handle, dct, get_leafs = True):	
	#Given a file handle, returns a dictionary corresponding to the file structure
	#get_leaves determines whether leaves are included in the structure (often redundant w Branch)
		obj_type = str(type(obj_handle)) 
		print obj_type
		obj_name = obj_handle.GetName()
		dct[obj_name] = {}
		
		#Choose iterator and method of calling subfile based on obj type
		#Recursively constructs dictionaries if obj is an iterator
		if obj_type in ["<class 'ROOT.TFile'>", "<class '__main__.TKey'>", "<class 'ROOT.TDirectoryFile'>"]:
			obj_iter = obj_handle.GetListOfKeys()
			for sub in obj_iter:
				sub = obj_handle.Get(sub.GetName())
				self.file_explorer(sub, dct[obj_name], get_leafs=get_leafs)

		elif obj_type == "<class 'ROOT.TTree'>":
			obj_iter = obj_handle.GetListOfBranches()
			for sub in obj_iter:
				sub = obj_handle.GetBranch(sub.GetName())
				self.file_explorer(sub, dct[obj_name], get_leafs=get_leafs)
				
		elif obj_type == "<class '__main__.TBranch'>" and get_leafs == True:
			print "!!!"
			obj_iter = obj_handle.GetListOfLeaves()
			for sub in obj_iter:
				sub = obj_handle.GetLeaf(sub.GetName())
				self.file_explorer(sub, dct[obj_name], get_leafs=get_leafs)
		
		return dct
		
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -	
	def master_range_finder(self, index, quiet=False):
	#Parses over all trees in all inputs and populates leaf_master_ranges with the range extrema
		
		#Only find master ranges once
		if self.master_ranges_found:
			return
		self.safe_print("Running master_range_finder()", quiet)
		
		#Initialize leaf_master_ranges, and br_master_rangesincluding branches from every file input
		for fileset in self.br_hist_master_dct:
			for br_name in self.br_hist_master_dct[fileset]:
				if br_name in self.leaf_master_ranges:
					continue
				#br_master_ranges only have the maxima for an entire 'leaf' (no subleafs)
				self.leaf_master_ranges[br_name] = {}
				for leaf_num in self.br_hist_master_dct[fileset][br_name]:
					#Put a list of sublists at each leaf, we will repopulate with the min/max of these sublists
					self.leaf_master_ranges[br_name][leaf_num] = [[],[]]
				#We use placeholders for later comparison to leaf_master_ranges running max
				self.br_master_ranges[br_name] = [[1e10], [-1e10]]

		
		#Populate leaf_master_ranges			
		for fileset in self.br_hist_master_dct:
			for br_name in self.br_hist_master_dct[fileset]:
				for leaf_num in self.br_hist_master_dct[fileset][br_name]:
					file_min_max_list = [[], []]
					#Find the range extrema for trees within this single file...
					for i_tree, tree in enumerate(self.tree_dct):
						file_min_max_list[0].append(self.br_hist_master_dct[fileset][br_name][leaf_num][i_tree][0])
						file_min_max_list[1].append(self.br_hist_master_dct[fileset][br_name][leaf_num][i_tree][1])				
					file_min_max = [min(file_min_max_list[0]), max(file_min_max_list[1])]
					#Append the most extreme extrema for this file to the master_range
					self.leaf_master_ranges[br_name][leaf_num][0].append(file_min_max[0])
					self.leaf_master_ranges[br_name][leaf_num][1].append(file_min_max[1])
		
		#Replace entry in each leaf with most extreme range
		#Also populate br_master_ranges with running best range
		for br_name in self.leaf_master_ranges:
			for leaf_num in self.leaf_master_ranges[br_name]:
				self.leaf_master_ranges[br_name][leaf_num] = \
				[min(self.leaf_master_ranges[br_name][leaf_num][0]),max(self.leaf_master_ranges[br_name][leaf_num][1])]
				if (self.leaf_master_ranges[br_name][leaf_num][0] < self.br_master_ranges[br_name][0]):
					self.br_master_ranges[br_name][0] = self.leaf_master_ranges[br_name][leaf_num][0] 
				if (self.leaf_master_ranges[br_name][leaf_num][1] > self.br_master_ranges[br_name][1]):
					self.br_master_ranges[br_name][1] = self.leaf_master_ranges[br_name][leaf_num][1]
					
		self.master_ranges_found = True
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -	
	def hist_manager(self, index, tree_index_list = [None], quiet=True):
	#Creates histogram of branches for trees in the current file
	#tree_index_list is a list of the positions of the desired trees in the file [0, 1, ...]
		self.safe_print("Runnning hist_manager(%s)" % index, quiet)
		
		#By default, make histograms for all trees in the current fileset
		if tree_index_list == [None]:
			tree_list = self.tree_dct[index]["handle"]
		else:
			tree_list = [self.tree_dct[index]["handle"][i] for i in tree_index_list]
		
		#Currently NOT iterating over all trees, just all branches found in all files
		#Run different histogram routines depending on the type of branch
		for i_tree, tree in enumerate(tree_list):
			for br_name in self.br_master_ranges:
				N_subleafs = self.bigleaf_master_dct[index][br_name][i_tree]
		
				#For Leafs with no subleafs
				if N_subleafs == 1:
					continue
					string = "(100, %s, %s)" % \
					(self.br_master_ranges[br_name][0], self.br_master_ranges[br_name][1])
					#Setting weights
					if self.options.wgt[index]:	
						string = "%s >> h1%s wgt * (%s)" % (br_name, string, br_name)
						tree.Draw(string)
					elif not self.options.wgt[index]:
						string = "%s >> h1%s" % (br_name, string)
						tree.Draw(string)
					htemp = R.gPad.GetPrimitive("h1")
					#Currently naming based on which TTree is being checked
					htemp.SetName("TTree%i:%s" % (i_tree+1, br_name))
					htemp.SetTitle("TTree%i:%s" % (i_tree+1, br_name))
					#Checking for empty histrogram
					integral = htemp.Integral()
					if integral == 0:
						self.safe_print("histogram_manager(%s) Warning: %s is empty, may need regeneration" % (index, br_name), quiet)
						continue
					else:
						if i_tree == 0:
							self.hout.mkdir(br_name)
						self.hout.cd(br_name)
						htemp.Write()
						self.safe_print("File %s:TTree %s: TBranch %s was written" % \
						(index, i_tree, br_name), quiet)
									
				#Profile for specific Leafs with self.NUM_CHANNELS Channels
				if N_subleafs == self.NUM_CHANNELS and self.options.profile == True:
					if i_tree == 0:
						self.hout.mkdir(br_name)
					self.hout.cd(br_name)
					string = "%s:%s >> h1(%s, 0, %s)" % (br_name, N_subleafs, N_subleafs, N_subleafs)
					print string
					print tree
					
					tree.Draw(string, "", "profs")
					profile = R.gPad.GetPrimitive("h1")
					profile.SetName("TTree%i:%s Profile" % (i_tree,br_name))
					profile.SetTitle("%s Profile" % br_name)
					profile.Write()
					

			
					#Construct subleaf histograms (over all trees, etc)
					for i in xrange(N_subleafs):
						n_s = "TTree%i:%s_%i" % (i_tree, br_name, i)
						string = "(100, %s, %s)" % (self.leaf_master_ranges[br_name][i][0], self.leaf_master_ranges[br_name][i][1])
						if self.options.wgt[index]:
							string = "%s[%i] >> h1%s wgt * (%s[%i])" % (br_name, i, string, br_name, i)
							tree.Draw(string)
						elif not self.options.wgt[index]:
							string = "%s[%i] >> h1%s" % (br_name, i, string)
							tree.Draw(string)
							htemp = R.gPad.GetPrimitive("h1")
						integral = htemp.Integral()
						if integral == 0:
							self.safe_print("histogram_manager(%s) Warning: %s[%i] is empty, may need regeneration" % (index, br_name, i), quiet)
							continue
						else:
							
							htemp.SetName(n_s)
							htemp.SetTitle(n_s)
							integral = htemp.Integral()
							htemp.Write()
							self.safe_print("File %s:TTree %s: Subleaf %s[%i] was written" % \
							(index, i_tree, br_name, i), quiet)	
						#DEBUG
						if self.debug and i > 20:
							break	
		self.hout.Close()	
					
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -		
	def execute(self):
	#Executes methods consistent with user inputs
	
		#Prevent canvases from popping up		
		R.gROOT.SetBatch(True)  
		
		#Populate master dictionary with entries for each data input
		for i, data in enumerate(self.data_list):
			self.file_organizer(i)
			self.leaf_range_finder(i, fast=False)
			self.range_reader(i)
			self.master_range_finder(i)
			self.hist_manager(i)

					
		#sys.exit()
		
		

#Plotting module
run = datadatadata(["debug"])		
		
		
		
'''
				
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -		
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -	
	def range_finder(self, auto=True, obj=None, quiet=False):
	#Takes a SINGLE branch instance as input
	#Creates a dictionary of TTree branches and corresponding [minimum, maximum]
	#If auto, method runs on all files in self.file_dct
		if not quiet and not self.global_quiet:
			print "range_finder initializing"
			
		if auto:
			obj = self.br_hist_dct
		elif not auto and obj != None:
			obj = {obj.GetName(): obj}
		else:
			print "Error: Method was not passed valid ROOT object"
			return
		
		if self.ranges_found:
			print "br_hist_dct already populated; range_finder cancelled"
			return
			
		#Approach: Create Python arrays ('minilists') of a given length, use built-ins to find max/min
		
		x = 0
		
		t0 = time.time()
		#Initialize mini-lists for each branch
		for branch in self.tree.GetListOfBranches():
			br_name = branch.GetName()
			self.br_hist_dct[br_name]["minilist"] = {}
			leaf = self.tree.GetLeaf(br_name)
			#Create dictionary of leaflet minilists for each branch
			for k in xrange(leaf.GetNdata()):
				self.br_hist_dct[br_name]["minilist"][k] = [0]*self.opt_len
				
			#Parse through leafs and subleafs to create mini-arrays to find minima, maxima
			for p, i in enumerate(xrange(x, x+self.opt_len)):
				if i >= self.tree_size:
					break
					self.tree.GetEntry(i)
					for k in xrange(leaf.GetNdata()):
						#populate minilist with appropriate leaflet, at appropriate position
						self.br_hist_dct[br_name]["minilist"][k][p] = leaf.GetValue(k)
		
		print "Time: ", time.time()-t0
		self.ranges_found = True
		return

		#IN PROGRESS; NEED MORE EFFICIENT METHOD OF READING IN LEAFS
			
		
		
		
		#Parse through branches in br_hist_dct; item is either hist or branch
		for key, subdct in obj.iteritems():
			item = obj[key]["handle"]
			item_type = str(type(item))
			item_name = key
			if item_type == "<class '__main__.TBranch'>":
				br_name = item_name
				for i in range(10):
					self.tree.GetEntry(i)
					leaf = self.tree.GetLeaf(br_name)

				print br_name
				for i in range(10):
					print leaf.GetValue(i)
					print leaf.GetMaximum()
					#print "  ",leaf.GetValue()
					#print "  ",leaf.GetLen()
					#self.br_hist_dct[br_name]["range"] = [item.GetMinimum(br_name), item.GetMaximum(br_name)]
					
							
			#Trying to histogram branches and return empty bins...		
			elif item_type == "hist?":
				hist_name = item_name
				pass
			
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -	

#Current tasks:
#Write out dictionaries to external text file
#~OR~ keep i==100000 and delete trailing 0's on vectors

#Fix TNtuple writing script:
		#Write out tuples
		#print self.hout




#This is all scratch, some of it useful
def scratch():
	#Auto-range using histogram function - use it for writing rfiles, if necessary.
	for br_name in self.bigleaf_master_dct[index]:
		for i_tree, tree in enumerate(self.tree_dct[index]["handle"]):
			N_sl = self.bigleaf_master_dct[index][br_name][i_tree]
			#Only do histogram ranging for BIG leafs
			#Double check that this branch wasn't assigned values in range_reader
			if N_sl == self.NUM_CHANNELS and self.br_hist_master_dct[index][br_name]["range"][i_tree] == [None, None]:
				if not quiet and not self.global_quiet:
					print "Auto-ranging file%s:tree%s:branch: %s" % (index, i_tree, br_name)
				tree.Draw("%s: %s >> h1(%s, 0, %s)" % (br_name, N_sl, N_sl, N_sl), "", "profs")
				h1 = R.gPad.GetPrimitive("h1")
				self.br_hist_master_dct[index][br_name]["range"][i_tree] =[h1.GetMinimum(), h1.GetMaximum()]
		else:
			continue


	#Scans trees for empty leaves, yes/no weights, 
		bottom = False
		prev = main
		current = main
		for sub in main.GetListOfKeys():
			current = sub
			current_type = current.GetClassName()
			current_obj = prev.Get(sub.GetName())
			current_name = current.GetName()
			
			print current, current_type, current_obj, current_name
			print "Current Folder:", current_name
			
			while bottom == False:
				try:
					current_iter = current.GetListOfKeys()
				except: 
					print "Bottom Found"
					bottom = True
				else:
					for subsub in current_iter:
						current = subsub

				
		if dim == 1:
			if br not in ranged3.keys():
				continue
			string = "(100,%s,%s)" % (ranged3[br][0],ranged3[br][1])
			if eval("options.wgt" + str(num)) == True:
				ntuple.Draw(br + " >> h1" + string,"wgt * (" + br + ")")
			else:
				ntuple.Draw(br + " >> h1" + string)
			htemp = gPad.GetPrimitive("h1")
			htemp.SetName(br)
			htemp.SetTitle(br)
			integral = htemp.Integral()
			if integral == 0:
				print br + " has no data points within the dictionaries bounds; may need to regenerate the dicitonary!"
				continue
			if options.exclude == True:
				if htemp.GetMinimum() == htemp.GetMaximum():
					continue
			else:
				htemp.Write()
				print br + " has been written."					
					
				#hout.mkdir(br)
				#hout.cd(br)
'''	
