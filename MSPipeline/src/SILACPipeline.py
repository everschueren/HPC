#!/usr/bin/env python
# encoding: utf-8
"""
SILACPipeline.py

Created by erik verschueren on 2013-07-06.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

import os
import sys
import getopt
import ConfigParser
import subprocess


help_message = '''
The help message goes here.
'''

MS_PIPELINE_PATH = "/Users/everschueren/Projects/HPCKrogan/Scripts/MSPipeline/"
src_dir = MS_PIPELINE_PATH + "/src/"
bin_dir = MS_PIPELINE_PATH + "/bin/"
os.environ['MS_PIPELINE_PATH'] = MS_PIPELINE_PATH


def readConfig(filename):
	file = open(filename)
	config = ConfigParser.RawConfigParser()
	config.readfp(file)
	
	return config

def pipeline(config):
	print(">> KROGAN PTM PIPELINE")

	## read files
	dir = config.get("files","dir")	
	data_file = dir+config.get("files","data")
	keys_file = dir+config.get("files", "keys")	

	##derived
	bname = os.path.splitext(os.path.basename(data_file))[0]
	tmp_in_file = data_file
	
	#################################
	## PREPARE CONTAINERS FOR RESULTS

	output_dir = dir + config.get("files","output_dir") + '/'
	if not os.path.exists(output_dir):
		os.mkdir(output_dir)

	# out_commands = open(output_dir+"commands.sh",'w')	
	
	## ################################ 
	## 1. read maxquant data and filter
	
	oxidations = False
	unique = False
	contaminants = False
	modification = None

	if config.getboolean("filters","enabled"):
		print(">> FILTERING DATA")
		bname = bname+"_FLT"
		tmp_out_file = output_dir+bname+".txt"
		oxidations = config.get("filters","oxidations")
		unique = config.get("filters","unique")
		contaminants = config.get("filters","contaminants")
		modification = config.get("filters","modification")

		subprocess.call([src_dir+'Conversion/MaxQFilter.R', '-d', tmp_in_file, '-o', tmp_out_file, '--contaminants_filter', contaminants, '-u', unique, '-x', oxidations, '-s', modification])
		tmp_in_file = tmp_out_file

	## #############################################
	## 2. convert to matrix of normalized log-ratios
	
	bname = bname+"_MAT"
	tmp_out_file = output_dir+bname+".txt"

	if config.get("normalization", "enabled"):
		print(">> NORMALIZING BETWEEN REPLICATES")
		normalization = config.get("normalization", "method")
	else:
		normalization = "none"

	if config.getboolean("general","use_ratios"):
		rep_treatment = config.get("general","ratio_rep_treatment")
		print(">> CONVERTING TO MATRIX:\tRATIOS")
		subprocess.call([src_dir+'Conversion/MaxQToMatrix.R', '-i', keys_file, '-d', tmp_in_file, '-o', tmp_out_file, '-r', rep_treatment, '-m', 'Ratio_H_L', '-n', normalization])

	if config.getboolean("general","use_intensities"):
		print(">> CONVERTING TO MATRIX:\tINTENSITIES")
		rep_treatment = config.get("general","intensity_rep_treatment")
		subprocess.call([src_dir+'Conversion/MaxQToMatrix.R', '-i', keys_file, '-d', tmp_in_file, '-o', tmp_out_file, '-r', rep_treatment, '-m', 'Intensity', '-n', normalization])
	tmp_in_file = tmp_out_file

	## #############################
	## 3. perform limma differential 

	bname = bname+"_LIM"
	tmp_out_file = output_dir+bname+".txt"

	if config.getboolean("limma","enabled"):
		print(">> LIMMA DIFFRENTIAL CHANGES")
		design_file = dir+config.get("limma","design")
		if config.getboolean("limma","enable_contrasts"):
			contrast_file = dir+config.get("limma","contrast")
		else:
			contrast_file = "none"
		subprocess.call([src_dir+'Stats/Limma.R', '-d', tmp_in_file, '-o', tmp_out_file, '-m', design_file, '-c', contrast_file])
	tmp_in_file = tmp_out_file		

	## ######################
	## 4. flatten per protein

	bname = bname+"_PRT"
	tmp_out_file = output_dir+bname+".txt"	
	
	if config.getboolean("flatten","enabled"):
		print(">> FLATTEN TO PROTEINS")
		pvl_method = config.get("flatten","pvl_method")
		subprocess.call([src_dir+'Stats/Flatten.R', '-d', tmp_in_file, '-o', tmp_out_file, '-m', pvl_method])
	tmp_in_file = tmp_out_file		
		
	## ###########
	## 5. annotate

	bname = bname+"_RES"
	# tmp_out_file = output_dir+bname+".txt"
	tmp_out_file = dir+ config.get("files", "output_dir") + "/" + config.get("files", "output_file")

	if config.getboolean("annotate","enabled"):
		print(">> ANNOTATING WITH UNIPROT")
		uniprot_dir = config.get("annotate","file_dir")
		species = config.get("annotate","species")
		subprocess.call([src_dir+'Conversion/AnnotateWithUniprot.R', '-d', tmp_in_file, '-o', tmp_out_file, '-s', species, '-k', "uniprot_ac", '-u', uniprot_dir])
	tmp_in_file = tmp_out_file

class Usage(Exception):
	def __init__(self, msg):
		self.msg = msg		


def main(argv=None):

	output				= ""
	config 				= ""


	if argv is None:
		argv = sys.argv
	try:
		try:
			opts, args = getopt.getopt(argv[1:], "ho:crfa", ["help", "output=", "config=", "remove-carryover=", "format=", "alone="])
		except getopt.error, msg:
			raise Usage(msg)
	
		# option processing
		for option, value in opts:
			if option == "-v":
				verbose = True
			if option in ("-h", "--help"):
				raise Usage(help_message)
			if option in ("-c", "--config"):
				config = value
			if option in ("-o", "--output"):
				output = value

		# config = "/Users/everschueren/Projects/HPCKrogan/Scripts/MSPipeline/tests/ptm/ptm.cfg"
		# config = "/Users/everschueren/Projects/HPCKrogan/Scripts/MSPipeline/tests/apms_maxq/apms_maxq.cfg"
		## first read config from file 
		cfg = readConfig(config)
		## then overwrite with command line args
		## to implement

		## call pipeline with config
		pipeline(cfg)

	except Usage, err:
		print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
		print >> sys.stderr, "\t for help use --help"
		return 2


if __name__ == "__main__":
		sys.exit(main())
