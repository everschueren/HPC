#!/usr/bin/env python
# encoding: utf-8
"""
APMSPipeline.py

Created by erik verschueren on 2013-06-26.
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
	design_file = dir+config.get("files","design")	

	##derived
	bname = os.path.splitext(os.path.basename(data_file))[0]

	## ################################ 
	## 1. read maxquant data and filter
	
	oxidations = False
	unique = False
	contaminants = False
	modification = None

	if config.getboolean("filters","enabled"):
		bname = bname+"_FLT"
		tmp_out_file = dir+"/processed/"+bname+".txt"
		oxidations = config.get("filters","oxidations")
		unique = config.get("filters","unique")
		contaminants = config.get("filters","contaminants")
		modification = config.get("filters","modification")
		subprocess.call(['Conversion/MaxQFilter.R', '-d', data_file, '-o', tmp_out_file, '--contaminants_filter', contaminants, '-u', unique, '-x', oxidations, '-s', modification])

	## #############################################
	## 2. convert to matrix of normalized log-ratios
	
	if config.get("normalization", "enabled"):
		normalization = config.get("normalization", "method")
	else:
		normalization = "none"

	if config.get("general","use_ratios"):
		rep_treatment = config.get("general","ratio_rep_treatment")
		#Conversion/MaxQToMatrix.R -i keys_file -d data_file -o tmp_out_file -r rep_treatment -m Ratio_H_L -f $FILTER -u $UNIQUE_ONLY

	if config.get("general","use_intensities"):
		rep_treatment = config.get("general","intensity_rep_treatment")
		subprocess.call(['Conversion/MaxQToMatrix.R', '-i', keys_file, '-d', data_file, '-o', tmp_out_file, '-r', rep_treatment, '-m', 'Intensities', '-n', normalization])
		
	## limma

	## flatten per protein 



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

		config = "/Users/everschueren/Projects/HPCKrogan/Scripts/MSPipeline/tests/ptm/ptm.cfg"		
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
