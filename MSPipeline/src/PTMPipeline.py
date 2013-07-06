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

	##derived
	bname = os.path.splitext(os.path.basename(data_file))[0]



class Usage(Exception):
	def __init__(self, msg):
		self.msg = msg		


def main(argv=None):

	remove_carryover 	= True
	format 				= "simple"
	alone 				= True 
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
			if option in ("-r", "--remove-carryover"):
				remove_carryover = value
			if option in ("-f", "--format"):
				format = value
			if option in ("a", "--alone"):
				alone = value

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
