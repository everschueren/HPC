#!/usr/bin/env python
# encoding: utf-8
"""
SaintWrapper.py

Created by erik verschueren on 2013-07-15.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

import os
import sys
import getopt
import subprocess
import csv

help_message = '''
The help message goes here.
'''

def skip_comments(iterable):
    for line in iterable:
        if not line.startswith('#'):
            yield line

def collapseFileTodict(collapse_file):
	d = {}
	with open(collapse_file,"r") as cf:
		cf_reader = csv.reader(cf, delimiter="\t")
		for line in cf_reader:
			key = line[0]
			val = line[1]
			d[key] = val
	return d

def removeFileTodict(remove_file):
	d = {}
	with open(remove_file,"r") as rf:
		rf_reader = csv.reader(skip_comments(rf), delimiter="\t")
		for line in rf_reader:
			key = line[0]
			d[key] = True
	return d

def prospectorToSaintFormat(prospector_file, dir, collapse_file, remove_file, format):

	interactions_file = dir + "/" + "saint_interactions.txt"
	bait_file = dir + "/" + "saint_baits.txt"
	prey_file = dir + "/" + "saint_preys.txt"

	collapse_dict = collapseFileTodict(collapse_file)
	remove_dict = removeFileTodict(remove_file)
	
	if format =="S":
		ip_idx = 1
		bait_idx = 0
		prey_idx = 4
		tsc_idx = 5
		mw_idx = 9
	elif format == "F":
		ip_idx = 2
		bait_idx = 8
		prey_idx = 22
		tsc_idx = 23
		mw_idx = 27

	## write the interactions file
	with open(prospector_file,"r") as csvfile:
		reader = csv.reader(csvfile, delimiter="\t")
		reader.next()
		reader.next()
		with open(interactions_file,"w") as interactions_out:
			with open(bait_file,"w") as bait_out:
				with open(prey_file,"w") as prey_out:
			
					prev_ip = ""
					prey_set = set()
		
					for line in reader:
						ip = line[ip_idx]
						
						## check if the IP is in the 'remove'  list
						if ip not in remove_dict.keys():
							bait = line[bait_idx]
							if bait in collapse_dict.keys():
								bait = collapse_dict[bait] ## give the bait the name from the collapse list
							prey = line[prey_idx]
							tsc = line[tsc_idx]
							mw = line[mw_idx]
		
							## 1. write out ALL interactions
							interactions_out.write(ip+"\t"+bait+"\t"+prey+"\t"+tsc+"\n")
		
							## 2. for the bait file only write out info if we're dealing with a new IP
							## assuming all IPs are SORTED in the prospectorfile
							if ip != prev_ip:
								bait_l = bait.lower()	
								if(bait_l == "negative" or bait_l == "negative_strep" or bait_l == "negative_flag" or bait_l == "control" or bait_l == "none" or bait_l == "vector"):
									control = "C" ## CONTROL
								else:
									control = "T" ## TEST
								bait_out.write((ip+"\t"+bait+"\t"+control+"\n"))
								prev_ip = ip
		
							## 3. write out prey info if it wasn't observed before
							if prey not in prey_set:
								prey_out.write(prey+"\t"+mw+"\t"+prey+"\n")
								prey_set.add(prey)
			
			interactions_out.close()
			bait_out.close()
			prey_out.close()


class Usage(Exception):
	def __init__(self, msg):
		self.msg = msg		


def main(argv=None):

	dir = os.getcwd()

	## DEFAULTS ########
	prospector_file = ""
	collapse_file = ""
	remove_file= ""
	format = "S"

	if argv is None:
		argv = sys.argv
	try:
		try:
			opts, args = getopt.getopt(argv[1:], "ho:dpcrf", ["help", "output=", "dir=", "prospector_file=", "collapse_file=", "remove_file=", "data_format="])
		except getopt.error, msg:
			raise Usage(msg)
	
		# option processing
		for option, value in opts:
			if option == "-v":
				verbose = True
			if option in ("-h", "--help"):
				raise Usage(help_message)
			if option in ("-d", "--dir"):
				dir = value
			if option in ("-p", "--prospector_file"):
				prospector_file = value
			if option in ("-c", "--collapse_file"):
				collapse_file = value
			if option in ("-r", "--remove_file"):
				remove_file = value
			if option in ("-f", "--data_format"):
				format = value

		print(collapse_file)
		prospectorToSaintFormat(prospector_file, dir, collapse_file, remove_file, format)

	except Usage, err:
		print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
		print >> sys.stderr, "\t for help use --help"
		return 2


if __name__ == "__main__":
		sys.exit(main())
