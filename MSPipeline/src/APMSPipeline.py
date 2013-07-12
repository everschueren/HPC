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

# mist_self_metrics_file = "HCV-293T-Andy-results_wMT_MIST_SELF_metrics.txt"
# mist_self_score_file = "HCV-293T-Andy-results_wMT_MIST_SELF_scores.txt"
# mist_hiv_score_file = "HCV-293T-Andy-results_wMT_MIST_HIV_scores.txt"
# comppass_score_file = "HCV-293T-Andy-results_wMT_COMPPASS.txt"


help_message = '''
The help message goes here.
'''

src_dir = "/Users/everschueren/Projects/HPCKrogan/Scripts/MSPipeline/src/"
bin_dir = "/Users/everschueren/Projects/HPCKrogan/Scripts/MSPipeline/bin/"


def readConfig(filename):
	file = open(filename)
	config = ConfigParser.RawConfigParser()
	config.readfp(file)
	
	return config

def pipeline(config):
	print(">> KROGAN APMS PIPELINE")

	## read filesÏ
	dir = config.get("files","dir")	
	data_file = dir+config.get("files","data")	
	keys_file = dir+config.get("files","keys")
	remove_file = dir+config.get("files","remove")
	collapse_file = dir+config.get("files","collapse")
	exclusions_file = dir+config.get("files","exclusions")
	output_dir = config.get("files","output_dir") + '/'
	
	##derived
	bname = os.path.splitext(os.path.basename(data_file))[0]

	## 

	mist_metrics_file=""
	mist_self_score_file=""
	mist_hiv_score_file=""
	comppass_score_file=""

	## ################## 
	## 1. merge additions

	tmp_in_file = data_file
	print(">> MERGING KEYS TO DATA\t\t" + tmp_in_file)
	tmp_out_file = output_dir+bname+".txt"
	f = open(tmp_out_file,'w')
	subprocess.call([src_dir+'FileIO/merge_additions.pl', keys_file, tmp_in_file, config.get("general","data_format")], stdout=f)
	f.close()

	## ##########################
	## 1.B merge with mastertable

	if config.getboolean('general','master_include'):
		tmp_in_file = tmp_out_file
		master_file = config.get('general','master_file')
		
		print(">> APPENDING DATA FILE TO MASTERFILE\t\t"+master_file)

		if(config.get('general','master_format')!=config.get('general','data_format')):
			print(">> DATA FORMAT != MASTER FORMAT")
			print(">> EXITING ...")
			exit(1)
		else:
			bname = bname+'_wMT'
			tmp_out_file = output_dir+bname+".txt"
			f = open(tmp_out_file,'w')
			with open(tmp_in_file, 'r') as d:
				lines = d.readlines()[3:]
				f.writelines(lines)
			with open(master_file, 'r') as m:
				lines = m.readlines()[2:]
				f.writelines(lines)	
			f.close()	

	

	## ################### 
	## 2. remove carryover

	if config.get("general","remove_carryover"):
		tmp_in_file = tmp_out_file
		print(">> REMOVING CARRYOVER FROM\t\t" + tmp_in_file)
		tmp_out_file = output_dir+bname+'_NoC'+".txt"
		f = open(tmp_out_file,'w')
		subprocess.call([src_dir+'FileIO/carryover_removal_comprehensive.pl', tmp_in_file, config.get("general","data_format")], stdout=f)
		f.close()

	## ####################### 
	## 3. converting to matrix
	# $data = $ARGV[0];
	# $format = $ARGV[1];
	# $inputdir = $ARGV[2];
	# $remove = $ARGV[3];
	# $collapse = $ARGV[4];
	# $exclusions = $ARGV[5];

	tmp_in_file = tmp_out_file
	print(">> CONVERTING TO MATRIX\t\t" + tmp_in_file)
	tmp_out_file = output_dir+bname+'_MAT'+".txt"
	f = open(tmp_out_file,'w')
	subprocess.call([src_dir+'FileIO/convert_to_matrix_format_comprehensive.pl', tmp_in_file, config.get("general","data_format"), remove_file, collapse_file, exclusions_file], stdout=f)
	f.close()

	## ########################## 
	## 4. checking file integrity

	tmp_in_file = tmp_out_file
	print(">> CHECKING MATRIX INTEGRITY\t\t" + tmp_in_file)
	subprocess.call([src_dir+'FileIO/FileCheck.py', tmp_in_file])

	##################################
	## 5. MiST scoring with HIV params

	if config.getboolean('mist','enabled'):
		print(">> SCORING MiST [HIV PARAMS]\t\t" + tmp_in_file)
		tmp_out_file = output_dir+bname+'_MIST_HIV'+".txt"
		subprocess.call([src_dir+'MiST/MiST.py', tmp_in_file, tmp_out_file, config.get('mist','filter'), config.get('mist','training')]) 
		mist_hiv_score_file = output_dir+bname+'_MIST_HIV'+"_scores.txt"

	###################################
	## 6. MiST scoring with SELF params

	if config.getboolean('mist_self','enabled'):
		print(">> SCORING MiST [SELF PARAMS]\t\t" + tmp_in_file)
		tmp_out_file = output_dir+bname+'_MIST_SELF'+".txt"
		subprocess.call([src_dir+'MiST/MiST.py', tmp_in_file, tmp_out_file, config.get('mist_self','filter'), config.get('mist_self','training')]) 
		mist_self_score_file = output_dir+bname+'_MIST_SELF'+"_scores.txt"
		mist_metrics_file = output_dir+bname+'_MIST_SELF'+"_metrics.txt"

	######################
	## 7. COMPPASS scoring 

	## replaced by newer version
	# if config.getboolean('comppass','enabled'):
	# 	print(">> SCORING COMPPASS\t\t" + tmp_in_file)
	# 	tmp_out_file = output_dir+bname+'_COMPPASS'+".txt"
	# 	f = open(tmp_out_file,'w')
	# 	subprocess.call(['Comppass/GetComppass.py', tmp_in_file], stdout=f)
	# 	f.close()

	if config.getboolean('comppass','enabled'):
		print(">> SCORING COMPPASS\t\t" + tmp_in_file)
		tmp_out_file = output_dir+bname+'_COMPPASS'+".txt"
		subprocess.call([src_dir+'Comppass/Comppass.R', '-d', tmp_in_file, '-o', tmp_out_file])
		comppass_score_file = tmp_out_file

	###################
	## 7. SAiNT scoring

	if config.getboolean('saint','enabled'):
		print(">> SCORING SAINT\t\t" + tmp_in_file)
		tmp_out_file = output_dir+bname+'_SAINT'+".txt"
		subprocess.call([bin_dir+'saint/saint-spc-noctrl-matrix', tmp_in_file, tmp_out_file, config.get('saint','nburnin'), config.get('saint','niter'), config.get('saint','ff')])

	####################
	## 8. COLLECT SCORES

	print(">> COLLECT SCORES AND ANNOTATE\t\t")
	tmp_out_file = output_dir+bname+'_ALLSCORES'+".txt"

	## replaced by newer version
	# f = open(tmp_out_file,'w')
	# file_dir = config.get("collect","file_dir")
	# subprocess.call([src_dir+'FileIO/convert_out_to_final_comprehensive.pl', output_dir+bname, config.get("general","data_format"), collapse_file, file_dir + config.get("collect","gene_names"), file_dir + config.get("collect","uniprot_entrez")], stdout=f )
	# f.close()

	uniprot_dir = config.get("collect", "uniprot_dir")
	species = config.get("collect", "species") 
	subprocess.call([src_dir+'FileIO/CompileResults.R', '-d', "", '-o', tmp_out_file, '-f', mist_hiv_score_file, '-t', mist_self_score_file, '-m', mist_metrics_file, '-c', comppass_score_file, '-u', uniprot_dir, '-s', species, '-i', src_dir ])


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
