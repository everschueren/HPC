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

# MS_PIPELINE_PATH = "/netapp/home/erikv/Scripts/MSPipeline/"
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
	print(">> KROGAN APMS PIPELINE")

	## read filesÏ
	dir = config.get("files","dir")	
	data_file = dir+config.get("files","data")	
	keys_file = dir+config.get("files","keys")
	remove_file = dir+config.get("files","remove")
	collapse_file = dir+config.get("files","collapse")
	exclusions_file = dir+config.get("files","exclusions")
	output_dir = dir + config.get("files","output_dir") + '/'
	
	##derived
	bname = os.path.splitext(os.path.basename(data_file))[0]
	
	#################################
	## PREPARE CONTAINERS FOR RESULTS

	if not os.path.exists(output_dir):
		os.mkdir(output_dir)
	mist_metrics_file=""
	mist_self_score_file=""
	mist_hiv_score_file=""
	comppass_score_file=""
	saint_score_file=""
	ip_out_file=""


	## ################## 
	## filter contaminants
	tmp_in_file = data_file
	
	if config.getboolean('general','filter_contaminants'): ## simple format so merge keys
		print ">> FILTERING COMMON CONTAMINANTS\t\t"
		tmp_out_file = output_dir+bname+"_filtered.txt"
		subprocess.call([src_dir+'FileIO/filterContaminants.R', '-k', keys_file, '-d', tmp_in_file, '-o', tmp_out_file, '-c', config.get("general","filter_contaminants")])

	else:
		tmp_out_file = tmp_in_file

	## ################## 
	## merge additions

	tmp_in_file = tmp_out_file
	
	if config.get('general','data_format') == 'S': ## simple format so merge keys
		print(">> MERGING KEYS TO DATA\t\t" + tmp_in_file)
		bname = bname + "_wKEYS"
		tmp_out_file = output_dir+bname+".txt"
		f = open(tmp_out_file,'w')
		subprocess.call([src_dir+'FileIO/merge_additions.pl', keys_file, tmp_in_file, config.get("general","data_format")], stdout=f)
		f.close()
		
	else:
		print(">> FULL FORMAT: NO KEYS WERE MERGED\t\t" + tmp_in_file)
		tmp_out_file = tmp_in_file ## no action

	## ##########################
	## merge with mastertable

	tmp_in_file = tmp_out_file

	if config.getboolean('general','master_include'):
		
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
				
				lines = d.readlines()[1:]
				f.writelines(lines)
			with open(master_file, 'r') as m:
				
				lines = m.readlines()[2:]
				f.writelines(lines)	
			f.close()	

	
	## ################### 
	## remove carryover

	tmp_in_file = tmp_out_file

	if config.getboolean("general","remove_carryover"):
		
		print(">> REMOVING CARRYOVER FROM\t\t" + tmp_in_file)
		bname = bname+'_NoC'
		tmp_out_file = output_dir+bname+".txt"
		f = open(tmp_out_file,'w')
		subprocess.call([src_dir+'FileIO/carryover_removal_comprehensive.pl', tmp_in_file, config.get("general","data_format")], stdout=f)
		f.close()
		## SET THESE UP TO BE (OPTIONALLY) USED BY SAINT LATER
		prospector_noc = tmp_out_file
	else: 
		prospector_noc = tmp_in_file

	## ####################### 
	## converting to matrix

	tmp_in_file = tmp_out_file

	print(">> CONVERTING TO MATRIX\t\t" + tmp_in_file)
	tmp_name = bname
	bname = bname + "_MAT"
	tmp_out_file = output_dir+bname+".txt"
	f = open(tmp_out_file,'w')
	subprocess.call([src_dir+'FileIO/convert_to_matrix_format_comprehensive.pl', tmp_in_file, config.get("general","data_format"), remove_file, collapse_file, exclusions_file], stdout=f)
	f.close()
	tmp_in_file = tmp_out_file
	
	print(">> CHECKING MATRIX INTEGRITY\t\t" + tmp_in_file)
	subprocess.call([src_dir+'FileIO/FileCheck.py', tmp_in_file])
	
	if(config.getboolean("general","ip_info")):
		print(">> COLLECTING PREY-IP INFORMATION\t\t" + tmp_in_file)
		tmp_out_file = output_dir + tmp_name + "_IPs" + ".txt"
		subprocess.call([src_dir+'Stats/Report.R', '-d', tmp_in_file, '-o', tmp_out_file])
		ip_out_file = tmp_out_file

	####################
	## CLUSTERING IPS

	if(config.getboolean("cluster","enabled")):
		print(">> CLUSTER IP MATRIX\t\t" + tmp_in_file)
		tmp_out_file = output_dir+bname+'_CLUSTERED.pdf'
		font_scale = config.get("cluster", "font_scale")
		subprocess.call([src_dir+'Stats/Cluster.R', '-d', tmp_in_file, '-o', tmp_out_file, '-s', font_scale])


	###################################
	## MiST scoring with HIV/SELF params

	mist_dir = output_dir+"/MIST/"
	if not os.path.exists(mist_dir):
		os.mkdir(mist_dir)

	if config.getboolean('mist','enabled'):
		print(">> SCORING MiST [HIV PARAMS]\t\t" + tmp_in_file)
		tmp_out_file = mist_dir+bname+'_MIST_HIV'+".txt"
		subprocess.call([src_dir+'MiST/MiST.py', tmp_in_file, tmp_out_file, config.get('mist','filter'), config.get('mist','training')]) 
	
	mist_hiv_score_file = mist_dir+bname+'_MIST_HIV'+"_scores.txt"
	
	if not os.path.exists(mist_hiv_score_file):
		mist_hiv_score_file = ""

	if config.getboolean('mist_self','enabled'):
		print(">> SCORING MiST [SELF PARAMS]\t\t" + tmp_in_file)
		tmp_out_file = mist_dir+bname+'_MIST_SELF'+".txt"
		subprocess.call([src_dir+'MiST/MiST.py', tmp_in_file, tmp_out_file, config.get('mist_self','filter'), config.get('mist_self','training')])

	mist_self_score_file = mist_dir+bname+'_MIST_SELF'+"_scores.txt"
	mist_metrics_file = mist_dir+bname+'_MIST_SELF'+"_metrics.txt"

	if not os.path.exists(mist_self_score_file):
		mist_self_score_file = ""
		mist_metrics_file = ""		
		
	#######################
	## COMPPASS scoring 

	comppass_dir = output_dir+"/COMPPASS/"
	if not os.path.exists(comppass_dir):
		os.mkdir(comppass_dir)

	if config.getboolean('comppass','enabled'):
		print(">> SCORING COMPPASS\t\t" + tmp_in_file)
		tmp_out_file = comppass_dir+bname+'_COMPPASS'+".txt"
		if config.getboolean('comppass','resampling'):
			resampling='TRUE'
		else:
			resampling='FALSE'
		subprocess.call([src_dir+'Comppass/Comppass.R', '-d', tmp_in_file, '-o', tmp_out_file, '-r', resampling])

	comppass_score_file = comppass_dir+bname+'_COMPPASS'+".txt"
		
	if not os.path.exists(comppass_score_file):
		comppass_score_file = ""

	###################
	## SAiNT scoring

	saint_dir = output_dir+"SAINT/"
	if not os.path.exists(saint_dir):
		os.mkdir(saint_dir)

	if config.getboolean('saint','enabled'):
		os.chdir(saint_dir)	
		print(">> SCORING SAINT\t\t" + prospector_noc)
		control_reduction = config.get("saint", "control_reduction")
		burnin = config.get("saint", "burnin")
		iter = config.get("saint", "iter")
		lowmode = config.get("saint", "lowmode")
		minfold = config.get("saint", "minfold")
		normalize = config.get("saint", "normalize")
		
		subprocess.call([src_dir+"Saint/SaintInput.py", '--dir', saint_dir, '--prospector_file', prospector_noc, '--collapse_file', collapse_file, '--remove_file', remove_file, '--data_format', config.get('general','data_format')])
		subprocess.call([bin_dir+"saint/saint-reformat", saint_dir+"saint_interactions.txt", saint_dir+"saint_preys.txt", saint_dir+"saint_baits.txt", control_reduction])
		subprocess.call([bin_dir+'saint/saint-spc-ctrl', saint_dir+"interaction.new", saint_dir+"prey.new", saint_dir+"bait.new", burnin, iter, lowmode, minfold, normalize])
		# for cases with NO controls:
		#subprocess.call([bin_dir+'saint/saint-spc-noctrl', saint_dir+"interaction.new", saint_dir+"prey.new", saint_dir+"bait.new", burnin, iter, "0.2", "0.1", "0", normalize])
	
	saint_score_file = saint_dir + "/RESULT/unique_interactions" 		

	if not os.path.exists(saint_score_file):
		saint_score_file = ""

	####################
	## COLLECT SCORES

	if config.getboolean('collect','enabled'):
		print(">> COLLECT SCORES AND ANNOTATE\t\t")
		tmp_out_file = output_dir+bname+'_ALLSCORES'+".txt"

		uniprot_dir = config.get("collect", "uniprot_dir")
		species = config.get("collect", "species")
		annotate =  config.get("collect", "annotate")

		# print(src_dir+'FileIO/CompileResults.R'+  ' -o '+ tmp_out_file+ ' -f '+ mist_hiv_score_file+ ' -t '+ mist_self_score_file+ ' -m '+ mist_metrics_file+ ' -c '+ comppass_score_file+ ' -s '+ saint_score_file+ ' -u '+ uniprot_dir+ ' -n '+ species+ ' -i '+ src_dir + '-a ' + annotate)
		subprocess.call([src_dir+'FileIO/CompileResults.R', '-d', "", '-o', tmp_out_file, '-f', mist_hiv_score_file, '-t', mist_self_score_file, '-m', mist_metrics_file, '-c', comppass_score_file, '-s', saint_score_file, '-p', ip_out_file,'-u',uniprot_dir, '-n', species, '-i', src_dir, '-a', annotate ])


	####################
	## QUALITY CONTROL PLOTS

	if config.getboolean('qc','enabled'):
		print(">> PLOTTING FOR QUALITY CONTROL\t\t")
		
		tmp_out_file = dir+"../summary"
		tmp_data_file = output_dir + os.path.splitext(os.path.basename(data_file))[0] + "_wKEYS_NoC.txt"

		# print(src_dir+'FileIO/CompileResults.R'+  ' -o '+ tmp_out_file+ ' -f '+ mist_hiv_score_file+ ' -t '+ mist_self_score_file+ ' -m '+ mist_metrics_file+ ' -c '+ comppass_score_file+ ' -s '+ saint_score_file+ ' -u '+ uniprot_dir+ ' -n '+ species+ ' -i '+ src_dir + '-a ' + annotate)
		subprocess.call([src_dir+'Stats/qualityCheck.R', '-d', tmp_data_file, '-o', tmp_out_file])


	####################
	## ENRICHMENT ANALYSIS

	if config.getboolean('enrich','enabled'):
		print(">> PERFORMING ENRICHMENT ANALYSIS\t\t")
		
		enriched_dir = output_dir+"enriched/"
		print enriched_dir
		if not os.path.exists(enriched_dir):
		    os.makedirs(enriched_dir)

		all_scores = output_dir+bname+'_ALLSCORES'+".txt"
		
		# if we want to only enrich the scores above a threshold: 
		if config.getboolean('enrich','threshold'):
			print(">> PERFORMING THRESHHOLDING ON SCORES\t\t")
			comppass_thresh = config.get('enrich','comppass_thresh')
			mist_thresh = config.get('enrich','mist_thresh')
			thresh_file = output_dir+bname+'_ALLSCORES_'+'c'+comppass_thresh.replace(".","")+'m'+mist_thresh.replace(".","")+".txt"
			subprocess.call([src_dir+'FileIO/ThreshholdResults.R', '-d', all_scores, '-o', thresh_file, '-c', comppass_thresh, '-m', mist_thresh, '-s', config.get('enrich','species')])
			all_scores = thresh_file
		#
		
		bait_col = config.get("enrich", "bait_col")
		prey_col = config.get("enrich", "prey_col")
		
		subprocess.call([src_dir+'Meta/Enrichment.R', '-d', all_scores, '-o', enriched_dir, '-s', bait_col, '-p', prey_col])


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
			opts, args = getopt.getopt(argv[1:], "hoc:rfa", ["help", "output=", "config=", "remove-carryover=", "format=", "alone="])
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
		## config = "/Users/everschueren/Projects/HPCKrogan/Scripts/MSPipeline/tests/hcv/TestConfig.txt"
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
