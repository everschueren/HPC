#scores datasets not in mastertable alone or against the mastertable

#ARG1: file name
#ARG2: against mastertable (M) or alone (A)
#ARG3: format: with experimental details (Full) or just prospector output (Simple)
#ARG4: carryover remove (1) or not (0) 
# e.g. bash scoring_script_comprehensive.sh VIF_Ctrl_VP.txt A F 0

rundir=pwd
prospectorfile=$1
pipeline=~/Projects/HPCKrogan/Scripts/Pipeline
scorebg=$2 #1 to score against mastertable or alone
format=$3 #F/S

###Score against mastertable with full format

if [ "$scorebg" = "M" ]
then
	if [ "$format" = "F" ]
	then

	#Merge Mastertable with individual dataset
    echo "\tMerging Mastertable with individual set"
	cat 2013_Apr24_MasterTable.txt "$prospectorfile"  > "$prospectorfile"_merged.txt
	prospectorfile="$prospectorfile"_merged.txt

	elif [ "$format" = "S" ] 
	then
	#Create the latest mastertable in the right format
	echo ">> Create the latest mastertable in the right format"
	awk -F "\t" '{print $3"\t"$9"\t"$21"\t"$22"\t"$23"\t"$24"\t"$25"\t"$26"\t"$27"\t"$28"\t"$29"\t"$30}' 2013_Apr24_MasterTable.txt.nocarryover > Mastertable.txt	
	cat Mastertable.txt "$prospectorfile"  > "$prospectorfile"_merged.txt
	prospectorfile="$prospectorfile"_merged.txt

	fi
fi

###############
######CARRYOVER
###############

carryover=$4 #1 to remove carryover
if [ $carryover -eq "1" ]
then
	echo ">> Removing carryover: " $prospectorfile
	perl $pipeline/carryover_removal_comprehensive.pl "$prospectorfile" "$format" > "$prospectorfile".nocarryover
	prospectorfile="$prospectorfile".nocarryover
fi
 
# 
# #######################
# ######CONVERT TO MATRIX
# #######################

echo ">> Converting to matrix: " $prospectorfile
perl $pipeline/convert_to_matrix_format_comprehensive.pl "$prospectorfile" "$format" > "$prospectorfile"_input.txt
echo ">> File Check: " $prospectorfile
python $pipeline/FileCheck.py "$prospectorfile"_input.txt


# ###########
# ######SCORE
# ###########

######MiST no training/using HIV metrics

echo ">> MiST with HIV metrics: " $prospectorfile"_input.txt"
python $pipeline/MiST/MiST.py "$prospectorfile"_input.txt "$prospectorfile"_notraining 0 0 
    
# ######MiST with training

echo ">> MiST with training: " $prospectorfile"_input.txt" 
python $pipeline/MiST/MiST_invert.py "$prospectorfile"_input.txt "$prospectorfile"_training 0 1

# ######CompASS

echo ">> Running ComppASS"
python $pipeline/GetComppass.py "$prospectorfile"_input.txt > "$prospectorfile"_compass
# echo -e "Done with compass"

######SAINT

echo ">> Running SAiNT" 
$pipeline/SAINT_v2.3.4/bin/saint-spc-noctrl-matrix "$prospectorfile"_input.txt "$prospectorfile"_out_saint 1500 10000 0.9    

# ############
# ######OUTPUT
# ############
 
echo ">> Convert all the scores to a final table with Entrez Gene Symbols and IDs"
perl $pipeline/convert_out_to_final_comprehensive.pl "$prospectorfile" "$format" | awk -F "\t" '$1' > "$prospectorfile"_scoring.txt
echo ">> Scores can be found in ""$prospectorfile"_scoring.txt
tr -d "\015" < "$prospectorfile"_scoring.txt > x ; mv x "$prospectorfile"_scoring.txt
 
