#!/bin/bash

# Command structure:
#                                           		verify-input.py
#  [0] ./easymap.sh										.
#  [1] $project_name									.
#  [2] $workflow[ins/snp]								.                      Maybe add a 3rd workflow: Analysis of SNPs
#  [3] $data_source[exp/sim]							.
#  [4] $ref_seq											*
#  [5] $ins_seq											*
#  [6] $gff_file										*
#  [7] $ann_file										*
#  [8] $read_s											*
#  [9] $read_f											*
# [10] $read_r											*
# [11] $lib_type_sample[se/pe]							.
# [12] $read_s_ctrl										*
# [13] $read_f_ctrl										*
# [14] $read_r_ctrl										*
# [15] $lib_type_ctrl [se/pe]							.
# [16] $is_ref_strain [ref/noref]						.                      Only for linkage analysis mapping
# [17] $cross_type [oc/bc]								.                      Only for linkage analysis mapping
# [18] $snp_analysis_type [par/f2wt]					.
# [19] $control_parental [mut/nomut] 					.                      Only for linkage analysis mapping
# [20] $sim_mut											.                      nbr+mod
# [21] $sim_recsel										.                      rfd+pos+mod+nre
# [22] $sim_seq											.                      rd+rl+fl+ber+gbs
# [23] $stringency

# sim-mut.py
# nbr:		${20}[0]
# mod:		${20}[1]
# con:		$4
# out:		constant
#
# sim-recsel.py
# rfd:		${21}[0]
# pos:		${21}[1]
# mod:		${21}[2]
# nre:		${21}[3]
# mut:		constant
# pol:		constant
# out:		constant
#
# sim-seq.py
# if:		constant
# out:		constant
# mod:		${11} or ${15}
# rd:		${22}[0]
# rl:		${22}[1]
# fl:		${22}[2]
# ber:		${22}[3]
# gbs:		${22}[4]
#
# The command has three levels of arguments:
# 1st: Main arguments are separated by whitespace.
# 2nd: Each simulator script receives an argument composed of a string of arguments reparated
#      by the '+' character.
# 3rd: Some of the 2nd level arguments are in turn composed of a string of values separated by
#      the ',' character.

############################################################
# Activate Python virtual environment
source src/virtualenv-15.1.0/easymap-env/bin/activate

############################################################
# Obtain and store date and time in format with no spaces
timestamp=$(date "+%F-%T")

############################################################
# Get command arguments and assign them to variables

project_name=user_projects/$timestamp"_"$1
workflow=$2
data_source=$3
ref_seq=$4
ins_seq=$5
gff_file=$6
ann_file=$7
read_s=$8
read_f=$9
read_r=${10}
lib_type_sample=${11}
read_s_ctrl=${12}
read_f_ctrl=${13}
read_r_ctrl=${14}
lib_type_ctrl=${15}
is_ref_strain=${16}
cross_type=${17}
snp_analysis_type=${18}
control_parental=${19}
sim_mut=${20}
sim_recsel=${21}
sim_seq=${22}
stringency=${23}

############################################################
# Several necessary checking/preparation steps before actually running easymap

# If easymap is being executed from server, set path for Python unzipped eggs in easymap/tmp dir
#if [ $env == 'server' ]; then
#	export PYTHON_EGG_CACHE=./tmp
#fi

# Declare a flag variable that will be used as exit code, and set it to 0 (no error)
exit_code=0

# Create 'user_projects' folder if the user has deleted it by mistake
[ -d user_projects ] || mkdir user_projects

# Store the location of each folder in variables
f0=user_data
f1=1_intermediate_files
f2=2_logs
f3=3_workflow_output

# Create project folder and subfolders
mkdir $project_name
mkdir $project_name/$f1
mkdir $project_name/$f2
mkdir $project_name/$f3

# Change permisssion so www-data can read and write in all folders of the project
#chmod -R 777 $project_name

# Deprecated
# Delete intermediate and final files from previous executions, For that, check whether dirs have any
# have content (folders or files) and, if so, remove it
[ "$(ls -A $project_name/$f1)" ] && rm --recursive $project_name/$f1/*
[ "$(ls -A $project_name/$f2)" ] && rm --recursive $project_name/$f2/*
[ "$(ls -A $project_name/$f3)" ] && rm --recursive $project_name/$f3/*

# Define path of log file and create it
my_log_file=$project_name/$f2/log.log
touch $my_log_file

# Define path of status file and create it
my_status_file=$project_name/$f2/status
touch $my_status_file
chmod 666 $my_status_file
echo 'status:running' >> $my_status_file
echo 'pid easymap '$$ >> $my_status_file


# Check that the folders /user_data and /user_data/gnm_ref exist. If they do not, 
if ! [ -d $f0 ]; then
	echo $(date "+%F > %T")": Execution could not start because folder user_data could not be found. Please, create the folder and use it to place the files you want to analyze." > $my_log_file
	echo 'status:error' >> $my_status_file
	exit_code=1
	#echo $exit_code
	exit
fi

############################################################
# Start easymap

echo $(date "+%F > %T")": Execution of project {" $project_name "} started." > $my_log_file
echo "" >> $my_log_file

echo "Program:										" $0 >> $my_log_file
echo "Project name:									" $timestamp"_"$1 >> $my_log_file
echo "Workflow:										" $2 >> $my_log_file
echo "Data source:									" $3 >> $my_log_file
echo "Reference sequence:							" $4 >> $my_log_file
echo "Insertion sequence:							" $5 >> $my_log_file
echo "GFF file:										" $6 >> $my_log_file
echo "Annotation file:								" $7 >> $my_log_file
echo "Single-end reads (test sample):			" $8 >> $my_log_file
echo "Forward reads (test sample):				" $9 >> $my_log_file
echo "Reverse reads (test sample):				" ${10} >> $my_log_file
echo "Library type (test sample):				" ${11} >> $my_log_file
echo "Single-end reads (control sample):			" ${12} >> $my_log_file
echo "Forward reads (control sample):				" ${13} >> $my_log_file
echo "Reverse reads (control sample):				" ${14} >> $my_log_file
echo "Library type (control sample):				" ${15} >> $my_log_file
echo "Mutant strain [ref/noref]:					" ${16} >> $my_log_file
echo "Type of cross [bc/oc]:						" ${17} >> $my_log_file
echo "SNP analysis type [par/f2wt]:					" ${18} >> $my_log_file
echo "Parental used as control [mut/nomut/np]:		" ${19} >> $my_log_file
echo "Simulator (sim-mut.py) command:				" ${20} >> $my_log_file
echo "Simulator (sim-recsel.py) command:			" ${21} >> $my_log_file
echo "Simulator (sim-seq.py) command:				" ${22} >> $my_log_file
echo "Stringency:									" ${23} >> $my_log_file

echo "" >> $my_log_file
echo "######################################################" >> $my_log_file
echo $(date "+%F > %T")": Project data directories created." >> $my_log_file


############################################################
# Overwrite read_s, read_f and read_r
# 'test': F2 recessive phenotype (mutant phenotype in recessive mutations)
# 'control': if snp_analysis_type=par_mut/par_nomut, reads from one of the parentals used in cross; if snp_analysis_type=f2wt, F2 dominant phenotype (wildtype phenotype in recessive mutations)

if [ $data_source == 'sim' ]; then
	if [ $lib_type_sample == 'se' ]; then
		read_s=$project_name/$f1/sim_data/sim_seq_output/sample/se_reads.fq
		read_s_ctrl=$project_name/$f1/sim_data/sim_seq_output/control/se_reads.fq
	else
		read_f=$project_name/$f1/sim_data/sim_seq_output/sample/pe-for_reads.fq
		read_r=$project_name/$f1/sim_data/sim_seq_output/sample/pe-rev_reads.fq
		read_f_ctrl=$project_name/$f1/sim_data/sim_seq_output/control/pe-for_reads.fq
		read_r_ctrl=$project_name/$f1/sim_data/sim_seq_output/control/pe-rev_reads.fq
	fi
fi

if [ $data_source == 'exp' ]; then
	if [ $lib_type_sample == 'se' ]; then
		read_s=$f0/$8
	else
		read_f=$f0/$9
		read_r=$f0/${10}
	fi

	if [ $lib_type_ctrl == 'se' ]; then
		read_s_ctrl=$f0/${12}
	else
		read_f_ctrl=$f0/${13}
		read_r_ctrl=$f0/${14}
	fi
fi


############################################################
# Run 'process-input.sh'

echo $(date "+%F > %T")": STARTING INPUT PROCESSING..." >> $my_log_file

process_input=`./process_input/process-input.sh $my_log_file $project_name $workflow $data_source $lib_type_sample $ins_seq $read_s $read_f $read_r $gff_file $ann_file $read_s_ctrl $read_f_ctrl $read_r_ctrl $ref_seq $lib_type_ctrl`

if [ $process_input == 0 ]; then
	echo $(date "+%F > %T")": All inputs correct." >> $my_log_file
else 
	echo $(date "+%F > %T")": One or more inputs incorrect (see details above in this log). Quit." >> $my_log_file
	echo 'status:error' >> $my_status_file
	echo "Easymap analysis failed. See log file for more info"
	exit
fi


############################################################
# Run 'simulator.sh'

if [ $data_source == 'sim' ]; then
	echo $(date "+%F > %T")": STARTING DATA SIMULATION..." >> $my_log_file
	simulator=`./simulator/simulator.sh $my_log_file $project_name $workflow $lib_type_sample $ins_seq $sim_mut $sim_recsel $sim_seq $cross_type $is_ref_strain $control_parental $snp_analysis_type`
	
	if [ $simulator == 0 ]; then
		echo $(date "+%F > %T")": Simulation completed." >> $my_log_file
	else 
		echo $(date "+%F > %T")": Simulation failed (see details above in this log). Quit." >> $my_log_file
		echo 'status:error' >> $my_status_file
		echo "Easymap analysis failed. See log file for more info"
		exit
	fi
fi

############################################################
# Run the chosen analysis workflow

if [ $workflow == 'ins' ]; then
	workflow_result=`./workflows/workflow-ins.sh $my_log_file $project_name $workflow $data_source $lib_type_sample $ins_seq $read_s $read_f $read_r $gff_file $ann_file`

	if [ $workflow_result == 0 ]; then
		echo $(date "+%F > %T")": Analysis workflow finished correctly." >> $my_log_file
	else 
		echo $(date "+%F > %T")": Analysis workflow failed (see details above in this log)." >> $my_log_file
		echo "Easymap analysis failed. See log file for more info"
		echo 'status:error' >> $my_status_file
		exit
	fi
fi

if [ $workflow == 'snp' ]; then
	workflow_result=`./workflows/workflow-snp.sh $my_log_file $project_name $workflow $data_source $lib_type_sample $ins_seq $read_s $read_f $read_r $gff_file $ann_file $read_s_ctrl $read_f_ctrl $read_r_ctrl $cross_type $is_ref_strain $control_parental $snp_analysis_type $lib_type_ctrl $stringency`

	if [ $workflow_result == 0 ]; then
		echo $(date "+%F > %T")": Analysis workflow finished correctly." >> $my_log_file
	else 
		echo $(date "+%F > %T")": Analysis workflow failed (see details above in this log)." >> $my_log_file
		echo "Easymap analysis failed. See log file for more info"
		echo 'status:error' >> $my_status_file
		exit
	fi
fi

echo $(date "+%F > %T")": Execution of project {" $project_name "} finished." >> $my_log_file
echo 'status:finished' >> $my_status_file

# This message must remain as is because install.sh relies on it to know whether the installation was successfull or not.
echo "Easymap analysis properly completed."
