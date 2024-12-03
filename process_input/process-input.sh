#!/bin/bash


# Check all the files and report to log if they are correct or not, but always continue
# ultil all input files checked. This way, the user can know in a single run all the files
# that are not correct.
# At the same time, if any file is incorrect, set exit_code=1. master.sh retrieves this value and
# stops the workflow. Then, user can see in log stream which file was incorrect.
#
# TO DO: Maybe add error checks after the execution of each python program, as in the rest of the
# workflows.
#
#


#######################################################################################################
# This block of code does some steps previous to the analysis
# 
# Description:
# 
# 
# 
#######################################################################################################

#######################################################################################################
# TO DO:
# 
# If mode is SNP and source exp, control fastQ files (parental or F2 control) have to be checked as well.
# However, wait until SE/PE combinations issue is sorted out.
# 
# 
#######################################################################################################


# Set 'exit_code' (flag variable) to it's initial value (0)
exit_code=0

# This is the command sent by 'master.sh':
# ./process-input.sh $my_log_file $project_name $analysis_type $data_source $lib_type_sample $ins_seq $read_s $read_f $read_r $gff_file $ann_file`
# example: ./process-input.sh project/log.log project ins sim se ins.fa fq_se.fq fq_1.fq fq_2.fq gff.gff ann.ann

# Store the location of each folder in a variable
f0=user_data
f1=1_intermediate_files
f3=3_workflow_output

# Get command arguments and assign them to variables
# The reference to the template genome is not passed from 'master.sh' because it is static
my_log_file=$1
project_name=$2
analysis_type=$3
data_source=$4
lib_type_sample=$5
ins_seq=$f0/$6
read_s=$7
read_f=$8
read_r=$9
gff_file=$f0/${10}
ann_file=$f0/${11}
ann_option=${11}
read_s_ctrl=${12}
read_f_ctrl=${13}
read_r_ctrl=${14}
ref_seq=${15}
lib_type_ctrl=${16}

# Establish locations of reference genome
ref_seqs_dir=$f0/gnm_ref
ref_seqs_merged_dir=$project_name/$f1/gnm_ref_merged
ref_seqs_merged_file=$project_name/$f1/gnm_ref_merged/genome.fa

#######################################################################################################
# This block of code does the initial processing of the input data provided by the user
#                                                                                        
# Description:
# This workflow calls several times 'process-input/verify-input.py', each time with a single option.
# When called this way, 'process-input/verify-input.py' prints a single value. In most cases: 
# 0=input is ok, 1=input is not ok. In other cases: 0=input ok, 1=input has problem,
# 2=input has another problem. This bash program interprets the value and decides whether
# to continue or to exit prematurely.
# Options are: -gnm, -ins, -fq, -gff, -ann, -match
# This workflow also calls 'fasta-merger.py'. The script gets all the fasta files in the input
# folder and concatenates their content, which is written to a new file. The purpose of creating
# an all-contigs-in-one-place file is multiple: First, it eases comparing contig names between
# fasta and gff inputs; second, this file s required for varanalyzer, which is called downstream;
# third, it is easier to run bowtie2-build.
#######################################################################################################


# Check fasta input(s)

fa=`python3 process_input/verify-input.py -gnm $ref_seq  2>> $my_log_file`

if [ $fa == 0 ]; then
	echo $(date "+%F > %T")": Genome fasta input check passed." >> $my_log_file
else
	echo $(date "+%F > %T")": Genome fasta input check failed. One or more genome fasta inputs are empty or have an incorrect format. Please provide new file(s)." >> $my_log_file
	exit_code=1
	echo $exit_code
	exit # Exit because fasta files are required for fasta-concat.py and for fasta-gff comparison
	     # In any other check I do not exit prematurely of process-input.sh so the program checks
	     # all the files and therefore several incorrect files can be revealed to user in a single run
fi


# Concatenate fasta input and store in a single file
{
	python3 process_input/fasta-concat.py -gnm $ref_seq -out_dir $ref_seqs_merged_dir 2>> $my_log_file
} || {
	echo  $(date "+%F > %T")": Processing of fasta genome input failed: fasta-concat.py could not concatenate fasta files into one file." >> $my_log_file
	exit_code=1
}
echo $(date "+%F > %T")": Processing of fasta genome input completed." >> $my_log_file

# Check fasta input with insertion sequence
# Do this only if analyzing insertions

if [ $analysis_type == 'ins' ]; then

	#Insertion fasta cleanup
	{
		python3 process_input/clean-fasta.py -in $ins_seq -out $project_name/$f1/clean-ins.fa 2>> $my_log_file
	} || {
		echo  $(date "+%F > %T")": Insertion fasta cleanup failed." >> $my_log_file
		exit_code=1
	}
	echo $(date "+%F > %T")": Insertion fasta cleanup succeeded." >> $my_log_file

	fa=`python3 process_input/verify-input.py -ins $project_name/$f1/clean-ins.fa 2>> $my_log_file`
	
	if [ $fa == 0 ]; then
		echo $(date "+%F > %T")": Insertion fasta input check passed." >> $my_log_file
	else
		echo $(date "+%F > %T")": Insertion fasta input check failed. File is empty or has an incorrect format. Please provide a new file." >> $my_log_file
		exit_code=1
	fi
fi


# Check fastq input file(s)
# Do this only if data source is exp (reads provided by the user, not simulated)

if [ $data_source == 'exp' ]; then
	if [ $lib_type_sample == 'se' ]; then
		
		fq=`python3 process_input/verify-input.py -fq $read_s 2>> $my_log_file` 
		
		if [ $fq == 0 ]; then
			echo $(date "+%F > %T")": Single-end fastq input (problem reads) passed." >> $my_log_file
		else
			echo $(date "+%F > %T")": Single-end fastq input (problem reads) failed. File is empty or has an incorrect format. Please provide a new file." >> $my_log_file
			exit_code=1
		fi

		fq_qual=`python3 ./graphic_output/fastq-stats.py -fasq $read_s -out $project_name/$f3/single-end-problem-reads-qual-stats.png 2>> $my_log_file` 
		
		if [ $fq_qual == 0 ]; then
			echo $(date "+%F > %T")": Single-end fastq quality (problem reads) encoding is Phred +33. Passed." >> $my_log_file
		else
			echo $(date "+%F > %T")": Single-end fastq quality (problem reads) encoding is not Phred +33. See documentatation to learn how to fix this issue." >> $my_log_file
			exit_code=1
		fi
	fi

	if [ $lib_type_sample == 'pe' ]; then
		fq_for=`python3 process_input/verify-input.py -fq $read_f` 2>> $my_log_file
		
		if [ $fq_for == 0 ]; then
			echo $(date "+%F > %T")": Paired-end forward fastq (problem reads) input passed." >> $my_log_file
		else
			echo $(date "+%F > %T")": Paired-end forward fastq (problem reads) input failed. File is empty or has an incorrect format. Please provide a new file." >> $my_log_file
			exit_code=1
		fi

		fq_rev=`python3 process_input/verify-input.py -fq $read_r` 2>> $my_log_file
		
		if [ $fq_rev == 0 ]; then
			echo $(date "+%F > %T")": Paired-end reverse fastq (problem reads) input passed." >> $my_log_file
		else
			echo $(date "+%F > %T")": Paired-end forward fastq (problem reads) input failed. File is empty or has an incorrect format. Please provide a new file." >> $my_log_file
			exit_code=1
		fi

		fq_qual_for=`python3 ./graphic_output/fastq-stats.py -fasq $read_f -out $project_name/$f3/paired-end-problem-forward-reads-qual-stats.png 2>> $my_log_file` 

		if [ $fq_qual_for == 0 ]; then
			echo $(date "+%F > %T")": Paired-end forward fastq quality (problem reads) encoding is Phred +33. Passed." >> $my_log_file
		else
			echo $(date "+%F > %T")": Paired-end forward fastq quality (problem reads) encoding is not Phred +33. See documentatation to learn how to fix this issue." >> $my_log_file
			exit_code=1
		fi

		fq_qual_rev=`python3 ./graphic_output/fastq-stats.py -fasq $read_r -out $project_name/$f3/paired-end-problem-reverse-reads-qual-stats.png 2>> $my_log_file` 
		
		if [ $fq_qual_rev == 0 ]; then
			echo $(date "+%F > %T")": Paired-end reverse fastq quality (problem reads) encoding is Phred +33. Passed." >> $my_log_file
		else
			echo $(date "+%F > %T")": Paired-end reverse fastq quality (problem reads) encoding is not Phred +33. See documentatation to learn how to fix this issue." >> $my_log_file
			exit_code=1
		fi
	fi

	# If workflow includes control reads, analyze them

	if [ $analysis_type == 'snp' ]; then
		if [ $lib_type_ctrl == 'se' ]; then
			
			fq=`python3 process_input/verify-input.py -fq $read_s_ctrl 2>> $my_log_file` 
			
			if [ $fq == 0 ]; then
				echo $(date "+%F > %T")": Single-end fastq input (control reads) passed." >> $my_log_file
			else
				echo $(date "+%F > %T")": Single-end fastq input (control reads) failed. File is empty or has an incorrect format. Please provide a new file." >> $my_log_file
				exit_code=1
			fi

			{
				fq_qual=`python3 ./graphic_output/fastq-stats.py -fasq $read_s_ctrl -out $project_name/$f3/single-end-control-reads-qual-stats.png 2>> $my_log_file` 
				
				if [ $fq_qual == 0 ]; then
					echo $(date "+%F > %T")": Single-end fastq quality (control reads) encoding is Phred +33. Passed." >> $my_log_file
				else
					echo $(date "+%F > %T")": Single-end fastq quality (control reads) encoding is not Phred +33. See documentatation to learn how to fix this issue." >> $my_log_file
					exit_code=1
				fi
			} || {
				echo $(date "+%F > %T")': fastq-stats.py failed on single-end control reads.' >> $my_log_file
				exit_code=1
			}
		fi

		if [ $lib_type_ctrl == 'pe' ]; then
			fq_for=`python3 process_input/verify-input.py -fq $read_f_ctrl 2>> $my_log_file` 
			
			if [ $fq_for == 0 ]; then
				echo $(date "+%F > %T")": Paired-end forward fastq input (control reads) passed." >> $my_log_file
			else
				echo $(date "+%F > %T")": Paired-end forward fastq input (control reads) failed. File is empty or has an incorrect format. Please provide a new file." >> $my_log_file
				exit_code=1
			fi

			fq_rev=`python3 process_input/verify-input.py -fq $read_r_ctrl 2>> $my_log_file` 
			
			if [ $fq_rev == 0 ]; then
				echo $(date "+%F > %T")": Paired-end reverse fastq input (control reads) passed." >> $my_log_file
			else
				echo $(date "+%F > %T")": Paired-end forward fastq input (control reads) failed. File is empty or has an incorrect format. Please provide a new file." >> $my_log_file
				exit_code=1
			fi

			{
				fq_qual_for=`python3 ./graphic_output/fastq-stats.py -fasq $read_f_ctrl -out $project_name/$f3/paired-end-control-forward-reads-qual-stats.png 2>> $my_log_file` 

				if [ $fq_qual_for == 0 ]; then
					echo $(date "+%F > %T")": Paired-end forward fastq quality (control reads) encoding is Phred +33. Passed." >> $my_log_file
				else
					echo $(date "+%F > %T")": Paired-end forward fastq quality (control reads) encoding is not Phred +33. See documentatation to learn how to fix this issue." >> $my_log_file
					exit_code=1
				fi
			} || {
				echo $(date "+%F > %T")': fastq-stats.py failed on control forward reads file.' >> $my_log_file
				exit_code=1
			}

			{
				fq_qual_rev=`python3 ./graphic_output/fastq-stats.py -fasq $read_r_ctrl -out $project_name/$f3/paired-end-control-reverse-reads-qual-stats.png 2>> $my_log_file` 
				
				if [ $fq_qual_rev == 0 ]; then
					echo $(date "+%F > %T")": Paired-end reverse fastq quality (control reads) encoding is Phred +33. Passed." >> $my_log_file
				else
					echo $(date "+%F > %T")": Paired-end reverse fastq quality (control reads) encoding is not Phred +33. See documentatation to learn how to fix this issue." >> $my_log_file
					exit_code=1
				fi
			} || {
				echo $(date "+%F > %T")': fastq-stats.py failed on reverse reads file.' >> $my_log_file
				exit_code=1
			}

		fi
	fi
fi


# Check gff input
gff=`python3 process_input/verify-input.py -gff $gff_file 2>> $my_log_file` 

if [ $gff == 0 ]; then
	echo $(date "+%F > %T")": GFF3 input check passed." >> $my_log_file
else
	echo $(date "+%F > %T")": GFF3 input check failed. File is empty or has an incorrect format. Please provide a new file." >> $my_log_file
	exit_code=1
fi

# The format or the gene functional annotation file is very flexible. Do not check anything from here for now.
# Check gene funtional annotation input
# Do this only if user has provided file (it's optional)
#if [ $ann_option != 'n/p' ]; then
#	ann=`python3 process_input/verify-input.py -ann $ann_file 2>> $my_log_file` 
#	echo 'Checking gene functional annotation input...' >> $my_log_file
#	if [ $ann == 0 ]; then
#		echo $(date "+%F > %T")": Gene annotation file check passed." >> $my_log_file
#	else
#		echo $(date "+%F > %T")": Gene annotation file check failed. File is empty or has an incorrect format. Please replace input gene functional annotation file or turn off the gene annotation option." >> $my_log_file
#		exit_code=1
#	fi
#fi


# Check contigs match between fasta and gff3 files
match=`python3 process_input/verify-input.py -fa_match $ref_seqs_merged_file -gff_match $gff_file 2>> $my_log_file` 

if [ $match == 0 ]; then
		echo $(date "+%F > %T")": Contigs match check between FASTA and GFF3 inputs passed." >> $my_log_file
elif [ $match == 1 ]; then
		echo $(date "+%F > %T")": Contigs match check between FASTA and GFF3 inputs failed. FASTA input, GFF3 input, or both are empty. Please provide new files." >> $my_log_file
		exit_code=1
elif [ $match == 2 ]; then
		echo $(date "+%F > %T")": Contigs match check between FASTA and GFF3 inputs failed. One or more contig names in the FASTA file are not in the GFF3 file. Please provide new files." >> $my_log_file
		exit_code=1
elif [ $match == 3 ]; then
		echo $(date "+%F > %T")": Contigs match check between FASTA and GFF3 inputs. Warning: some contig names in the GFF3 file are not in the FASTA file. The execution will continue. Please provide new files if considered necessary." >> $my_log_file	
fi

echo $exit_code
