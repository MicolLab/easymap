#!/bin/bash

# When user has selected to simualte data, this bash script is run, which in turn runs
# python programs: sim-mut.py, sim-recsel.py, sim-seq.py.
# If the user chooses to simulate data, hse is not asked to provide reads files.
# If the user chooses workflow insertion analysis, only sim-mut.py and sim-seq.py are
# required because this analysis does not involve a mapping population.
# If the user chooses snp mapping through a recombinant mapping population, in this case
# sim-recsel.py is also executed to create the mapping population from where the reads are
# generated.
#
# There is quite a lot of code repetition in the snp section due to the great amount of possible
# workflows that have to be simulated and to the fact that these workflows have many differences
# between them. I found easier this way to handle this complexity. 
#
# This is the command sent by 'master.sh':
# ./simulator.sh $my_log_file $project_name $analysis_type $lib_type $ins_seq $sim-mut{nbr+mod} $sim-recsel{} $sim-seq{rd+rl+fl+ber+gbs}
# example: ./simulator.sh project/log.log project ins se ins.fa 1+ins {} 10+30,0+0,0+1+50

# Set 'exit_code' (flag variable) to it's initial value (0)
exit_code=0

# Store the location of each folder in a variable
f0=user_data
f1=1_intermediate_files
f2=2_logs

# Get command arguments and assign them to variables
my_log_file=$1
project_name=$2
analysis_type=$3
lib_type=$4
ins_seq=$project_name/$f1/clean-ins.fa
snp_control=${12}

# Write PID to status file
my_status_file=$project_name/$f2/status
echo 'pid simulator '$$ >> $my_status_file

# Get the string that contains the parameters for sim-mut.py and extract them by splitting the string by the '+' character
sim_mut_statement=$6
IFS='+' read -ra sim_mut_array <<< "$sim_mut_statement"
nbr_muts=${sim_mut_array[0]}
mut_mode=${sim_mut_array[1]}

# Get the string that contains the parameters for sim-recsel.py and extract them by splitting the string by the '+' character
sim_recsel_statement=$7
IFS='+' read -ra sim_recsel_array <<< "$sim_recsel_statement"
rec_freq_distr=${sim_recsel_array[0]} # Recombination frequency distribution. Pass it to program as a string and analyze it with python
mut_pos=${sim_recsel_array[1]} # This parameter is also used by in sim_mut
if [[ $mut_pos == *"-"* ]]; then
	second_site_mutagenesis=true
else
	second_site_mutagenesis=false
fi
sel_mode=${sim_recsel_array[2]} # I assume that mutation is always recessive and I select the HM individuals.

# In f2wt mode, for the sample I select the HM individuals and for the control the HZ+WT individuals
nbr_rec_chrs=${sim_recsel_array[3]}

# Get the string that contains the parameters for sim-seq.py and extract them by splitting the string by the '+' character
sim_seq_statement=$8
IFS='+' read -ra sim_seq_array <<< "$sim_seq_statement"
read_depth=${sim_seq_array[0]}
basecalling_error_rate=${sim_seq_array[3]}
gc_bias_strength=${sim_seq_array[4]}

# Split '$rl' (= read length) by the ',' character because it contains a pair of values (mean, sd)
rl=${sim_seq_array[1]}
IFS=',' read -ra rl_array <<< "$rl"
read_length_mean=${rl_array[0]}
read_length_sd=${rl_array[1]}

# Split '$fl' (= fragment length) by the ',' character because it contains a pair of values (mean, sd)
fl=${sim_seq_array[2]}
IFS=',' read -ra fl_array <<< "$fl"
fragment_length_mean=${fl_array[0]}
fragment_length_sd=${fl_array[1]}

# Establish the location of the reference sequence
ref_seqs_merged_file=$project_name/$f1/gnm_ref_merged/genome.fa

# Establish location of some folders required for the simulation
sim_mut_output_folder_ref_lab=$project_name/$f1/sim_data/sim_mut_output/ref-lab_strain
sim_mut_output_folder_noref_lab=$project_name/$f1/sim_data/sim_mut_output/noref-lab_strain
sim_mut_output_folder_mutantstrain=$project_name/$f1/sim_data/sim_mut_output/mutant_strain
# mutated sequence is in: $project_name/$f1/sim_data/sim_mut_output/mutated_genome/mutated_genome.fa
# info is in: $project_name/$f1/sim_data/sim_mut_output/info/info_all_mutations.txt
sim_recsel_output_folder_recessive=$project_name/$f1/sim_data/sim_recsel_output_recessive
sim_recsel_output_folder_dominant=$project_name/$f1/sim_data/sim_recsel_output_dominant
sim_seq_output_folder_sample=$project_name/$f1/sim_data/sim_seq_output/sample
sim_seq_output_folder_control=$project_name/$f1/sim_data/sim_seq_output/control

# Determine if user wants backcross or outcross, the background of the mutant (ref/noref), and the parental to sequence and use as control in 'snp' mode
cross_type=${9}
is_ref_strain=${10}
parental_genome_to_sequence=${11}


# Simulate data for $analisys_type=ins
if [ $analysis_type == 'ins' ]; then
	# Run sim-mut.py
	{
		python3 simulator/sim-mut.py -nbr $nbr_muts -mod $mut_mode -con $ref_seqs_merged_file -ins $ins_seq -out $sim_mut_output_folder_mutantstrain 2>> $my_log_file

	} || {
		echo $(date "+%F > %T")": Simulation of mutagenesis failed. Quit." >> $my_log_file
		exit_code=1; echo $exit_code; exit
	}
	echo $(date "+%F > %T")": Simulation of mutagenesis completed." >> $my_log_file
	
	# Run sim-seq.py. The input is a folder becasuse the program works with all the fasta files that finds in a folder. This is necessary to simulate the sequencing of bulked DNA.
	{
		python3 simulator/sim-seq.py -input_folder $sim_mut_output_folder_mutantstrain/mutated_genome -out $sim_seq_output_folder_sample -mod $lib_type -rd $read_depth -rlm $read_length_mean -rls $read_length_sd -flm $fragment_length_mean -fls $fragment_length_sd -ber $basecalling_error_rate -gbs $gc_bias_strength 2>> $my_log_file

	} || {
		echo $(date "+%F > %T")": Simulation of high-throughput sequencing failed. Quit." >> $my_log_file
		exit_code=1; echo $exit_code; exit
	}
	echo $(date "+%F > %T")": Simulation of high-throughput sequencing completed." >> $my_log_file
fi


# Simulate data for $analisys_type=snp
if [ $analysis_type == 'snp' ]; then

	# Calculate genome length to know how many natural SNPs to introduce
	{
		genome_length=`python3 simulator/calculate-genome-length.py -gnm $ref_seqs_merged_file 2>> $my_log_file`

	} || {
		echo $(date "+%F > %T")": simulator/calculate-genome-length.py failed. Quit." >> $my_log_file
		exit_code=1; echo $exit_code; exit
	}
	echo $(date "+%F > %T")": simulator/calculate-genome-length.py finished." >> $my_log_file

	# Calculate how many natural SNPs to introduce in the ref-lab and the noref-lab strains.
	# The amount is based on observed natural mutations between NCBI reference sequence and
	#	- Columbia-0 (reference strain): 0.001% bases of the genome contain an SNP
	#	- Landsberg (non-reference strain): 0.42% bases of the genome contain an SNP
	nbr_natural_mutations_ref=`python3 -c "print(int(round($genome_length * 0.000010084)))" 2>> $my_log_file`
	nbr_natural_mutations_noref=`python3 -c "print(int(round($genome_length * 0.0042016807)))" 2>> $my_log_file`

	# Using as input the reference sequence provided by user, simulate ref-lab and noref-lab sequences 
	
	# Run sim-mut.py to create ref-lab strain. Mutate 0.001% of bases.
	{
		python3 simulator/sim-mut.py -nbr $nbr_natural_mutations_ref -mod d -con $ref_seqs_merged_file -out $sim_mut_output_folder_ref_lab 2>> $my_log_file

	} || {
		echo $(date "+%F > %T")": Simulation of mutagenesis to ref-lab strain failed. Quit." >> $my_log_file
		exit_code=1; echo $exit_code; exit
	}
	echo $(date "+%F > %T")": Simulation of mutagenesis to create ref-lab strain completed." >> $my_log_file

	# Run sim-mut.py to create noref-lab strain. Mutate 0.4% of bases.
	{
		python3 simulator/sim-mut.py -nbr $nbr_natural_mutations_noref -mod d -con $ref_seqs_merged_file -out $sim_mut_output_folder_noref_lab 2>> $my_log_file

	} || {
		echo $(date "+%F > %T")": Simulation of mutagenesis to noref-lab strain failed. Quit." >> $my_log_file
		exit_code=1; echo $exit_code; exit
	}
	echo $(date "+%F > %T")": Simulation of mutagenesis to create noref-lab strain completed." >> $my_log_file
	
	# Create the mutant sequence starting from one of the lab strains
	if [ $is_ref_strain == 'ref' ]; then
		template=$sim_mut_output_folder_ref_lab/mutated_genome/mutated_genome.fa
	else
		template=$sim_mut_output_folder_noref_lab/mutated_genome/mutated_genome.fa
	fi
	
	# Run sim-mut.py to create mutant strain
	mut_pos_1=$(echo $mut_pos | cut -d'-' -f 1) #I'm adding this just in case we are dealing with a second site mutagenesis. Fist mutagenesis happens here, while the second is below.
	 
	{
		python3 simulator/sim-mut.py -nbr $nbr_muts -mod $mut_mode -con $template -out $sim_mut_output_folder_mutantstrain -causal_mut $mut_pos_1 2>> $my_log_file

	} || {
		echo $(date "+%F > %T")": Simulation of mutagenesis to create the mutant strain failed. Quit." >> $my_log_file
		exit_code=1; echo $exit_code; exit
	}
	echo $(date "+%F > %T")": Simulation of mutagenesis to create the mutant strain completed." >> $my_log_file

	
	# Define location of mutant genome
	parmut_sample=$sim_mut_output_folder_mutantstrain/mutated_genome/mutated_genome.fa


	# Perform second site mutagenesis
	sim_mut_output_folder_mutantstrain2=$project_name/$f1/sim_data/sim_mut_output/mutant_strain_2
	
	if [ $second_site_mutagenesis == true ]; then

		mut_pos=$(echo $mut_pos | cut -d'-' -f 2)

		{
			python3 simulator/sim-mut.py -nbr $nbr_muts -mod $mut_mode -con $parmut_sample -out $sim_mut_output_folder_mutantstrain2 -causal_mut $mut_pos 2>> $my_log_file

		} || {
			echo $(date "+%F > %T")": Simulation of second site mutagenesis to create the mutant strain failed. Quit." >> $my_log_file
			exit_code=1; echo $exit_code; exit
		}
		echo $(date "+%F > %T")": Simulation of second site mutagenesis to create the mutant strain completed." >> $my_log_file

		parmut_sample=$sim_mut_output_folder_mutantstrain2/mutated_genome/mutated_genome.fa

		sim_mut_output_folder_ref_lab=$sim_mut_output_folder_mutantstrain
		sim_mut_output_folder_noref_lab=$sim_mut_output_folder_mutantstrain
		sim_mut_output_folder_mutantstrain=$sim_mut_output_folder_mutantstrain2
	fi


	# Create recombinant F2 population(s)
	# If control sample is going to be one of the parentals, simply simulate the recessive phenotype F2 population.
	# If the control sample is going to be the dominant phenotype population, simulate this one also.
	
	# Define locations of polymorphic samples, necessary to run sim-recsel.py
	if [ $is_ref_strain == 'ref' ] && [ $cross_type == 'bc' ]; then
		parpol_sample=$sim_mut_output_folder_ref_lab/mutated_genome/mutated_genome.fa
	elif [ $is_ref_strain == 'ref' ] && [ $cross_type == 'oc' ]; then
		parpol_sample=$sim_mut_output_folder_noref_lab/mutated_genome/mutated_genome.fa
	elif [ $is_ref_strain == 'noref' ] && [ $cross_type == 'bc' ]; then
		parpol_sample=$sim_mut_output_folder_noref_lab/mutated_genome/mutated_genome.fa
	elif [ $is_ref_strain == 'noref' ] && [ $cross_type == 'oc' ]; then
		parpol_sample=$sim_mut_output_folder_ref_lab/mutated_genome/mutated_genome.fa
	fi
	
	if [ $snp_control == 'par' ]; then
		
		# Run sim-recsel.py to create recombinant chromosomes selected to carry the mutation
		{
			python3 simulator/sim-recsel.py -outdir $sim_recsel_output_folder_recessive -rec_freq_distr $rec_freq_distr -parmut $parmut_sample -parpol $parpol_sample -mutpos $mut_pos -smod $sel_mode -nrec $nbr_rec_chrs 2>> $my_log_file 

		} || {
			echo $(date "+%F > %T")": Simulation of recombination and phenotype selection failed. Quit." >> $my_log_file
			exit_code=1; echo $exit_code; exit
		}
		echo $(date "+%F > %T")": Simulation of recombination and phenotype selection completed." >> $my_log_file
		
	else #f2wt
		
		# Run sim-recsel.py to create the F2 recessive population
		{
			python3 simulator/sim-recsel.py -outdir $sim_recsel_output_folder_recessive -rec_freq_distr $rec_freq_distr -parmut $parmut_sample -parpol $parpol_sample -mutpos $mut_pos -smod $sel_mode -nrec $nbr_rec_chrs 2>> $my_log_file
		
		} || {
			echo $(date "+%F > %T")": Simulation of recombination and phenotype selection to create the F2 recessive population failed. Quit." >> $my_log_file
			exit_code=1; echo $exit_code; exit
		}
		echo $(date "+%F > %T")": Simulation of recombination and phenotype selection to create the F2 recessive population completed." >> $my_log_file
		
		# Run sim-recsel.py to create the F2 dominant population
		{
			python3 simulator/sim-recsel.py -outdir $sim_recsel_output_folder_dominant -rec_freq_distr $rec_freq_distr -parmut $parmut_sample -parpol $parpol_sample -mutpos $mut_pos -smod wt -nrec $nbr_rec_chrs 2>> $my_log_file
		
		} || {
			echo $(date "+%F > %T")": Simulation of recombination and phenotype selection to create the F2 dominant population failed. Quit." >> $my_log_file
			exit_code=1; echo $exit_code; exit
		}
		echo $(date "+%F > %T")": Simulation of recombination and phenotype selection to create the F2 dominant population completed." >> $my_log_file
	fi
	

	# Simulate reads from the sample and the control templates
	if [ $snp_control == 'par' ]; then
		if [ $is_ref_strain == 'ref' ] && [ $parental_genome_to_sequence == 'mut' ]; then
			input_folder_control=$sim_mut_output_folder_ref_lab/mutated_genome
		elif [ $is_ref_strain == 'ref' ] && [ $parental_genome_to_sequence == 'nomut' ]; then
			input_folder_control=$sim_mut_output_folder_noref_lab/mutated_genome
		elif [ $is_ref_strain == 'noref' ] && [ $parental_genome_to_sequence == 'mut' ]; then
			input_folder_control=$sim_mut_output_folder_noref_lab/mutated_genome
		elif [ $is_ref_strain == 'noref' ] && [ $parental_genome_to_sequence == 'nomut' ]; then
			input_folder_control=$sim_mut_output_folder_ref_lab/mutated_genome
		fi
		
	else #f2wt
		input_folder_control=$sim_recsel_output_folder_dominant
	fi
	
	
	# Run sim-seq.py on control genome. The input is a folder because the program works with all the fasta files that finds in a folder. This is necessary to simulate the sequencing of bulked DNA.
	{
		python3 simulator/sim-seq.py -input_folder $input_folder_control -out $sim_seq_output_folder_control -mod $lib_type -rd $read_depth -rlm $read_length_mean -rls $read_length_sd -flm $fragment_length_mean -fls $fragment_length_sd -ber $basecalling_error_rate -gbs $gc_bias_strength 2>> $my_log_file

	} || {
		echo $(date "+%F > %T")": Simulation of high-throughput sequencing reads on control genome failed. Quit." >> $my_log_file
		exit_code=1; echo $exit_code; exit
	} 
	echo $(date "+%F > %T")": Simulation of high-throughput sequencing reads on control genome completed." >> $my_log_file

	# Run sim-seq.py on F2 recombinant population. The input is a folder becasuse the program works with all the fasta files that finds in a folder. This is necessary to simulate the sequencing of bulked DNA.
	{
		python3 simulator/sim-seq.py -input_folder $sim_recsel_output_folder_recessive -out $sim_seq_output_folder_sample -mod $lib_type -rd $read_depth -rlm $read_length_mean -rls $read_length_sd -flm $fragment_length_mean -fls $fragment_length_sd -ber $basecalling_error_rate -gbs $gc_bias_strength 2>> $my_log_file

	} || {
		echo $(date "+%F > %T")": Simulation of high-throughput sequencing on F2 recombinant population failed. Quit." >> $my_log_file
		exit_code=1; echo $exit_code; exit
	}
	echo $(date "+%F > %T")": Simulation of high-throughput sequencing reads on F2 recombinant population completed." >> $my_log_file

fi

echo $exit_code
