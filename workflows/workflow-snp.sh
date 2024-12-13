#!/bin/bash
#
# This is the command sent by 'master.sh':
# ./workflow-x.sh $my_log_file $project_name $workflow $data_source $lib_type $ins_seq $read_s $read_f $read_r $gff_file $ann_file
#
# This command is always the same, regardless of the workflow
#
# Fields that can be equal to 'n/p' (= not provided) are the following: $ins_seq, $read_s, $read_f, $read_r, $ann_file
# If $ins_seq = n/p, that is because $workflow = snp. Therefore, no insertion seq is needed and $ins_seq is ignored by hte program.
# If $read_s = n/p, that is because $lib_type = pe, so it is ignored by the program.
# If $read_f and $read_r = n/p, that is because $lib_type = se, so it is ignored by the program.
# If $ann_file = n/p, this is because user did no have it. This program has the deal with this: if data not provided, simply do not
# include gene annotated info to the report.
#
# my_log_file		>	$1
# project_name		>	$2
# workflow			>	$3
# data_source		>	$4
# lib_type			>	$5
# ins_seq			>	$6
# read_s			>	$7
# read_f			>	$8
# read_r			>	$9
# gff_file			>	${10}
# ann_file			>	${11}
# read_s_par							>	${12}
# read_f_par							>	${13}
# read_r_par							>	${14}
# cross_type							>	${15} 
# is_ref_strain							>	${16} 
# parental_reads_provided				>	${17}
# snp_analysis_type [par/f2wt]			>	${18}
# lib_type_control						>	${19}
# stringency							>	${20}


# Some initial parameters
start_time=`date +%s`	
exit_code=0 				# Set 'exit_code' (flag variable) to 0
my_log_file=$1 				# Set location of log file
export location="$PWD" 			#Save path to hisat2-build and hisat2 in variable BT2

# Create input variables
my_log_file=$1
project_name=$2
my_sample_mode=$5 												#[pe, se], paired/single  
my_control_mode=${19}											#TEMPORAL, in future independent						
my_rd=$7											 			#reads (single)
my_rf=$8 														#forward reads
my_rr=$9												 		#reverse reads 			
my_p_rd=${12}											 		#reads (single) parent	
my_p_rf=${13} 													#forward reads parent	
my_p_rr=${14}											 		#reverse reads parent	
my_gs=gnm_ref_merged/genome.fa 									#genome sequence
my_ix=genome_index 							
my_gff=${10}													#Genome feature file
my_ann=${11}													
my_rrl=250 														#Regulatory region length
my_log_file=$1
my_mut=snp  													#my_mut takes the values 'snp' in this workflow and 'lin' in the large insertions workflow, for the execution of the graphic output module
my_cross=${15}													#oc / bc : f2 obtained by outcross or backcross 									
my_mutbackgroud=${16}											#ref / noref : genetic background of the mutation									
my_pseq=${17}													#mut / nomut : sequenced parental provided is the mutagenized one or the other		
snp_analysis_type=${18}
stringency=${20}

#Set number of maximum CPU for steps compatible with multithreading, default = 1 
threads=1

# Set internal variables according to the SNP validation stringency chosen by the user
if [ $stringency == high_stringency ]; then
	#problemSample_bowtie_mp="--mp 6,2"
	problemSample_mpileup_C="-C50"
	problemSample_snpQualityTheshold="120"
else
	#problemSample_bowtie_mp="--mp 3,2"
	problemSample_mpileup_C=""
	problemSample_snpQualityTheshold="30"
fi

# Define the folders in the easymap directory 
f0=user_data
f1=$project_name/1_intermediate_files
f2=$project_name/2_logs
f3=$project_name/3_workflow_output

# Write PID to status file
my_status_file=$f2/status
echo 'pid workflow '$$ >> $my_status_file

#Check genome size to set interval_width
{
	interval_width=`python3 $location/scripts_snp/set-interval.py -a $f1/$my_gs`
} || {
	interval_width=4000001
	echo $(date "+%F > %T")': set-interval.py failed.' >> $my_log_file
}
echo $(date "+%F > %T")': set-interval.py finished, interval set at: '$interval_width   >> $my_log_file


##################################################################################################################################################################################
#																																												 #
#																																												 #
#																	F2 FQ PROCESSING FUNCTION																					 #
#																																												 #
#																																												 #
##################################################################################################################################################################################

#Run hisat2-build on genome sequence 
{
	$location/hisat2/hisat2-build $f1/$my_gs $f1/$my_ix 1> $f2/hisat2-build_std1.txt 2> $f2/hisat2-build_std2.txt

} || {
	echo $(date "+%F > %T")': hisat2-build on genome sequence returned an error. See log files.' >> $my_log_file
	exit_code=1
	echo $exit_code
	exit
}
echo $(date "+%F > %T")': hisat2-build finished.' >> $my_log_file

function get_problem_va {  
	if [ $my_sample_mode == se ] 
	then
		#Run hisat2 unpaired to align raw F2 reads to genome 
		{
			$location/hisat2/hisat2 -p $threads -x $f1/$my_ix -U $my_rd -S $f1/alignment1.sam 2> $f2/hisat2_problem-sample_std2.txt

		} || {
			echo $(date "+%F > %T")': hisat2 returned an error during the aligment of F2 reads. See log files.' >> $my_log_file
			exit_code=1
			echo $exit_code
			exit
		}
		echo $(date "+%F > %T")': hisat2 finished the alignment of F2 reads to genome.' >> $my_log_file
	fi

	if [ $my_sample_mode == pe ] 
	then
		#Run hisat2 paired to align raw F2 reads to genome 
		{
			$location/hisat2/hisat2  -p $threads  -x $f1/$my_ix -1 $my_rf -2 $my_rr -S $f1/alignment1.sam 2> $f2/hisat2_problem-sample_std2.txt

		} || {
			echo $(date "+%F > %T")': hisat2 returned an error during the aligment of F2 reads. See log files.' >> $my_log_file
			exit_code=1
			echo $exit_code
			exit
		}
		echo $(date "+%F > %T")': hisat2 finished the alignment of F2 reads to genome.' >> $my_log_file
	fi

	#SAM to BAM
	{
		$location/samtools1/samtools sort  -@ $threads  $f1/alignment1.sam > $f1/alignment1.bam 2> $f2/sam-to-bam_problem-sample_std2.txt
		rm -rf ./user_projects/$project_name/1_intermediate_files/alignment1.sam

	} || {
		echo 'Error transforming SAM to BAM.' >> $my_log_file
		exit_code=1
		echo $exit_code
		exit
	}
	echo $(date "+%F > %T")': SAM to BAM finished.' >> $my_log_file

	#Variant calling
	{

		$location/samtools1/samtools mpileup  -t DP,ADF,ADR $problemSample_mpileup_C -uf $f1/$my_gs $f1/alignment1.bam 2> $f2/mpileup_problem-sample_std.txt | $location/bcftools-1.3.1/bcftools call -mv -Ov > $f1/raw_variants.vcf 2> $f2/call_problem-sample_std.txt
		# -B: Disables probabilistic realignment for the computation of base alignment quality (BAQ). Applying this argument reduces the number of false negatives during the variant calling
		# -t DP,ADF,ADR: output VCF file contains the specified optional columns: read depth (DP), allelic depths on the forward strand (ADF), allelic depths on the reverse strand (ADR)
		# -uf: uncompressed vcf output / fasta imput genome file
		# -mv: include only polymorphic sites in output
		# -Ov: uncompressed VCF output file 
		# -C50: reduce the effect of reads with excessive mismatches. This aims to fix overestimated mapping quality

	} || {
		echo $(date "+%F > %T")': Error during variant-calling of F2 data.' >> $my_log_file
		exit_code=1
		echo $exit_code
		exit
	}
	echo $(date "+%F > %T")': F2 data variant calling finished.' >> $my_log_file

	#Groom vcf
	{
		python3 $location/scripts_snp/vcf-groomer.py -a $f1/raw_variants.vcf -b $f1/F2_raw.va  2>> $my_log_file

	} || {
		echo $(date "+%F > %T")': Error during execution of vcf-groomer.py with F2 data.' >> $my_log_file
		exit_code=1
		echo $exit_code
		exit
	}
	echo $(date "+%F > %T")': VCF grooming of F2 data finished.' >> $my_log_file

	#RD graphics
	depth_alignment $f1/alignment1.bam $f3/frequence_depth_alignment_distribution_sample.png

	#Run vcf filter
	if [ $my_cross == bc ]; then mut_type=EMS ; fi
	if [ $my_cross == oc ]; then mut_type=all ; fi

	if [ $av_rd -gt 25 ]; then dp_min=15 ; fi
	if [ $av_rd -le 25 ]; then dp_min=10 ; fi

	dp_max=$(($av_rd * 3))
	if [ $dp_max -le 40 ]; then dp_max=100 ; fi


	{
		python3 $location/scripts_snp/variants-filter.py -a $f1/F2_raw.va -b $f1/F2_filtered.va -step 3 -fasta $f1/$my_gs -dp_min $dp_min -dp_max $dp_max -qual_min $problemSample_snpQualityTheshold -mut_type $mut_type  2>> $my_log_file

	} || {
		echo 'Error during execution of variants-filter.py with F2 data.' >> $my_log_file
		exit_code=1
		echo $exit_code
		exit
	}
	echo $(date "+%F > %T")': First VCF filtering step of F2 data finished.' >> $my_log_file

	#Intermediate files cleanup
	rm -f $f1/*.sam
	rm -f $f1/*.vcf

}


##################################################################################################################################################################################
#																																												 #
#																																												 #
#																	control FQ PROCESSING FUNCTION																				 #
#																																												 #
#																																												 #
##################################################################################################################################################################################

function get_control_va { 
	if [ $my_control_mode == se ] 
	then
		#Run hisat2 unpaired to align raw F2 reads to genome 
		{
			$location/hisat2/hisat2  -p $threads  -x $f1/$my_ix -U $my_p_rd -S $f1/alignment1P.sam 2> $f2/hisat2_control-sample_std2.txt

		} || {
			echo $(date "+%F > %T")': hisat2 returned an error during the aligment of control reads. See log files.' >> $my_log_file
			exit_code=1
			echo $exit_code
			exit
		}
		echo $(date "+%F > %T")': hisat2 finished the alignment of control reads to genome.' >> $my_log_file
	fi

	if [ $my_control_mode == pe ] 
	then
		#Run hisat2 paired to align raw F2 reads to genome 
		{
			$location/hisat2/hisat2  -p $threads  -x $f1/$my_ix -1 $my_p_rf -2 $my_p_rr -S $f1/alignment1P.sam 2> $f2/hisat2_control-sample_std2.txt

		} || {
			echo $(date "+%F > %T")': hisat2 returned an error during the aligment of control reads. See log files.' >> $my_log_file
			exit_code=1
			echo $exit_code
			exit
		}
		echo $(date "+%F > %T")': hisat2 finished the alignment of control reads to genome.' >> $my_log_file
	fi

	#SAM to BAM
	{
		$location/samtools1/samtools sort  -@ $threads $f1/alignment1P.sam > $f1/alignment1P.bam 2> $f2/sam-to-bam_control-sample_std2.txt

		rm -rf ./user_projects/$project_name/1_intermediate_files/alignment1P.sam

	} || {
		echo $(date "+%F > %T")': Error transforming SAM to BAM' >> $my_log_file
		exit_code=1
		echo $exit_code
		exit
	}
	echo $(date "+%F > %T")': SAM to BAM finished' >> $my_log_file

	#Variant calling
	{

		$location/samtools1/samtools mpileup  -t DP,ADF,ADR  -uf $f1/$my_gs $f1/alignment1P.bam 2> $f2/mpileup_control-sample_std.txt | $location/bcftools-1.3.1/bcftools call -mv -Ov > $f1/raw_p_variants.vcf 2> $f2/call_control-sample_std.txt

	} || {
		echo $(date "+%F > %T")': Error during variant-calling of control data' >> $my_log_file
		exit_code=1
		echo $exit_code
		exit
	}
	echo $(date "+%F > %T")': Control data variant calling finished' >> $my_log_file

	#Groom vcf
	{
		python3 $location/scripts_snp/vcf-groomer.py -a $f1/raw_p_variants.vcf -b $f1/control_raw.va  2>> $my_log_file

	} || {
		echo $(date "+%F > %T")': Error during execution of vcf-groomer.py with control data.' >> $my_log_file
		exit_code=1
		echo $exit_code
		exit
	}
	echo $(date "+%F > %T")': VCF grooming of control data finished.' >> $my_log_file

	#RD graphics
	depth_alignment $f1/alignment1P.bam $f3/frequence_depth_alignment_distribution_control.png

	#Run vcf filter

	if [ $av_rd -gt 25 ]; then dp_min=15 ; fi
	if [ $av_rd -le 25 ]; then dp_min=10 ; fi

	dp_max=$(($av_rd * 3))
	if [ $dp_max -le 40 ]; then dp_max=100 ; fi

	{
		python3 $location/scripts_snp/variants-filter.py -a $f1/control_raw.va -b $f1/control_filtered.va -step 3 -fasta $f1/$my_gs -dp_min 10 -dp_max $dp_max -qual_min 20  2>> $my_log_file

	} || {
		echo $(date "+%F > %T")': Error during execution of variants-filter.py with control data.' >> $my_log_file
		exit_code=1
		echo $exit_code
		exit
	}
	echo $(date "+%F > %T")': First VCF filtering step of control data finished.' >> $my_log_file

	#Intermediate files cleanup
	rm -f $f1/*.sam
	rm -f $f1/*.vcf
}


##################################################################################################################################################################################
#																																												 #
#																																												 #
#																			DEPTH ALIGNMENT ANALYSIS FUNCTION																	 #
#																																												 #
#																																												 #
##################################################################################################################################################################################

function depth_alignment {
	{
		python3 $location/scripts_snp/depth_measures_generation.py -genome $f1/$my_gs -bam $1 -out $f1/coverage_alignment1.txt  2>> $my_log_file
		rm -rf ./user_projects/$project_name/1_intermediate_files/alignment1.bam

	} || {
		echo $(date "+%F > %T")': Error during obtaining of alignment depth .' >> $my_log_file
		#exit_code=1
		#echo $exit_code
		#exit
	}

	{
		av_rd=`python3 $location/graphic_output/graphic-alignment.py -coverages $f1/coverage_alignment1.txt   -out $2  2>> $my_log_file `

	} || {
		echo $(date "+%F > %T")': Error during Graphic_alignment execution in sample alignment.' >> $my_log_file
		av_rd=20			#SDL - simple bypass
		#exit_code=1
		#echo $exit_code
		#exit
	}
}


##################################################################################################################################################################################
#																																												 #
#																																												 #
#																	CANDIDATE REGION ANALYSIS FUNCTION																			 #
#																																												 #
#																																												 #
##################################################################################################################################################################################

function cr_analysis {

	# Run vcf filter, selecting snps in the candidate region defined by map-mutation.py, with an alelic frequence > 0.8 and corresponding to EMS mutations
	{
		python3 $location/scripts_snp/variants-filter.py -a $f1/F2_control_comparison.va -b $f1/final_variants.va -step 2 -cand_reg_file $f1/map_info.txt -af_min 0.8 -mut_type EMS  2>> $my_log_file

	} || {
		echo $(date "+%F > %T")': Error during the second execution of variants-filter.py .' >> $my_log_file
		exit_code=1
		echo $exit_code
		exit
	}
	echo $(date "+%F > %T")': Second VCF filtering step finished.' >> $my_log_file

	# Create input for varanalyzer and run varanalyzer.py (one file for the candidate region and one for the whole genome)
	{
		python3 $location/scripts_snp/snp-to-varanalyzer.py -a $f1/final_variants.va -b $f1/snp-to-varanalyzer.txt  2>> $my_log_file
		python3 $location/scripts_snp/snp-to-varanalyzer.py -a $f1/F2_control_comparison.va -b $f1/snp-to-varanalyzer-total.txt  2>> $my_log_file

	} || {
		echo $(date "+%F > %T")': Error during execution of snp-to-varanalyzer.py .' >> $my_log_file
		exit_code=1
		echo $exit_code
		exit
	}
	echo $(date "+%F > %T")': Input for varanalyzer finished.' >> $my_log_file
	# Varanalyzer
	{
		python3 $location/varanalyzer/varanalyzer.py -itp snp -con $f1/$my_gs -gff $f0/$my_gff -var $f1/snp-to-varanalyzer.txt -rrl $my_rrl -pname $project_name -ann $f0/$my_ann -out $f1/varanalyzer_output.txt 2>> $my_log_file
		python3 $location/varanalyzer/varanalyzer.py -itp snp -con $f1/$my_gs -gff $f0/$my_gff -var $f1/snp-to-varanalyzer-total.txt -rrl $my_rrl -pname $project_name -ann $f0/$my_ann -out $f1/varanalyzer_output_total.txt  2>> $my_log_file

	} || {
		echo $(date "+%F > %T")': Error during execution of varanalyzer.py .' >> $my_log_file
		exit_code=1
		echo $exit_code
		exit
	}
	echo $(date "+%F > %T")': Varanalyzer finished.' >> $my_log_file

	# Run primer generation script
	{
		python3 $location/primers/primer-generation.py -file $f1/varanalyzer_output.txt -fasta $f1/$my_gs -out $f1/primer_generation_output.txt  -mode 2   2>> $my_log_file
		python3 $location/primers/primer-generation.py -file $f1/varanalyzer_output_total.txt -fasta $f1/$my_gs -out $f1/primer_generation_output_total.txt  -mode 2   2>> $my_log_file

	}|| {
		echo $(date "+%F > %T")': primer-generation.py failed.'>> $my_log_file
	}
	echo $(date "+%F > %T")': primer-generation.py finished.' >> $my_log_file
	
	# Run extend-snp-variants-info                                              --project-name $project_name
	result_extend_snp_info=`python3 $location/scripts_snp/extend-snp-variants-info.py  --variants $f1/primer_generation_output.txt --snp-info $f1/snp-to-varanalyzer.txt --project-name $project_name --map-info $f1/map_info.txt --output-file $f3/candidate_variants.txt --region CR 2>> $my_log_file`
	result_extend_snp_info=`python3 $location/scripts_snp/extend-snp-variants-info.py  --variants $f1/primer_generation_output_total.txt --snp-info $f1/snp-to-varanalyzer-total.txt --project-name $project_name --map-info $f1/map_info.txt --output-file $f3/candidate_variants_total.txt --region total 2>> $my_log_file`
	
	if [ $result_extend_snp_info == 'success' ]; then
		echo $(date "+%F > %T")": extend-snp-variants-info.py finished." >> $my_log_file
	elif [ $result_extend_snp_info == 'error' ]; then
		echo $(date "+%F > %T")": Error: extend-snp-variants-info.py failed." >> $my_log_file
		exit_code=1
		echo $exit_code
		exit
	fi
	
	# Filter SNPs to draw
	{
		python3 $location/scripts_snp/variants-filter.py -a $f1/$1 -b $f1/F2_control_comparison_drawn.va -step 1 -af_min $2   2>> $my_log_file

	} || {
		echo $(date "+%F > %T")': Error during third execution of variants-filter.py . ' >> $my_log_file
		exit_code=1
		echo $exit_code
		exit
	}
	echo $(date "+%F > %T")': Third VCF filtering step finished.' >> $my_log_file

	# Draw candidates 
	{
		python3 $location/graphic_output/graphic-output.py -my_mut af_candidates -asnp $f1/F2_control_comparison.va -bsnp $f1/$my_gs -rrl $my_rrl -iva $2/1_intermediate_files/varanalyzer_output.txt -gff $f0/$my_gff -pname $project_name  -cross $my_cross -snp_analysis_type $snp_analysis_type  2>> $my_log_file
		
	} || {
		echo $(date "+%F > %T")': Error during execution of graphic-output.py .' >> $my_log_file
		exit_code=1
		echo $exit_code
		exit
	}
	echo $(date "+%F > %T")': Graphic output created.' >> $my_log_file

	# python3 ./graphic_output/graphic-output.py -my_mut snp -asnp ./user_projects/project/1_intermediate_files/F2_control_comparison_drawn.va -bsnp ./user_projects/project/1_intermediate_files/gnm_ref_merged/genome.fa -rrl 150 -iva ./user_projects/project/1_intermediate_files/varanalyzer_output.txt -gff ./user_data/complete.gff -pname user_projects/project  -cross bc -snp_analysis_type par  
	# (6) Create graphic output
	{
		python3 $location/graphic_output/graphic-output.py -my_mut $my_mut  -interval_width $interval_width  -asnp $f1/F2_control_comparison_drawn.va -bsnp $f1/$my_gs -rrl $my_rrl -iva $project_name/1_intermediate_files/varanalyzer_output.txt -gff $f0/$my_gff -pname $project_name  -cross $my_cross -snp_analysis_type $snp_analysis_type  2>> $my_log_file
		
	} || {
		echo $(date "+%F > %T")': Error during execution of graphic-output.py .' >> $my_log_file
		exit_code=1
		echo $exit_code
		exit
	}
	echo $(date "+%F > %T")': Graphic output created.' >> $my_log_file

	# (7) Create report
	#Some arrangements

	{
		cp $location/fonts/legend.png $f3/legend.png
		cp $location/fonts/gene_legend_snp.png $f3/gene_legend_snp.png
		zip $f3/report_images.zip $f3/*.png > $f2/zip.txt
	} || {
		echo $(date "+%F > %T")': Error during zip compression of report images. Please check that the zip program is installed in your system. Continuing anyway.' >> $my_log_file
	}

	{
		python3 $location/graphic_output/report.py -files_dir $f3 -variants $f3/candidate_variants.txt -log $f2/log.log -output_html $f3/report.html -project $project_name -mut_type $my_mut  2>> $my_log_file
		
	} || {
		echo $(date "+%F > %T")': Error during report generation.' >> $my_log_file
		exit_code=1
		echo $exit_code
		exit
	}
	echo $(date "+%F > %T")': Report file created.' >> $my_log_file
}



##################################################################################################################################################################################
#																																												 #
#																																												 #
#																			DATA ANALYSIS 																						 #
#																																												 #
#																																												 #
##################################################################################################################################################################################

#_________________________________________________________________________________________________________________________________________________________________________________

#_________________________________Case 1 and 5: Mutant in ref/noref background, backcross, mutant parental control (ref/noref bc mut)_______________________________________________________________
#_________________________________________________________________________________________________________________________________________________________________________________

if [ $my_pseq == mut ] && [ $my_cross == bc ]  && [ $snp_analysis_type == par ]
then

	# (1) Get problem and control VA files
	get_problem_va 
	get_control_va

	#draw snps
	python3 $location/graphic_output/graphic-output.py -my_mut af_control -asnp $f1/control_filtered.va -bsnp $f1/$my_gs -rrl $my_rrl -iva $2/1_intermediate_files/varanalyzer_output.txt -gff $f0/$my_gff -pname $2  -cross $my_cross -snp_analysis_type $snp_analysis_type  2>> $my_log_file
	python3 $location/graphic_output/graphic-output.py -my_mut af_sample -asnp $f1/F2_filtered.va -bsnp $f1/$my_gs -rrl $my_rrl -iva $2/1_intermediate_files/varanalyzer_output.txt -gff $f0/$my_gff -pname $2  -cross $my_cross -snp_analysis_type $snp_analysis_type  2>> $my_log_file

	# (2) Run VA operations: Remove control SNPs from problem file
	my_operation_mode=A
	{
		python3 $location/scripts_snp/variants-operations.py -a $f1/F2_filtered.va -b $f1/control_filtered.va -c $f1/F2_control_comparison.va -mode $my_operation_mode -primary 1  2>> $my_log_file

	} || {
		echo $(date "+%F > %T")': Error during first execution of variants-operations.py .' >> $my_log_file
		exit_code=1
		echo $exit_code
		exit
	}
	echo $(date "+%F > %T")': VCF operations finished.' >> $my_log_file

	# (3) Run mapping analysis
	my_analysis_mode=back

	if [ $(wc -l < $f1/F2_control_comparison.va) -gt 1 ]
	then 
		{
			python3 $location/scripts_snp/map-mutation.py -file $f1/F2_control_comparison.va -fasta $f1/$my_gs -mode $my_analysis_mode -window_size 600000 -window_space 500000 -output $f1/map_info.txt -control_modality $my_mutbackgroud -interval_width $interval_width -snp_analysis_type $snp_analysis_type  2>> $my_log_file

		} || {
			echo $(date "+%F > %T")': Error during execution of map-mutation.py .' >> $my_log_file
			exit_code=1
			echo $exit_code
			exit
		}
		echo $(date "+%F > %T")': Mutation mapping module finished.' >> $my_log_file

		# (4) Candidate region analysis function
		cr_analysis F2_control_comparison.va 0.1

	else 
		echo $(date "+%F > %T")': No informative polymorphisms can be identified for mutation mapping. Please review your input data.' >> $my_log_file
		exit_code=1
		echo $exit_code
		exit

	fi
fi


#_________________________________________________________________________________________________________________________________________________________________________________
#
#_________________________________Case 2 and case 6: Backcross, f2wt control (bc mut)_________________________________________________________________________________________
#_________________________________________________________________________________________________________________________________________________________________________________

#if [ $my_mutbackgroud == ref ] && [ $my_cross == bc ]  && [ $snp_analysis_type == f2wt ]
if [ $my_cross == bc ]  && [ $snp_analysis_type == f2wt ]
then

	# (1) Get problem and control VA files
	get_problem_va
	get_control_va

	#draw snps
	python3 $location/graphic_output/graphic-output.py -my_mut af_control -asnp $f1/control_filtered.va -bsnp $f1/$my_gs -rrl $my_rrl -iva $2/1_intermediate_files/varanalyzer_output.txt -gff $f0/$my_gff -pname $2  -cross $my_cross -snp_analysis_type $snp_analysis_type  2>> $my_log_file
	python3 $location/graphic_output/graphic-output.py -my_mut af_sample -asnp $f1/F2_filtered.va -bsnp $f1/$my_gs -rrl $my_rrl -iva $2/1_intermediate_files/varanalyzer_output.txt -gff $f0/$my_gff -pname $2  -cross $my_cross -snp_analysis_type $snp_analysis_type  2>> $my_log_file

	# (2) Run VA filter: eliminate SNPs with FA > 0.5 from control reads
	{
		python3 $location/scripts_snp/variants-filter.py -a $f1/control_filtered.va -b $f1/control_filtered2.va -step 3 -fasta $f1/$my_gs -af_max 0.5 2>> $my_log_file

	} || {
		echo $(date "+%F > %T")': Error during execution of variants-filter.py with control data.' >> $my_log_file
		exit_code=1
		echo $exit_code
		exit
	}
	echo $(date "+%F > %T")': First VCF filtering step of control data finished.' >> $my_log_file

	# (3) Run af-comparison: Intersection of filtered control SNPs with problem reads: outputs VA file with 4 columns of allele absolute frequence
	{
		python3 $location/scripts_snp/af-comparison.py -mode $my_mutbackgroud -f2_mut $f1/F2_filtered.va -f2_wt $f1/control_filtered2.va -out $f1/F2_control_comparison.va -f_input $f1/$my_gs -step 1 2>> $my_log_file 

	} || {
		echo $(date "+%F > %T")': Error during execution of af_comparison.py .' >> $my_log_file
		exit_code=1
		echo $exit_code
		exit
	}
	echo $(date "+%F > %T")': Allelic frequence comparison finished.' >> $my_log_file

	# (4) Run mapping analysis
	my_analysis_mode=back

	if [ $(wc -l < $f1/F2_control_comparison.va) -gt 1 ]
	then 
		{
			python3 $location/scripts_snp/map-mutation.py -file $f1/F2_control_comparison.va -fasta $f1/$my_gs -mode $my_analysis_mode -window_size 350000 -window_space 250000 -output $f1/map_info.txt -control_modality $my_mutbackgroud -interval_width $interval_width -snp_analysis_type par  2>> $my_log_file

		} || {
			echo $(date "+%F > %T")': Error during execution of map-mutation.py .' >> $my_log_file
			exit_code=1
			echo $exit_code
			exit
		}
		echo $(date "+%F > %T")': Mutation mapping module finished.' >> $my_log_file

		# (5) Re-write F2_control_comparison.va
		{
			python3 $location/scripts_snp/af-comparison.py -mode $my_mutbackgroud -f2_mut $f1/F2_filtered.va -f2_wt $f1/control_filtered2.va -out $f1/F2_control_comparison.va -f_input $f1/$my_gs -step 2 2>> $my_log_file 

		} || {
			echo $(date "+%F > %T")': Error during execution of af_comparison.py .' >> $my_log_file
			exit_code=1
			echo $exit_code
			exit
		}
		echo $(date "+%F > %T")': Allelic frequence comparison finished.' >> $my_log_file

		# Filler SNPs: AFs between 0.2 and 0.8 present in both samples, only for drawing 

		{
			python3 $location/scripts_snp/af-comparison.py -mode $my_mutbackgroud -f2_mut $f1/F2_filtered.va -f2_wt $f1/control_filtered.va -out $f1/filler_variants.va -f_input $f1/$my_gs -step 3 2>> $my_log_file 

		} || {
			echo $(date "+%F > %T")': Error during execution of af_comparison.py .' >> $my_log_file
			exit_code=1
			echo $exit_code
			exit
		}
		echo $(date "+%F > %T")': Allelic frequence comparison finished.' >> $my_log_file

		# (7) Candidate region analysis funtion
		cr_analysis F2_control_comparison.va 0.1

	else 
		echo $(date "+%F > %T")': No informative polymorphisms can be identified for mutation mapping. Please review your input data.' >> $my_log_file
		exit_code=1
		echo $exit_code
		exit
	fi	
fi


# _________________________________________________________________________________________________________________________________________________________________________________

# _________________________________Case 3: Mutant in ref background, outcross, mutant parental control (ref oc mut)________________________________________________________________
# _________________________________________________________________________________________________________________________________________________________________________________


if [ $my_mutbackgroud == ref ] && [ $my_pseq == mut ] && [ $my_cross == oc ]  && [ $snp_analysis_type == par ]
then

	# (1) Get problem and control VA files
	get_problem_va
	get_control_va

	#draw snps
	python3 $location/graphic_output/graphic-output.py -my_mut af_control -asnp $f1/control_filtered.va -bsnp $f1/$my_gs -rrl $my_rrl -iva $2/1_intermediate_files/varanalyzer_output.txt -gff $f0/$my_gff -pname $2  -cross $my_cross -snp_analysis_type $snp_analysis_type  2>> $my_log_file
	python3 $location/graphic_output/graphic-output.py -my_mut af_sample -asnp $f1/F2_filtered.va -bsnp $f1/$my_gs -rrl $my_rrl -iva $2/1_intermediate_files/varanalyzer_output.txt -gff $f0/$my_gff -pname $2  -cross $my_cross -snp_analysis_type $snp_analysis_type  2>> $my_log_file

	# (2) Run VA operations: Remove control SNPs from problem file
	my_operation_mode=A
	{
		python3 $location/scripts_snp/variants-operations.py -a $f1/F2_filtered.va -b $f1/control_filtered.va -c $f1/F2_control_comparison.va -mode $my_operation_mode -primary 1  2>> $my_log_file

	} || {
		echo $(date "+%F > %T")': Error during first execution of variants-operations.py .' >> $my_log_file
		exit_code=1
		echo $exit_code
		exit
	}
	echo $(date "+%F > %T")': VCF operations finished.' >> $my_log_file

	# (3) Run mapping analysis
	my_analysis_mode=out

	if [ $(wc -l < $f1/F2_control_comparison.va) -gt 1 ]
	then 
		{
			python3 $location/scripts_snp/map-mutation.py -file $f1/F2_control_comparison.va -fasta $f1/$my_gs -mode $my_analysis_mode -window_size 250000 -window_space 25000 -output $f1/map_info.txt -control_modality $my_mutbackgroud -interval_width $interval_width -snp_analysis_type $snp_analysis_type  2>> $my_log_file


		} || {
			echo $(date "+%F > %T")': Error during execution of map-mutation.py .' >> $my_log_file
			exit_code=1
			echo $exit_code
			exit
		}
		echo $(date "+%F > %T")': Mutation mapping module finished.' >> $my_log_file

		# (4) Candidate region analysis function
		cr_analysis F2_control_comparison.va 0.1

	else 
		echo $(date "+%F > %T")': No informative polymorphisms can be identified for mutation mapping. Please review your input data.' >> $my_log_file
		exit_code=1
		echo $exit_code
		exit
	fi
	
fi


#_________________________________________________________________________________________________________________________________________________________________________________

#_________________________________Case 4: Mutant in ref background, outcross, wt parental control (ref oc nomut)________________________________________________________________
#_________________________________________________________________________________________________________________________________________________________________________________


if [ $my_mutbackgroud == ref ] && [ $my_pseq == nomut ] && [ $my_cross == oc ]  && [ $snp_analysis_type == par ]
then
	# (1) Get control VA file
	get_control_va

	# (2) Run vcf filter to get SNPs with af > 0.75
	{
		python3 $location/scripts_snp/variants-filter.py -a $f1/control_filtered.va -b $f1/control_filtered2.va -step 3 -fasta $f1/$my_gs -af_min 0.75  2>> $my_log_file
		#draw snps
		python3 $location/graphic_output/graphic-output.py -my_mut af_control -asnp $f1/control_filtered2.va -bsnp $f1/$my_gs -rrl $my_rrl -iva $2/1_intermediate_files/varanalyzer_output.txt -gff $f0/$my_gff -pname $2  -cross $my_cross -snp_analysis_type $snp_analysis_type  2>> $my_log_file

	} || {
		echo $(date "+%F > %T")': Error during the second execution of variants-filter.py .' >> $my_log_file
		exit_code=1
		echo $exit_code
		exit
	}
	echo $(date "+%F > %T")': Second VCF filtering step finished.' >> $my_log_file


	# (3) Change ref seq, generate a "noref genome"
	{
		python3 $location/scripts_snp/change-snp.py -var $f1/control_filtered2.va  -gnm_ref $f1/$my_gs -out $f1/gnm_ref_merged/genome2.fa  2>> $my_log_file
	
		rm -rf $f1/$my_gs
		mv $f1/gnm_ref_merged/genome2.fa $f1/gnm_ref_merged/genome.fa

	} || {
		echo $(date "+%F > %T")': Error during first execution of variants-operations.py .' >> $my_log_file
		exit_code=1
		echo $exit_code
		exit
	}
	echo $(date "+%F > %T")': VCF operations finished.' >> $my_log_file

	# (4) Get problem VA file
	get_problem_va

	#draw snps
	python3 $location/graphic_output/graphic-output.py -my_mut af_sample -asnp $f1/F2_filtered.va -bsnp $f1/$my_gs -rrl $my_rrl -iva $2/1_intermediate_files/varanalyzer_output.txt -gff $f0/$my_gff -pname $2  -cross $my_cross -snp_analysis_type $snp_analysis_type  2>> $my_log_file

	# (5) Run VA operations: Intersection to get SNPs for mapping the mutation
	my_operation_mode=I
	{
		python3 $location/scripts_snp/variants-operations.py -a $f1/F2_filtered.va -b $f1/control_filtered2.va -c $f1/F2_control_comparison_mapping.va -mode $my_operation_mode -primary 1  2>> $my_log_file

	} || {
		echo $(date "+%F > %T")': Error during first execution of variants-operations.py .' >> $my_log_file
		exit_code=1
		echo $exit_code
		exit
	}
	echo $(date "+%F > %T")': VCF operations finished.' >> $my_log_file

	# (6) Run mapping analysis
	my_analysis_mode=out

	if [ $(wc -l < $f1/F2_control_comparison_mapping.va) -gt 1 ]
	then 
		{
			python3 $location/scripts_snp/map-mutation.py -file $f1/F2_control_comparison_mapping.va -fasta $f1/$my_gs -mode $my_analysis_mode -window_size 250000 -window_space 25000 -output $f1/map_info.txt -control_modality noref -interval_width $interval_width -snp_analysis_type $snp_analysis_type  2>> $my_log_file

		} || {
			echo $(date "+%F > %T")': Error during execution of map-mutation.py .' >> $my_log_file
			exit_code=1
			echo $exit_code
			exit
		}
		echo $(date "+%F > %T")': Mutation mapping module finished.' >> $my_log_file

		# (7) Run VA operations: Remove control SNPs from problem file 
		my_operation_mode=A
		{
			python3 $location/scripts_snp/variants-operations.py -a $f1/F2_filtered.va -b $f1/control_filtered.va -c $f1/F2_control_comparison.va -mode $my_operation_mode -primary 1  2>> $my_log_file
			#draw snps
			#python3 $location/graphic_output/graphic-output.py -my_mut af_candidates -asnp $f1/F2_control_comparison.va -bsnp $f1/$my_gs -rrl $my_rrl -iva $2/1_intermediate_files/varanalyzer_output.txt -gff $f0/$my_gff -pname $2  -cross $my_cross -snp_analysis_type $snp_analysis_type  

		} || {
			echo $(date "+%F > %T")': Error during first execution of variants-operations.py .' >> $my_log_file
			exit_code=1
			echo $exit_code
			exit
		}
		echo $(date "+%F > %T")': VCF operations finished.' >> $my_log_file

		# (8) Candidate region analysis function
		cr_analysis F2_control_comparison_mapping.va 0.1
	else 
		echo $(date "+%F > %T")': No informative polymorphisms can be identified for mutation mapping. Please review your input data.' >> $my_log_file
		exit_code=1
		echo $exit_code
		exit
	fi
fi



#_________________________________________________________________________________________________________________________________________________________________________________

#_________________________________Case 7: Mutant in noref background, outcross, mutant parental control (noref oc mut)________________________________________________________________
#_________________________________________________________________________________________________________________________________________________________________________________

if [ $my_mutbackgroud == noref ] && [ $my_pseq == mut ] && [ $my_cross == oc ]  && [ $snp_analysis_type == par ]
then

	# (1) Get problem and control VA files
	get_problem_va
	get_control_va

	#draw snps
	python3 $location/graphic_output/graphic-output.py -my_mut af_control -asnp $f1/control_filtered.va -bsnp $f1/$my_gs -rrl $my_rrl -iva $2/1_intermediate_files/varanalyzer_output.txt -gff $f0/$my_gff -pname $2  -cross $my_cross -snp_analysis_type $snp_analysis_type  2>> $my_log_file
	python3 $location/graphic_output/graphic-output.py -my_mut af_sample -asnp $f1/F2_filtered.va -bsnp $f1/$my_gs -rrl $my_rrl -iva $2/1_intermediate_files/varanalyzer_output.txt -gff $f0/$my_gff -pname $2  -cross $my_cross -snp_analysis_type $snp_analysis_type  2>> $my_log_file

	# (2) Run VA operations: Intersection to get mapping SNPs
	my_operation_mode=I
	{
		python3 $location/scripts_snp/variants-operations.py -a $f1/F2_filtered.va -b $f1/control_filtered.va -c $f1/F2_control_comparison_mapping.va -mode $my_operation_mode -primary 1  2>> $my_log_file

	} || {
		echo $(date "+%F > %T")': Error during first execution of variants-operations.py .' >> $my_log_file
		exit_code=1
		echo $exit_code
		exit
	}
	echo $(date "+%F > %T")': VCF operations finished.' >> $my_log_file

	# (3) Run mapping analysis
	my_analysis_mode=out


	if [ $(wc -l < $f1/F2_control_comparison_mapping.va) -gt 1 ]
	then 
		{
			python3 $location/scripts_snp/map-mutation.py -file $f1/F2_control_comparison_mapping.va -fasta $f1/$my_gs -mode $my_analysis_mode -window_size 250000 -window_space 25000 -output $f1/map_info.txt -control_modality $my_mutbackgroud -interval_width $interval_width -snp_analysis_type $snp_analysis_type  2>> $my_log_file

		} || {
			echo $(date "+%F > %T")': Error during execution of map-mutation.py .' >> $my_log_file
			exit_code=1
			echo $exit_code
			exit
		}
		echo $(date "+%F > %T")': Mutation mapping module finished.' >> $my_log_file


		# (4) Run VA operations: Remove control SNPs from problem
		my_operation_mode=A
		{
			python3 $location/scripts_snp/variants-operations.py -a $f1/F2_filtered.va -b $f1/control_filtered.va -c $f1/F2_control_comparison.va -mode $my_operation_mode -primary 1  2>> $my_log_file
			#draw snps
			#python3 $location/graphic_output/graphic-output.py -my_mut af_candidates -asnp $f1/F2_control_comparison.va -bsnp $f1/$my_gs -rrl $my_rrl -iva $2/1_intermediate_files/varanalyzer_output.txt -gff $f0/$my_gff -pname $2  -cross $my_cross -snp_analysis_type $snp_analysis_type  

		} || {
			echo $(date "+%F > %T")': Error during first execution of variants-operations.py .' >> $my_log_file
			exit_code=1
			echo $exit_code
			exit
		}
		echo $(date "+%F > %T")': VCF operations finished.' >> $my_log_file

		# (5) Candidate region analysis function
		cr_analysis F2_control_comparison_mapping.va 0.1

	else 
		echo $(date "+%F > %T")': No informative polymorphisms can be identified for mutation mapping. Please review your input data.' >> $my_log_file
		exit_code=1
		echo $exit_code
		exit
	fi
fi


#Intermediate files cleanup
rm -f $f1/*.bam
rm -f $f1/*.bai
rm -f $f1/*.ht2
rm -rf $f1/gnm_ref_merged
if [ -d "$f1/sim_data" ]; then rm -Rf $f1/sim_data/; fi


echo $exit_code
