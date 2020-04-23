#!./src/Python-2.7.18/.localpython/bin/python2
import cgi, cgitb, os, subprocess
cgitb.enable() # For debugging only

# Set Content-type header so XMLHttpRequest in JS can understand the response
print "Content-Type: text/html"
print ""

arguments = cgi.FieldStorage()
projectName = str(arguments['p'].value).strip()

# Default PID values
pid_easymap = 0
pid_simulator = 0
pid_workflow = 0

# Fill PID values
with open('./user_projects/'+projectName+'/2_logs/status') as status_file:
	for line in status_file.readlines():
		if line.startswith('pid easymap'): 
			pid_easymap = str(line.split('easymap')[1]).strip()
			
		if line.startswith('pid simulator'): 
			pid_simulator = str(line.split('simulator')[1]).strip()
		
		if line.startswith('pid workflow'): 
			pid_workflow = str(line.split('workflow')[1]).strip()
		

simulator_children = [
	'sim-mut.py',
	'sim-recsel.py',
	'sim-seq.py'
	'calculate-genome-length.py'
	]

workflow_children = [
	'bowtie2-build-s',
	'bowtie2-align-s',
	'samtools',
	'bcftools',
	'draw.py',
	'draw.pyc',
	'graphic-output.py',
	'local-analysis.py',
	'paired-analysis.py',
	'filter1.py',
	'filter2.py',
	'sam-file-check.py',
	'sort.py',
	'ins-to-varanalyzer.py',
	'af-comparison.py',
	'map-mutation.py',
	'variants-filter.py',
	'vcf-groomer.py',
	'variants-operations.py',
	'snp-to-varanalyzer.py',
	'varanalyzer.py'
	'ins-primers.py'
	'primer-generation.py'
	'extend-ins-info.py'
	'depth_measures_generation.py'
	'graphic-alignment.py'
	'report.py'
	'snp-to-varanalyzer.py'
	'extend-snp-variants-info.py'
	'af-comparison.py'
	'change-snp.py'
	]


# Kill all children
if pid_simulator != 0: 
	for child in simulator_children:
		command = 'pkill -9 -P ' + pid_simulator + ' -f ' + child
		subprocess.call(command, shell=True)

if pid_workflow != 0: 
	for child in workflow_children:
		command = 'pkill -9 -P ' + pid_workflow + ' -f ' + child
		subprocess.call(command, shell=True)



# Kill .sh parents
if pid_workflow != 0: 
	command = 'pkill -9 -P ' + pid_workflow
	subprocess.call(command, shell=True)

if pid_simulator != 0: 
	command = 'pkill -9 -P ' + pid_simulator
	subprocess.call(command, shell=True)

if pid_easymap != 0: 
	command = 'pkill -9 -P ' + pid_easymap
	subprocess.call(command, shell=True)

# Append 'status: killed' to status file
with open('./user_projects/'+projectName+'/2_logs/status', 'a') as status_file:
	status_file.write('status:killed')

# Write to log
log_file = './user_projects/'+projectName+'/2_logs/log.log'
subprocess.call('echo $(date)": Project interrupted by the user." >> ' + log_file, shell=True)


































'''
<?php

/*
I am unable to kill all children processes by using the -TERM argument inside a
shell_exec(). I use a longer alternative instead:

Make that process-input.sh, simulator.sh, and the worklflows write their PIDs to file
status in the logs folder. Then with PHP I get the PIDs and use something like:

kill -9 <PID>                  (for easymap.sh)
pkill -9 -P <PPID> xxx.sh      (for the workflows) 
pkill -9 -P <PPID> program     (for bowtie, samtools, bcftools, python scripts...)

Decide whether to use signal '9' or not in the command kill
*/

// Create variables for the PIDs and fill the with default value 0
$pid_easymap = 0;
$pid_simulator = 0;
$pid_workflow = 0;

// Get the PIDs of the .sh scripts, which are stored in /2_logs/status
$project = $_GET['p'];

$status_file = '../user_projects/'. $project .'/2_logs/status';
//$status_contents = file_get_contents($status_file);
$status_contents = fopen($status_file, 'r');
while(!feof($status_contents)) {
	$line = fgets($status_contents);
	$line_fields = explode(' ', $line);
	if ($line_fields[0] == 'pid' AND $line_fields[1] == 'easymap') {
		$pid_easymap = $line_fields[2];
	}
	if ($line_fields[0] == 'pid' AND $line_fields[1] == 'simulator') {
		$pid_simulator = $line_fields[2];
	}
	if ($line_fields[0] == 'pid' AND $line_fields[1] == 'workflow') {
		$pid_workflow = $line_fields[2];
	}
}
fclose($status_contents);


//echo $pid_easymap .'<br>';
//echo $pid_simulator .'<br>';
//echo $pid_workflow .'<br>';

// Store in arrays all the programs/scripts that are direct children of simulator.sh and 
// of workflow-ins/snp.sh
$simulator_children = array(
	'sim-mut.py',
	'sim-recsel.py',
	'sim-seq.py'
);

$workflow_children = array(
	'bowtie2-build-s',
	'bowtie2-align-s',
	'samtools',
	'bcftools',
	'draw.py',
	'draw.pyc',
	'graphic-output.py',
	'local-analysis.py',
	'paired-analysis.py',
	'filter1.py',
	'filter2.py',
	'sam-file-check.py',
	'sort.py',
	'ins-to-varanalyzer.py',
	'af-comparison.py',
	'map-mutation.py',
	'variants-filter.py',
	'vcf-groomer.py',
	'variants-operations.py',
	'snp-to-varanalyzer.py',
	'varanalyzer_v0.py',
	'varanalyzer_v1.py'
);

// Kill all processes that are direct children of simulator.sh and of workflow-ins/snp.sh
// I don't kill children of process-input.sh because runtime and output are negligible
foreach ($simulator_children as $simulator_element) {
	//echo 'pkill -9 -P '. $pid_simulator .' '. $simulator_element .'<br>';
	shell_exec('pkill -9 -P '. $pid_simulator .' '. $simulator_element);
}

foreach ($workflow_children as $workflow_element) {
	//echo 'pkill -9 -P '. $pid_workflow .' '. $workflow_element .'<br>';
	shell_exec('pkill -9 -P '. $pid_workflow .' '. $workflow_element);
}

// Kill .sh scripts using their PPID
shell_exec('pkill -9 -P '. $pid_easymap .' simulator.sh');
//echo 'pkill -9 -P '. $pid_easymap .' simulator.sh<br>';
shell_exec('pkill -9 -P '. $pid_easymap .' workflow-ins.sh');
//echo 'pkill -9 -P '. $pid_easymap .' workflow-ins.sh<br>';
shell_exec('pkill -9 -P '. $pid_easymap .' workflow-snp.sh');
//echo 'pkill -9 -P '. $pid_easymap .' workflow-snp.sh<br>';

// Kill easymap.sh
shell_exec('pkill -9 '. $pid_easymap);
//echo 'pkill -9 '. $pid_easymap .'<br>';
'''
'''
// Append 'status:killed' to status file
$status = fopen($status_file, 'a');
fwrite($status, 'status:killed');
fclose($status);

// Write to log
$log_file = '../user_projects/'. $project .'/2_logs/log.log';
shell_exec('echo $(date)": Project interrupted by the user." >> '. $log_file);

?>
'''