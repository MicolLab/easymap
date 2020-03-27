<?php

// Run allow-new-project.py, get output, and send it to whichever script that calls the present script
$out  = shell_exec('cd ../config; python allow-new-project.py');
echo $out;


// The following commented code is the former way we used to check if data and compute limits
// specified in config/config had been reached

/*
// Get amount of data currently in user_projects directory
$folder_size_output = shell_exec('du -s ../user_projects');
$folder_size_output_array = explode('	', $folder_size_output);
$user_projects_size = $folder_size_output_array[0];

// Get number of projects currently running
$projects = array_slice(scandir('../user_projects'), 2);
$num_running_projects = 0;

foreach ($projects as $project) {
	$status_file = '../user_projects/'. $project .'/2_logs/status';		
	$status_contents = fopen($status_file, 'r');
	
	while(!feof($status_contents)) {
		$line = fgets($status_contents);
		$line_fields = explode(':', $line);
		
		if ($line_fields[0] == 'status') {
			$current_status = trim($line_fields[1]);
		}
	}
	
	if ($current_status == 'running') {
		$num_running_projects ++;
	}
	
	fclose($status_contents);
}


// Get config/config>user_projects-size-limit
// Get config/config>max-simultaneous-jobs
// (set by the system administrator)
// Compute % of allowed size in user_projects that user has used
// Compute if max simult. jobs has been reached

$config_contents = fopen('../config/config', 'r');

while (!feof($config_contents)) {
	$line = fgets($config_contents);
	$fields = explode(':', $line);
	
	if ($fields[0] == 'user_projects-size-limit') {
		$size_limit = trim($fields[1]);
		// If config files has "user_projects-size-limit:0", interpret it as unlimited
		if ($size_limit == 0) {
			$output1 = 1; // This simply tells javascript that user is using 1% (<100%) of the max allowed space
		} else {
			$output1 = floor(($user_projects_size / $size_limit) * 100);
		}
	}
	
	if ($fields[0] == 'max-simultaneous-jobs') {
		$max_jobs = trim($fields[1]);
		// If config files has "max-simultaneous-jobs:0", interpret it as unlimited
		if ($max_jobs == 0) {
			$output2 = 1; // This simply tell javascript that user is running 1% (<100%) of the max allowed simultaneous jobs
		} else {
			$output2 = floor(($num_running_projects / $max_jobs) * 100);
		}
	}
} 

fclose($config_contents);

if (!isset($output1)) {
	$output1 = 'user_projects-size-limit is not properly configured in config/config file';
}
if (!isset($output2)) {
	$output2 = 'max-simultaneous-jobs is not properly configured in config/config file';
}

echo $output1 .','. $output2 .','. $size_limit .','. $max_jobs;

*/

?>
