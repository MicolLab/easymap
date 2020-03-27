<?php

$projects = array_reverse(array_slice(scandir('../user_projects'), 2));

if (empty($projects)) {
	echo '<div class="emptyResult">There are no projects currently on the server. To start a new one go to "Run new project".</div>';
} else {
	foreach ($projects as $project) {	
		
		echo '<div style="background-color:#79e59d; border:solid green 1px; border-radius:5px; padding:10px;">';

		$status_file = '../user_projects/'. $project .'/2_logs/status';
		if (file_exists($status_file)) {
			$status_contents = fopen($status_file, 'r');
			while(!feof($status_contents)) {
				$line = fgets($status_contents);
				
				$line_fields = explode(':', $line);
				if ($line_fields[0] == 'status') {
					$current_status = trim($line_fields[1]);
				}
			}
			fclose($status_contents);
			
			// Project name and status
			echo '<h4>Project: '. $project .'</h4>';
			echo 'Project status: '. $current_status .'<br>';
			
			// Folder size
			if ($current_status != 'running') {
				$folder_size_output = shell_exec('du -sh ../user_projects/'. $project);
				$folder_size_output_array = explode('	', $folder_size_output);
				$folder_size = $folder_size_output_array[0];
				echo 'Project size: '. $folder_size .'<br>';
			}
			
			echo '<div class="managing-buttons">';
				
			// Log link
			echo '<a href="view-log.php?p='. $project .'" class="button">View log file</a>';
			
			if ($current_status == 'finished') {
				echo '<a href="view-report.php?p='. $project .'" class="button">View report</a>';
			}
			
			if ($current_status == 'running') {
				echo '<form><input type="button" class="button" onclick="displayStopProject(\''. $project .'\')" value="Stop execution" /></form>';
			}
			
			if ($current_status != 'running') {
				echo '<form><input type="button" class="button" onclick="displayRemoveProject(\''. $project .'\')" value="Remove from disk" /></form>';
			}
			
			echo '</div>';

			echo '<div id="stopPrj_'. $project .'" style="display:none;">
				Are you sure you want to stop running this project?
				<form>
					<input type="button" onclick="stopProject(\''. $project .'\')" value="Confirm"/>
					<input type="button" onclick="hideStopProject(\''. $project .'\')" value="Cancel"/>
				</form>
			</div>';

			echo '<div id="rmPrj_'. $project .'" style="display:none;">
				Are you sure you want to permanently remove this project from disk? (Your input data will not be removed).
				<form>
					<input type="button" onclick="removeProject(\''. $project .'\')" value="Confirm"/>
					<input type="button" onclick="hideRemoveProject(\''. $project .'\')" value="Cancel"/>
				</form>
			</div>';

		} else {
			// Project name and status
			echo '<h4>Project: '. $project .'</h4>';
			echo 'Project status: The status file of this project could not be found on the server. Please delete this project and rerun it.<br>';
			echo '<div class="managing-buttons">';
			echo '<form><input type="button" class="button" onclick="displayRemoveProject(\''. $project .'\')" value="Remove from disk" /></form>';
			echo '</div>';

			echo '<div id="rmPrj_'. $project .'" style="display:none;">
				Are you sure you want to permanently remove this project from disk? (Your input data will not be removed).
				<form>
					<input type="button" onclick="removeProject(\''. $project .'\')" value="Confirm"/>
					<input type="button" onclick="hideRemoveProject(\''. $project .'\')" value="Cancel"/>
				</form>
			</div>';
		}

		echo '</div><br>';
	}
}


//echo $project .' --- Status --- View log file --- View report --- Stop execution --- Remove from disk<br>';
//                     running    +                 -               +                  + > -                
//                     finished   +                 +               -                  + > -
//                     killed     +                 -               -                  + > -
//                     error      +                 -               -                  + > -





?>
