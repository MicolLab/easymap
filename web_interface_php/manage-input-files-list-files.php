<?php

// Get in an array all items in 'user_data' directory excluding references
// to self (.) and parent (..) dirs
$user_input_files = array_slice(scandir('../user_data'), 2);

// DEPRECATED
// Remove folder 'gnm-ref' from array
// This folder is only useful for the command-line version
//foreach ($user_input_files as $key => $item) {
//	if (trim($item) == 'gnm_ref') {
//		unset($user_input_files[$key]);
//	}
//}


if (empty($user_input_files)) {
	echo '<div class="emptyResult">There are no input files currently on the server.</div>';
} else {
	echo '<hr class="files-separator"></hr>';

	foreach ($user_input_files as $user_input_file) {
		
		// Get file size to display it to user
		$file_size_output = shell_exec('du -sh  ../user_data/'. $user_input_file);
		$file_size_output_array = explode('	', $file_size_output);
		$file_size = $file_size_output_array[0];
		
		echo '
			<div id="file_'. $user_input_file .'">
				<div class="files-item left">
					<h4>'. $user_input_file .'</h4>
				</div>
				<div class="files-item center">
					<h4>'. $file_size .'</h4>
				</div>
				<div class="files-item right1">
					<a href="manage-input-files-retrieve-header.php?f='. $user_input_file .'" target="_blank" class="button" style="width: 100px">Preview</a>
				</div>
				<div class="files-item right2">
					<form>
						<input style="margin: 3 0 7px 0" type="button" class="button" onclick="ShowRemoveFile(\''. $user_input_file .'\')" value="Remove from disk" />
					</form>
				</div>
				<div id="rmFile_'. $user_input_file .'" style="display:none;">
					Are you sure you want to permanently remove this file from disk?
					<form>
						<input type="button" onclick="removeFile(\''. $user_input_file .'\')" value="Confirm" />
						<input type="button" onclick="HideRemoveFile(\''. $user_input_file .'\')" value="Cancel" />
					</form>
				</div>
				<hr class="files-separator"></hr>
			</div>
		';
	}
}

?>
