<?php

$project = $_GET['p'];

$command = 'rm -rf --recursive ../user_projects/'. $project;

shell_exec($command);

?>
