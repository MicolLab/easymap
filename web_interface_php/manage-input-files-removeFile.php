<?php

$file = $_GET['f'];

$command = 'rm ../user_data/'. $file;

shell_exec($command);

?>
