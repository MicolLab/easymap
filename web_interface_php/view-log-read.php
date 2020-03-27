<?php
$project_name = $_GET['p'];
$location = '../user_projects/'. $project_name .'/2_logs/log.log';
$log_contents = fopen($location, "r");
while (!feof($log_contents)) {
  $line = fgets($log_contents);
  echo $line .'<br>';
}
fclose($log_contents);
?>