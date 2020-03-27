<?php

/*
This file creates the command with info from javascript (user interface) and tells the shell to run
easymap.sh (the master .sh program to run easymap workflows)

Warning: If easymap.sh is run from php, the user is www-data (the apache server). Therefore www-data
must have permission to read and write in the appropriate folders.
*/

// Send 200 OK to close the client request before easymap execution finishes
// (https://stackoverflow.com/questions/15273570/continue-processing-php-after-sending-http-response)

//ignore_user_abort(true);
//ob_start();
// do initial processing here
//echo $response; // send the response
header("HTTP/1.1 200 OK");
//header('Connection: close');
//header('Content-Length: '.ob_get_length());
ob_end_flush();
ob_flush();
flush();

// Now handle user request

// Handle client data in JSON format
//header("Content-Type: application/json");

// build a PHP variable from JSON sent using POST method
//$cmdArray = json_decode(stripslashes(file_get_contents("php://input")));
$cmdArray = json_decode(file_get_contents("php://input"));

// Elaborate the command string
//$refSeqsString = implode("+", $cmdArray[4]);
//$cmdArray[4] = $refSeqsString;
$cmdString = implode(" ", $cmdArray);

// Run workflow
// Add another argument that equals to "server" at the end of the command
shell_exec('cd ..; '. $cmdString .' server');

// Only for development
//$file = fopen('file.txt', 'w');
//fwrite($file, $cmdString);
//fclose($file);


?>
