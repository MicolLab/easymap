<?php

$file = $_GET['f'];

$nbrOfLines = 1000;

$command = 'head -'. $nbrOfLines .' ../user_data/'. $file;

$preview = shell_exec($command);

echo '
	<p>File: '. $file .'.</p>
	<p>(showing first '. $nbrOfLines .' lines)</p>
	<p><pre>***********************************************************</pre></p>
	<pre>'. htmlentities($preview) .'</pre>
';

?>
