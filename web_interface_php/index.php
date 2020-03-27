<!DOCTYPE html>
<html>
<title>Easymap</title>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<link rel="stylesheet" href="w3c.css">
<link rel="stylesheet" href="style.css">
<link rel="stylesheet" href="https://fonts.googleapis.com/css?family=Poppins">
<script type="text/javascript" src="plupload-2.3.1/js/plupload.full.min.js"></script>
<!-- <script type="text/javascript" src="upload-files.js"></script> -->
<script type="text/javascript" src="manage-input-files.js"></script>
<style>
body,h1,h2,h3,h4,h5 {font-family: "Poppins", sans-serif}
body {font-size:16px;}
.w3-half img{margin-bottom:-6px;margin-top:16px;opacity:0.8;cursor:pointer}
.w3-half img:hover{opacity:1}
</style>
<body>

<!-- Sidebar/menu -->
<nav class="w3-sidebar w3-red w3-collapse w3-top w3-large w3-padding" style="z-index:3;width:300px;font-weight:bold;" id="mySidebar"><br>
  <a href="javascript:void(0)" onclick="w3_close()" class="w3-button w3-hide-large w3-display-topleft" style="width:100%;font-size:22px">Close Menu</a>
  <div class="w3-container">
    <h3 class="w3-padding-64"><b>Easymap</b></h3>
  </div>
  <div class="w3-bar-block">
    <a href="manage-input-files.php" onclick="w3_close()" class="w3-bar-item w3-button w3-hover-white">Manage input files</a> 
    <a href="manage-projects.php" onclick="w3_close()" class="w3-bar-item w3-button w3-hover-white">Manage projects</a> 
    <a href="run-new-project.php" onclick="w3_close()" class="w3-bar-item w3-button w3-hover-white">Run new project</a> 
    <a href="documentation.php" onclick="w3_close()" class="w3-bar-item w3-button w3-hover-white">Documentation</a>
  </div>
</nav>

<!-- Top menu on small screens -->
<header class="w3-container w3-top w3-hide-large w3-red w3-xlarge w3-padding">
  <a href="javascript:void(0)" class="w3-button w3-red w3-margin-right" onclick="w3_open()">☰</a>
  <span>Easymap</span>
</header>

<!-- Overlay effect when opening sidebar on small screens -->
<div class="w3-overlay w3-hide-large" onclick="w3_close()" style="cursor:pointer" title="close side menu" id="myOverlay"></div>

<!-- beginning of page content -->
<div class="w3-main" style="margin-left:340px;margin-right:40px">

  <div class="w3-container" style="margin-top:75px">
    <h1 class="w3-xxxlarge w3-text-red"><b>Access to easymap</b></h1>
    <hr style="width:50px;border:5px solid red" class="w3-round">
    <h3>Password:</h3>
     
		<?php

		// Get config/config>web-interface-password
		// (set by the system administrator)

		$config_contents = fopen('../config/config', 'r');

		while (!feof($config_contents)) {
			$line = fgets($config_contents);
			$fields = explode(':', $line);

			if ($fields[0] == 'web-interface-password') {
				$reference_password = trim($fields[1]);
			}
		}

		fclose($config_contents);

		// Start a session and store the reference password in $_SESSION[]
		session_start();
		$_SESSION['reference_password'] = $reference_password;
		?> 
      
    <form action="manage-input-files.php" method="post">
      <input type="password" name="password" style="margin:8px 10px 0 0; padding:8px; display:block; float:left;">
      <input type="submit" value="Submit" name="submit" class="button">
    </form>

  </div>
  
<!-- End of page content -->
</div>

<!-- W3.CSS Container -->
<div class="w3-light-grey w3-container w3-padding-32" style="margin-top:75px;padding-right:58px"><p class="w3-right">Easymap</a></p></div>

<script>
// Script to open and close sidebar
function w3_open() {
    document.getElementById("mySidebar").style.display = "block";
    document.getElementById("myOverlay").style.display = "block";
}
 
function w3_close() {
    document.getElementById("mySidebar").style.display = "none";
    document.getElementById("myOverlay").style.display = "none";
}
/*
// Modal Image Gallery
function onClick(element) {
  document.getElementById("img01").src = element.src;
  document.getElementById("modal01").style.display = "block";
  var captionText = document.getElementById("caption");
  captionText.innerHTML = element.alt;
}
*/
</script>

</body>
</html>

