<?php
// Start of code for dealing with password

session_start();
if (!isset($_SESSION['reference_password'])) {
	header("Location: index.php");
	exit();
}

$reference_password = $_SESSION['reference_password'];
if (isset($_POST['password'])) {
	$_SESSION['password'] = $_POST['password'];
	$password = $_SESSION['password'];
} else {
	@$password = $_SESSION['password'];
}

if ($password !== $reference_password) {

	echo '
		<!DOCTYPE html>
		<html>
			<head>
				<title></title>
			</head>
			<body>
				Wrong password. <a href="index.php">Try again</a>
			</body>
		</html>
	';
}

if ($password === $reference_password) {
// End of code for dealing with password
?>

<!DOCTYPE html>
<html>
<title>Easymap</title>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1">

<link rel="stylesheet" href="w3c.css">
<link rel="stylesheet" href="style.css">
<link rel="stylesheet" href="https://fonts.googleapis.com/css?family=Poppins">

<script type="text/javascript" src="run-new-project.js"></script>

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
    <h1 class="w3-xxxlarge w3-text-red"><b>Run new project</b></h1>
    <hr style="width:50px;border:5px solid red" class="w3-round">
    <h3>Design and execute a new project</h3>
    
    <!-- If the size of data or the number or currently running projects is over the limits established by machine administrator,
    display Warning messages and do not display div with id runNewProject -->
    
    <div id="sizeWarning" class="warningMessage"></div>
    <div id="simultWarning" class="warningMessage"></div>
    
    <div id="runNewProject">
			
		<form id="form1">
			
			<div class="formInfo" id="formInfoSim">General parameters</div>
			
			<p>
				Project name (only alphanumeric characters are allowed):<br>
				<input type="text" name="projectName" value="My project" />
				<div id="projectNameValidationInfo" class="warningMessage"></div>
			</p>
			
			Mapping-by-sequencing strategy:
			<div class="buttons-wrap" id="analysisTypeWrap">
				<div class="mx-button">
					<input type="radio" class="analysisType" name="mx12" id="button1" />
					<label for="button1" unselectable>Tagged sequence mapping</label>
				</div>
				<div class="mx-button">
					<input type="radio" class="analysisType" name="mx12" id="button2" />
					<label for="button2" unselectable>Linkage analysis mapping</label>
				</div>
				<div class="clear-floats"></div>
			</div>

			<div id="analysisTypeValidationInfo" class="warningMessage"></div>
			
			Data source:
			<div class="buttons-wrap">
				<div class="mx-button">
					<input type="radio" class="dataSource" name="mx34" id="button3" />
					<label for="button3" unselectable>Use my own reads</label>
				</div>
				<div class="mx-button">
					<input type="radio" class="dataSource" name="mx34" id="button4" />
					<label for="button4" unselectable>Simulate data</label>
				</div>
				<div class="clear-floats"></div>
			</div>

			<div id="dataSourceValidationInfo" class="warningMessage"></div>

			<hr style="width: 100%; border: 2px solid rgb(150,150,150)" class="w3-round">
			
			Reference sequence:<br>
			<select id="refFileSelector" style="display:block; margin-bottom:22px;"></select>
			<div id="refSeqValidationInfo" class="warningMessage"></div>
			
			<div id="insSeqField"> <!-- Not displayed by default -->
				Insertion sequence file:<br>
				<select id="insFileSelector" style="display:block; margin-bottom:22px;"></select>
				<div id="insFileValidationInfo" class="warningMessage"></div>
			</div>
			
			GFF3 file (gene structural annotation):<br>
			<select id="gffFileSelector" style="display:block; margin-bottom:22px;"></select>
			<div id="gffFileValidationInfo" class="warningMessage"></div>

			Gene functional annotation file [OPTIONAL]:<br>
			<select id="annFileSelector" style="display:block; margin-bottom:22px;"></select>
			<div id="annFileValidationInfo" class="warningMessage"></div>
			
			<div id="backgroundCrossCtype" style="display:none;"> <!-- Not displayed by default -->
				<hr style="width: 100%; border: 2px solid rgb(150,150,150)" class="w3-round">
				
				Mutant background:
				<div class="buttons-wrap">
					<div class="mx-button">
						<input type="radio" class="mutBackground" name="mxMB" id="button11" />
						<label for="button11" unselectable>Reference</label>
					</div>
					<div class="mx-button">
						<input type="radio" class="mutBackground" name="mxMB" id="button12" />
						<label for="button12" unselectable>Non-reference</label>
					</div>
					<div class="clear-floats"></div>
				</div>

				<div id="mutBackgroundValidationInfo" class="warningMessage"></div>

				Mapping cross performed:
				<div class="buttons-wrap">
					<div class="mx-button">
						<input type="radio" class="crossType" name="mxCR" id="button13" />
						<label for="button13" unselectable>Backcross</label>
					</div>
					<div class="mx-button">
						<input type="radio" class="crossType" name="mxCR" id="button14" />
						<label for="button14" unselectable>Outcross</label>
					</div>
					<div class="clear-floats"></div>
				</div>

				<div id="crossTypeValidationInfo" class="warningMessage"></div>

				Origin of the control reads:
				<div class="buttons-wrap">
					<div class="mx-button">
						<input type="radio" class="contType" name="mxCO" id="button15" />
						<label for="button15" unselectable>Mutant parental</label>
					</div>
					<div class="mx-button">
						<input type="radio" class="contType" name="mxCO" id="button16" />
						<label for="button16" unselectable>Polymorphic parental</label>
					</div>
					<div class="mx-button">
						<input type="radio" class="contType" name="mxCO" id="button17" />
						<label for="button17" unselectable>F2 wild-types</label>
					</div>
					<div class="clear-floats"></div>
				</div>

				<div id="contTypeValidationInfo" class="warningMessage"></div>

				<div id="backgroundCrossCtypeWarnMsg" class="warningMessage">Invalid combination. Easymap does not support this experimental design.</div>
			
				<hr style="width: 100%; border: 2px solid rgb(150,150,150)" class="w3-round">

				Use low stringency during SNP analyses?

				<div id="checkBoxStringency">
					<input type="checkbox" class="switch_1" id="stringency">
				</div>

				By default, Easymap only considers SNPs that pass certain quality checks. However, in some read datasets that are not optimal, 
				performing more lenient filtering can help to identify a candidate interval and the causal mutation. If you analyze your reads
				in the default mode and obtain very few SNPs or believe that the causal mutation could have been discarded, turn on this option
				and run the program again.

			</div>

			<div id="expDataInterface">
				<div class="formInfo" id="formInfoSim">Experimental reads</div>

				Problem reads (if your reads are paired-end, select both files while holding the Ctrl/Cmd key):<br>
				<select multiple id="readsProblemSelector" style="display:block; margin-bottom:22px;"></select>
				<div id="readsProblemWarnMsg" class="warningMessage">Please select one file for single-end reads and two files for paired-end reads.</div>
				<div id="readsProblemWarnMsg2" class="warningMessage">The problem and control read files must be different.</div>
				
				<div id="readsControl">
					Control reads (if your reads are paired-end, select both files while holding the Ctrl/Cmd key):<br>
					<select multiple id="readsControlSelector" style="display:block; margin-bottom:22px;"></select>
					<div id="readsControlWarnMsg" class="warningMessage">Please select one file for single-end reads and two files for paired-end reads.</div>
					<div id="readsControlWarnMsg2" class="warningMessage">The problem and control read files must be different.</div>
				</div>
			</div>
		
			<div id="simDataInterface">
				<div class="formInfo" id="formInfoSim">Simulation parameters</div>

				<div id="simMutInterface">
					Mutagenesis parameters:<br>
					<input type="text" name="simMut" class="sim" style="display:block;" value='{"numberMutations":"0"}' />
					<div id="simMutValMsg" class="warningMessage"></div>
				</div>
				
				<div id="simRecselInterface">
					Recombination and selection parameters:<br>
					<input type="text" name="simRecselA" class="sim" style="display:block;" value='{"contigCausalMut":"1", "posCausalMut":"1000000", "numRecChrs":"100"}' />
					<div id="simRecselAValMsg" class="warningMessage"></div>
					Contig recombination frequency distributions:<br>
					<input type="text" name="simRecselB" class="sim" style="display:block;" value='{"1":"0,14-1,31-2,33-3,15-4,5-5,2", "2":"0,24-1,42-2,25-3,6-4,1-5,2"}' />
					<div id="simReqselBValMsg" class="warningMessage"></div>
				</div>
		
				<div id="simSeqInterface">
					Sequencing parameters:<br>
					<input type="text" name="simSeq" class="sim" style="display:block;" value='{"lib":"-", "frSz":"500", "frSd":"100", "rdDepth":"50", "rdSz":"100", "rdSd":"0", "errRt":"1", "gcBias":"50"}' />
					<div id="simSeqValMsg" class="warningMessage"></div>
				</div>		
			</div>
			
			<div id="formButtons">
				<hr style="width: 100%; border: 2px solid rgb(150,150,150)" class="w3-round">

				<input type="button" class="button" id="checkFormButton" value="Check input and run project" style="display:block; width:400px;"/>
				
				<div class="checkout" id="checkout-error">
					Your input contains errors. Please review all the red messages in the form. Then, click on "Check input and run project" again.
					<p>If problems persist, please read the documentation.</p>
				</div>
				
				<div class="checkout" id="checkout-success">
					All inputs seem correct but you can review your project now.
					<p id="annReminderMsg" style="display:none;"></p>
					<p>if you are sure you want to start this project, click on the button below.</p>
					<input type="button" class="button" id="runProjectButton" value="Start project" style="display:block;"zzz/>
					<div style="clear:both;"></div>
				</div>
			</div>
		
		</form>
    </div>
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

<?php
// Code for dealing with password
}
?>
