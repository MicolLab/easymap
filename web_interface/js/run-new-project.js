/*
This .js file simply calls PHP files via AJAX
easymap.htm --> ajax.js --> xxx.py --> easymap.sh --> log.log --> read_log.php ... 

*/

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// THIS SECTION DEALS WITH GENERAL FUNCTIONALITY OF THE PAGE, SUCH AS CHECKING CONFIG FILE OR TRIGGERING A PROJECT
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Function to communicate html with php via AJAX to check if user_projects directory
// is over config/config>user_projects-size-limit
// This function has to be triggered whe page loads
function allowNewProject() {
	var xmlhttp = new XMLHttpRequest();
	xmlhttp.onreadystatechange = function() {
		if (this.readyState == 4 && this.status == 200) {
			var response = this.responseText;
			var responseFields = response.split(",");
			if (responseFields[0] >= 100) {
				document.getElementById("runNewProject").style.display = "none";
				document.getElementById("sizeWarning").style.display = "block";
				document.getElementById("sizeWarning").innerHTML = "WARNING: You are over the maximum space limit set by the machine administrator (" + responseFields[2] + " Gb). Delete unused input files and/or past projects to free disk space.";
			} else if (responseFields[0] >= 75 && responseFields[0] < 100) {
				document.getElementById("sizeWarning").style.display = "block";
				document.getElementById("sizeWarning").innerHTML = "WARNING: You are over 75% of the maximum space limit set by the machine administrator (" + responseFields[2] + " Gb).";
			}
			if (responseFields[1] >= 100) {
				document.getElementById("runNewProject").style.display = "none";
				document.getElementById("simultWarning").style.display = "block";
				document.getElementById("simultWarning").innerHTML = "WARNING: You have reached the maximum number of simultaneously running projects set by the machine administrator (" + responseFields[3] + " projects). Wait until currently running jobs finish.";
			}
			//console.log(response);
		}
	};
	xmlhttp.open("GET", "../cgi-bin/run-new-project-allow-new.py", true);
	xmlhttp.send();
}

// Function to list all genome reference files to be displayed inside <select multiple> tag
// This function has to be triggered whe page loads
function listInputFiles() {
	var xmlhttp = new XMLHttpRequest();
	xmlhttp.onreadystatechange = function() {
		if (this.readyState == 4 && this.status == 200) {

			// Parse response in JSON format
			var inputFilesresponse = JSON.parse(this.responseText);

			// Create interface to select reference sequence contigs and insertion sequence
			var fastaBasenames = inputFilesresponse[0];
			var refFiles = document.getElementById('refFileSelector');
			if ((fastaBasenames.length < 1)) {
				refFiles.options[refFiles.options.length] = new Option('There are no files with .fa extension', 'n/p');
			} else {
				refFiles.options[refFiles.options.length] = new Option('Select a basename', 'n/p');
				for (i = 0; i < fastaBasenames.length; i++) {
					refFiles.options[refFiles.options.length] = new Option(fastaBasenames[i], fastaBasenames[i]);
				}
			}

			// Create interface to select insertion sequence
			var fastaFiles = inputFilesresponse[1];
			var insFiles = document.getElementById('insFileSelector');
			if (fastaFiles.length < 1) {
				insFiles.options[insFiles.options.length] = new Option('There are no files with .fa extension', 'n/p');
			} else {
				insFiles.options[insFiles.options.length] = new Option('Select a file', 'n/p');
				for (i = 0; i < fastaFiles.length; i++) {
					insFiles.options[insFiles.options.length] = new Option(fastaFiles[i], fastaFiles[i]);
				}
			}

			// Create interface to select GFF and ANN files
			var otherFiles = inputFilesresponse[3];
			var gffFiles = document.getElementById('gffFileSelector');
			var annFiles = document.getElementById('annFileSelector');
			if (otherFiles.length < 1) {
				gffFiles.options[gffFiles.options.length] = new Option('There are no files with the appropriate extension', 'n/p');
				annFiles.options[annFiles.options.length] = new Option('There are no files with the appropriate extension');
			} else {
				gffFiles.options[gffFiles.options.length] = new Option('Select a file', 'n/p');
				annFiles.options[annFiles.options.length] = new Option('Select a file', 'n/p');
				for (i = 0; i < otherFiles.length; i++) {
					gffFiles.options[gffFiles.options.length] = new Option(otherFiles[i], otherFiles[i]);
					annFiles.options[annFiles.options.length] = new Option(otherFiles[i], otherFiles[i]);
				}
			}

			// Create interfaces to select fastq files
			var fastqFiles = inputFilesresponse[2];
			var readsProblem = document.getElementById('readsProblemSelector');
			var readsControl = document.getElementById('readsControlSelector');
			if (fastqFiles.length < 1) {
				readsProblem.options[readsProblem.options.length] = new Option('There are no files with .fq extension', 'n/p');
				readsControl.options[readsControl.options.length] = new Option('There are no files with .fq extension', 'n/p');
			} else {
				for (i = 0; i < fastqFiles.length; i++) {
					readsProblem.options[readsProblem.options.length] = new Option(fastqFiles[i], fastqFiles[i]);
					readsControl.options[readsControl.options.length] = new Option(fastqFiles[i], fastqFiles[i]);
				}
			}
		}
	};
	xmlhttp.open("GET", "../cgi-bin/run-new-project-list-input-files.py", true);
	xmlhttp.send();
}

// Function to check if a string has a valid JSON format
function isJsonString(str) {
	try {
		JSON.parse(str);
	} catch (e) {
		return false;
	}
	return true;
}

// Function to determine if an array has duplicated values
// Used to know if the user has chosen common read files for the problem and control samples
function HasDuplicates(array) {
    return (new Set(array)).size !== array.length;
}

// Trigger required functions when page loads
allowNewProject();
listInputFiles()

function resetTextField() {
	this.value = "";
}

function HideCheckoutBoxes() {
	var CheckoutBoxes = document.getElementsByClassName("checkout");
	for (var i=0; i<CheckoutBoxes.length; i++) {
		CheckoutBoxes[i].style.display = "none";
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// THIS SECTION DEALS WITH THE DYNAMIC FORMATTING OF THE FORM AND WITH FIELDS VALIDATION
///////////////////////////////////////////////////////////////////////////////////////////////////////

window.onload = function() {
	
	// Functions that need to be declared after document has completely loaded
	
	// Uptdate the command arguments in the screen. For development purposes.
	// This function is called in many other functions after updating the value of an argument.
	function updateCmd() {
			console.log(cmdArgs);
			//document.getElementById("commandString").innerHTML = cmdArgs;
	}
	
	// Check if project name is suitable for easymap
	function verifyProjectName(){
		HideCheckoutBoxes();

		var text = document.getElementById("form1").projectName.value;
		if(/[^a-zA-Z0-9]/.test( text) ) {
			cmdArgs[1] = 'n/p';
			projectNameValidationInfoMessage = 'Input is not alphanumeric';
			document.getElementById("projectNameValidationInfo").innerHTML = projectNameValidationInfoMessage;
			document.getElementById("projectNameValidationInfo").style.display = "block";
		} else if (text == '') {
			cmdArgs[1] = 'n/p';
			projectNameValidationInfoMessage = 'You must give a name to the project';
			document.getElementById("projectNameValidationInfo").innerHTML = projectNameValidationInfoMessage;
			document.getElementById("projectNameValidationInfo").style.display = "block";
		} else {
			cmdArgs[1] = document.getElementById("form1").projectName.value;
			document.getElementById("projectNameValidationInfo").style.display = "none";
		}
		//updateCmd();
	}
	
	// Determine option button selected and define the appropriate command argument
	function buttons_analysisType() {
		HideCheckoutBoxes();

		var options = document.getElementsByClassName("analysisType");
		for (var i=0; i<options.length; i++) {
			if (options[i].checked == true) {
				var checkedOption = options[i].id;
			}
		}
		if (checkedOption == 'button1') {
			cmdArgs[2] = 'ins';
			document.getElementById("insSeqField").style.display = "block";
			document.getElementById("readsControl").style.display = "none";
			document.getElementById("backgroundCrossCtype").style.display = "none";
			document.getElementById("simRecselInterface").style.display = "none";
			if (cmdArgs[3] == 'exp') {
				document.getElementById("expDataInterface").style.display = "block";
			}
			if (cmdArgs[3] == 'sim') {
				document.getElementById("simDataInterface").style.display = "block";
			}
		} else {
			cmdArgs[2] = 'snp';
			document.getElementById("insSeqField").style.display = "none";
			document.getElementById("readsControl").style.display = "block";
			document.getElementById("backgroundCrossCtype").style.display = "block";
			document.getElementById("simRecselInterface").style.display = "block";
			if (cmdArgs[3] == 'exp') {
				document.getElementById("expDataInterface").style.display = "block";
			}
			if (cmdArgs[3] == 'sim') {
				document.getElementById("simDataInterface").style.display = "block";
			}
		}
		//updateCmd();
		document.getElementById("analysisTypeValidationInfo").style.display = "none";		
	}
	
	// Determine option button selected and define the appropriate command argument
	function buttons_dataSource() {
		HideCheckoutBoxes();

		var options = document.getElementsByClassName("dataSource");
		for (var i=0; i<options.length; i++) {
			if (options[i].checked == true) {
				var checkedOption = options[i].id;
			}
		}
		if (checkedOption == 'button3') {
			cmdArgs[3] = 'exp';
			document.getElementById("simDataInterface").style.display = "none";
			if (cmdArgs[2] != 'n/p') {
				document.getElementById("expDataInterface").style.display = "block";
			}
		} else {
			cmdArgs[3] = 'sim';
			document.getElementById("expDataInterface").style.display = "none";
			if (cmdArgs[2] != 'n/p') {
				document.getElementById("simDataInterface").style.display = "block";
			}
		}
		//updateCmd();
		document.getElementById("dataSourceValidationInfo").style.display = "none";
	}
	
	/* OLD METHOD TO CHECK MULTIPLE REFERENCE SEQUENCE FILES
	// Determine all the reference file names selected, add them to array, and then to command argument
	function checkRefSeqs() {
		var contigs = document.getElementById("refFileSelector");
		var contigsList = [];
		
		for (var i=0; i<contigs.length; i++) {
			if (contigs[i].selected == true) {
				contigsList.push(contigs[i].value);
			}
		}
		cmdArgs[4] = contigsList;
		//updateCmd();
		document.getElementById("refSeqValidationInfo").style.display = "none";
	}
	*/

	// Update command arguments after each user interaction with sinlge selectors
	function processSingleSelectors() {
		HideCheckoutBoxes();

		if (this.id == 'refFileSelector') {
			cmdArgs[4] = this.value;
			document.getElementById("refSeqValidationInfo").style.display = "none";
		}
		if (this.id == 'insFileSelector') {
			cmdArgs[5] = this.value;
			document.getElementById("insFileValidationInfo").style.display = "none";
		}
		if (this.id == 'gffFileSelector') {
			cmdArgs[6] = this.value;
			document.getElementById("gffFileValidationInfo").style.display = "none";
		}
		if (this.id == 'annFileSelector') {
			cmdArgs[7] = this.value;
		}
		//updateCmd();
		document.getElementById("annReminderMsg").style.display = "none";	
	}

	// Mutant background: determine option button selected and define the appropriate command argument
	function buttons_mutBackground() {
		HideCheckoutBoxes();

		var options = document.getElementsByClassName("mutBackground");
		for (var i=0; i<options.length; i++) {
			if (options[i].checked == true) {
				var checkedOption = options[i].id;
			}
		}
		if (checkedOption == 'button11') {
			cmdArgs[16] = 'ref';
		} else {
			cmdArgs[16] = 'noref';
		}
		//updateCmd();
		document.getElementById("mutBackgroundValidationInfo").style.display = "none";
	}

	// Mapping cross preformed: determine option button selected and define the appropriate command argument
	function buttons_crossType() {
		HideCheckoutBoxes();

		var options = document.getElementsByClassName("crossType");
		for (var i=0; i<options.length; i++) {
			if (options[i].checked == true) {
				var checkedOption = options[i].id;
			}
		}
		if (checkedOption == 'button13') {
			cmdArgs[17] = 'bc';
		} else {
			cmdArgs[17] = 'oc';
		}
		//updateCmd();
		document.getElementById("crossTypeValidationInfo").style.display = "none";
	}

	// Origin of the control reads: determine option button selected and define the appropriate command argument
	function buttons_contType() {
		HideCheckoutBoxes();

		var options = document.getElementsByClassName("contType");
		for (var i=0; i<options.length; i++) {
			if (options[i].checked == true) {
				var checkedOption = options[i].id;
			}
		}
		if (checkedOption == 'button15') {
			cmdArgs[18] = 'par';
			cmdArgs[19] = 'mut';
		} else if (checkedOption == 'button16') {
			cmdArgs[18] = 'par';
			cmdArgs[19] = 'nomut';
		} else {
			cmdArgs[18] = 'f2wt';
			cmdArgs[19] = 'n/p';
		}
		//updateCmd();
		document.getElementById("contTypeValidationInfo").style.display = "none";
	}

	// Check if combination of mutant background, cross performed, and origin of control reads, is supported
	function checkBackgroundCrossCtypeIntermediateCheck() {
		HideCheckoutBoxes();

		var experimentalDesign = cmdArgs[16] + '_' + cmdArgs[17] + '_' + cmdArgs[18] + '_' + cmdArgs[19];

		var allowedExperimentalDesigns = [
			'ref_bc_par_mut',
			'ref_bc_f2wt_n/p',
			'ref_oc_par_mut',
			'ref_oc_par_nomut',
			'noref_bc_f2wt_n/p',
			'noref_oc_par_mut',
			'ref_bc_n/p_n/p', 'ref_oc_n/p_n/p', 'ref_n/p_n/p_n/p',
			'noref_bc_n/p_n/p', 'noref_oc_n/p_n/p', 'noref_n/p_n/p_n/p',
			'n/p_n/p_n/p_n/p'
		];

		//var found = allowedExperimentalDesigns.includes(experimentalDesign);

		if (allowedExperimentalDesigns.includes(experimentalDesign)) {
			document.getElementById("backgroundCrossCtypeWarnMsg").style.display = "none";
		} else {
			document.getElementById("backgroundCrossCtypeWarnMsg").style.display = "block";
		}
	}

	// Check reads selectors (max two files selected per sample) and update command arguments
	function checkProblemReads() {
		HideCheckoutBoxes();

		var reads = document.getElementById("readsProblemSelector");
		var readsList = [];
		
		for (var i=0; i<reads.length; i++) {
			if (reads[i].selected == true) {
				readsList.push(reads[i].value);
			}
		}

		if (readsList.length == 1) {
			cmdArgs[8] = readsList[0]; cmdArgs[9] = 'n/p'; cmdArgs[10] = 'n/p'; cmdArgs[11] = 'se';
			//updateCmd();
			// Hide error message
			document.getElementById("readsProblemWarnMsg").style.display = "none";
		} else if (readsList.length == 2) {
			cmdArgs[8] = 'n/p'; cmdArgs[9] = readsList[0]; cmdArgs[10] = readsList[1]; cmdArgs[11] = 'pe';
			//updateCmd();
			// Hide error message
			document.getElementById("readsProblemWarnMsg").style.display = "none";
		} else {
			cmdArgs[8] = 'n/p'; cmdArgs[9] = 'n/p'; cmdArgs[10] = 'n/p'; cmdArgs[11] = 'n/p';
			//updateCmd();
			// Display error message
			document.getElementById("readsProblemWarnMsg").style.display = "block";
		}
		document.getElementById("readsProblemWarnMsg2").style.display = "none";
	}

	function checkControlReads() {
		HideCheckoutBoxes();

		var readsC = document.getElementById("readsControlSelector");
		var readsListC = [];
		
		for (var i=0; i<readsC.length; i++) {
			if (readsC[i].selected == true) {
				readsListC.push(readsC[i].value);
			}
		}

		if (readsListC.length == 1) {
			cmdArgs[12] = readsListC[0]; cmdArgs[13] = 'n/p'; cmdArgs[14] = 'n/p'; cmdArgs[15] = 'se';
			//updateCmd();
			// Hide error message
			document.getElementById("readsControlWarnMsg").style.display = "none";
		} else if (readsListC.length == 2) {
			cmdArgs[12] = 'n/p'; cmdArgs[13] = readsListC[0]; cmdArgs[14] = readsListC[1]; cmdArgs[15] = 'pe';
			//updateCmd();
			// Hide error message
			document.getElementById("readsControlWarnMsg").style.display = "none";
		} else {
			cmdArgs[12] = 'n/p'; cmdArgs[13] = 'n/p'; cmdArgs[14] = 'n/p'; cmdArgs[15] = 'n/p';
			//updateCmd();
			// Show error message
			document.getElementById("readsControlWarnMsg").style.display = "block";
		}
		document.getElementById("readsControlWarnMsg2").style.display = "none";
	}

	function checkStringency() {
		if (this.checked) {
			cmdArgs[23] = 'low_stringency';
		} else {
			cmdArgs[23] = 'high_stringency';
		}
	}

	function verifySimMut() {
		HideCheckoutBoxes();

		var simMutInput = document.getElementById("form1").simMut.value;
		if (isJsonString(simMutInput) == false) {
			document.getElementById("simMutValMsg").innerHTML = 'The structure of the input is not correct.';
			document.getElementById("simMutValMsg").style.display = "block";
			cmdArgs[20] = 'n/p';
		} else {
			var simMutInput = JSON.parse(simMutInput);
			if (/[^0-9]/.test(simMutInput.numberMutations) || simMutInput.numberMutations < 1) {
				document.getElementById("simMutValMsg").innerHTML = '"numberMutations" must be a positive integer.';
				document.getElementById("simMutValMsg").style.display = "block";
				cmdArgs[20] = 'n/p';
			} else {
				document.getElementById("simMutValMsg").style.display = "none";
				if (cmdArgs[2] == 'ins') {
					cmdArgs[20] = simMutInput.numberMutations + "~li";
				} else {
					cmdArgs[20] = simMutInput.numberMutations + "~e";
				}
			}
		}
		//updateCmd();
	}

	function verifySimrecselFieldA() {
		HideCheckoutBoxes();

		var simRecselInputA = document.getElementById("form1").simRecselA.value;
		var simRecselErrors = ["Errors:"]
		var simRecselErr = false;

		if (isJsonString(simRecselInputA) == false) {
			document.getElementById("simRecselAValMsg").innerHTML = 'The structure of the input is not correct.';
			document.getElementById("simRecselAValMsg").style.display = "block";
			return false;
		} else {
			var simRecselInputA = JSON.parse(simRecselInputA);
			
			if (/[^0-9]/.test(simRecselInputA.contigCausalMut) || simRecselInputA.contigCausalMut < 1) {
				simRecselErr = true;
				simRecselErrors.push('"contigCausalMut" must be an integer equal or greater than 1.');
			}
			if (/[^0-9]/.test(simRecselInputA.posCausalMut) || simRecselInputA.posCausalMut < 1) {
				simRecselErr = true;
				simRecselErrors.push('"posCausalMut" must be an integer equal or greater than 1.');
			}
			if (/[^0-9]/.test(simRecselInputA.numRecChrs) || simRecselInputA.numRecChrs < 1) {
				simRecselErr = true;
				simRecselErrors.push('"numRecChrs" must be an integer equal or greater than 1.');
			}

			if (simRecselErr == true) {
				document.getElementById("simRecselAValMsg").innerHTML = simRecselErrors.join("<br>");
				document.getElementById("simRecselAValMsg").style.display = "block";
				return false;
			} else {
				document.getElementById("simRecselAValMsg").style.display = "none";
				return true;
			}
		}
	}

	function verifySimrecselFieldB() {
		HideCheckoutBoxes();

		var simRecselInputB = document.getElementById("form1").simRecselB.value;
		var simRecselErrors = ["Errors:"]
		var simRecselErr = false;

		if (isJsonString(simRecselInputB) == false) {
			document.getElementById("simReqselBValMsg").innerHTML = 'The structure of the input is not correct.';
			document.getElementById("simReqselBValMsg").style.display = "block";
			return false;
		} else {
			var simRecselInputB = JSON.parse(simRecselInputB);

			for (var key in simRecselInputB) {				
				var recFreqs = simRecselInputB[key];
				var ArrRecFreqs = recFreqs.split("-");
				if (ArrRecFreqs.length < 2) {
					simRecselErr = true;
					simRecselErrors.push('The recombination frequencies of contig ' + key + ' have an incorrect format.');
				} else {
					// check every element against the pattern nbr,nbr
					for (i = 0; i < ArrRecFreqs.length; i++) {
						if (!/^\d{1,2}\,\d{1,2}$/.test(ArrRecFreqs[i])) {
							simRecselErr = true;
							simRecselErrors.push('The following parameter has a wrong format: ' + ArrRecFreqs[i] + '.');
						}
					}
				}
			}

			if (simRecselErr == true) {
				document.getElementById("simReqselBValMsg").innerHTML = simRecselErrors.join("<br>");
				document.getElementById("simReqselBValMsg").style.display = "block";
				return false;
			} else {
				document.getElementById("simReqselBValMsg").style.display = "none";
				return true;
			}
		}
	}

	function verifySimSeq() {
		HideCheckoutBoxes();
		
		var simSeqInput = document.getElementById("form1").simSeq.value;
		var simSeqErrors = ["Errors:"]
		var simSeqErr = false;

		if (isJsonString(simSeqInput) == false) {
			document.getElementById("simSeqValMsg").innerHTML = 'The structure of the input is not correct.';
			document.getElementById("simSeqValMsg").style.display = "block";
			cmdArgs[22] = 'n/p'; cmdArgs[11] = 'n/p'; cmdArgs[15] = 'n/p';
		} else {
			var simSeqInput = JSON.parse(simSeqInput);
			
			if (["se","pe"].indexOf(simSeqInput.lib) === -1) {
				simSeqErr = true;
				simSeqErrors.push('"lib" must be equal to "se" (single-and library) or "pe" (paired-end library).');
			}
			if (/[^0-9]/.test(simSeqInput.frSz) || simSeqInput.frSz < 50) {
				simSeqErr = true;
				simSeqErrors.push('"frSz" must be an integer equal or greater than 50.');
			}
			if (/[^0-9]/.test(simSeqInput.frSd) || simSeqInput.frSd < 0 || simSeqInput.frSd > simSeqInput.frSz) {
				simSeqErr = true;
				simSeqErrors.push('"frSd" must be an integer in the interval [0-"frSz"].');
			}
			if (/[^0-9]/.test(simSeqInput.rdDepth) || simSeqInput.rdDepth < 1) {
				simSeqErr = true;
				simSeqErrors.push('"rdDepth" must be a positive integer.');
			}
			if (/[^0-9]/.test(simSeqInput.rdSz) || simSeqInput.rdSz < 20 || simSeqInput.rdSz > simSeqInput.frSz) {
				simSeqErr = true;
				simSeqErrors.push('"rdSz" must be an integer equal or greater than 20 and smaller or equal to "frSz".');
			}
			if (/[^0-9]/.test(simSeqInput.rdSd) || simSeqInput.rdSd < 0 || simSeqInput.rdSd > simSeqInput.rdSz) {
				simSeqErr = true;
				simSeqErrors.push('"rdSd" must be an integer in the interval [0-"rdSz"].');
			}
			if (/[^0-9]/.test(simSeqInput.errRt) || simSeqInput.errRt < 0 || simSeqInput.errRt > 10) {
				simSeqErr = true;
				simSeqErrors.push('"errRt" must be an integer in the interva [0-10].');
			}
			if (/[^0-9]/.test(simSeqInput.gcBias) || simSeqInput.gcBias < 0 || simSeqInput.gcBias > 100) {
				simSeqErr = true;
				simSeqErrors.push('"gcBias" must be an integer in the interval [0-100].');
			}

			if (simSeqErr == true) {
				cmdArgs[22] = 'n/p'; cmdArgs[11] = 'n/p'; cmdArgs[15] = 'n/p';
				var simSeqErrorsString = simSeqErrors.join("<br>");
				document.getElementById("simSeqValMsg").innerHTML = simSeqErrorsString;
				document.getElementById("simSeqValMsg").style.display = "block";
			} else {
				document.getElementById("simSeqValMsg").style.display = "none";
				cmdArgs[22] = simSeqInput.rdDepth + "~" + simSeqInput.rdSz + "," + simSeqInput.rdSd + "~" + simSeqInput.frSz + "," + simSeqInput.frSd + "~" + simSeqInput.errRt + "~" + simSeqInput.gcBias + "~" + simSeqInput.lib;
				cmdArgs[11] = simSeqInput.lib; cmdArgs[15] = simSeqInput.lib;
			}
		}
		//updateCmd();
	}

	// All intermediate checks require the user to interact with the different form fields. However, the user 
	// could skip some elements by mistake. Therefore, I force the user to click in a button that executes the 
	// following function, which checks the result of all form fields.
	function commandFinalCheck() {
		var userErrors = false;

		// Check that project name has been set
		if (cmdArgs[1] == 'n/p') {
			userErrors = true;
			projectNameValidationInfoMessage = 'You must give a name to the project.';
			document.getElementById("projectNameValidationInfo").innerHTML = projectNameValidationInfoMessage;
			document.getElementById("projectNameValidationInfo").style.display = "block";
		}

		// Check that mapping by sequencing strategy has been set
		if (cmdArgs[2] == 'n/p') {
			userErrors = true;
			analysisTypeValidationInfoMessage = 'You must choose a mapping by sequencing strategy.';
			document.getElementById("analysisTypeValidationInfo").innerHTML = analysisTypeValidationInfoMessage;
			document.getElementById("analysisTypeValidationInfo").style.display = "block";
		}

		// Check that data source has been set
		if (cmdArgs[3] == 'n/p') {
			userErrors = true;
			dataSourceValidationInfoMessage = 'You must choose a data source.';
			document.getElementById("dataSourceValidationInfo").innerHTML = dataSourceValidationInfoMessage;
			document.getElementById("dataSourceValidationInfo").style.display = "block";
		}

		// Check that reference sequence has been set
		if (cmdArgs[4] == 'n/p') {
			userErrors = true;
			refSeqValidationInfoMessage = 'You must select one or more reference sequence files.';
			document.getElementById("refSeqValidationInfo").innerHTML = refSeqValidationInfoMessage;
			document.getElementById("refSeqValidationInfo").style.display = "block";
		}

		// If user chose tagged sequence mapping, check that an insertion sequence file has been selected
		if (cmdArgs[2] == 'ins' && cmdArgs[5] == 'n/p') {
			userErrors = true;
			insFileValidationInfoMessage = 'You must select an insertion sequence file.';
			document.getElementById("insFileValidationInfo").innerHTML = insFileValidationInfoMessage;
			document.getElementById("insFileValidationInfo").style.display = "block";
		}

		// Check that gff file has been set
		if (cmdArgs[6] == 'n/p') {
			userErrors = true;
			gffFileValidationInfoMessage = 'You must select a GFF file that matches the names the reference sequence(s) selected.';
			document.getElementById("gffFileValidationInfo").innerHTML = gffFileValidationInfoMessage;
			document.getElementById("gffFileValidationInfo").style.display = "block";
		}

		// Determine if ann file has been set
		if (cmdArgs[7] == 'n/p') {
			document.getElementById("annReminderMsg").innerHTML = 'REMINDER: You did not select any gene functional annotation file. While easymap can run without it, this information can help to interpret the final results. To know more about where to find this information and how to format it for easymap, see the documentation';
			document.getElementById("annReminderMsg").style.display = "block";
		}

		// If user chose MbS analysis, check if all two-way selectors were clicked on
		if (cmdArgs[2] == 'snp') {
			if (cmdArgs[16] == 'n/p') {
				userErrors = true;
				document.getElementById("mutBackgroundValidationInfo").innerHTML = 'You must select a mutant background.';
				document.getElementById("mutBackgroundValidationInfo").style.display = "block";
			}
			if (cmdArgs[17] == 'n/p') {
				userErrors = true;
				document.getElementById("crossTypeValidationInfo").innerHTML = 'You must select the mapping cross performed.';
				document.getElementById("crossTypeValidationInfo").style.display = "block";
			}
			if (cmdArgs[18] == 'n/p') {
				userErrors = true;
				document.getElementById("contTypeValidationInfo").innerHTML = 'You must select the origin of the control reads.';
				document.getElementById("contTypeValidationInfo").style.display = "block";
			}
		}

		// If user chose own experimental data, check if problem reads file(s) have been specified
		if (cmdArgs[3] == 'exp' && cmdArgs[11] == 'n/p') {
			userErrors = true;
			document.getElementById("readsProblemWarnMsg").style.display = "block";
		}

		// If user chose own experimental data and MbS analysis, check if control reads file(s) have been specified
		if (cmdArgs[2] == 'snp' && cmdArgs[3] == 'exp' && cmdArgs[15] == 'n/p') {
			userErrors = true;
			document.getElementById("readsControlWarnMsg").style.display = "block";
		}

		// Check if user has selected, by mistake, one or more common files as the problem and control reads
		if (cmdArgs[2] == 'snp' && cmdArgs[3] == 'exp') {
			var AllReadArgs = [cmdArgs[8], cmdArgs[9], cmdArgs[10], cmdArgs[12], cmdArgs[13], cmdArgs[14]];
			var AllReads = AllReadArgs.filter(function(Read) {
			  return Read != 'n/p';
			});
			if (HasDuplicates(AllReads)) {
				userErrors = true;
				document.getElementById("readsProblemWarnMsg2").style.display = "block";
				document.getElementById("readsControlWarnMsg2").style.display = "block";
			} else {
				document.getElementById("readsProblemWarnMsg2").style.display = "none";
				document.getElementById("readsControlWarnMsg2").style.display = "none";
			}
		}

		// If user chose simulated data, check if simMut and simSeq parameters have been set properly
		if (cmdArgs[3] == 'sim' && cmdArgs[20] == 'n/p') {
			userErrors = true;
			document.getElementById("simMutValMsg").innerHTML = 'The input in this field is not correct.';
			document.getElementById("simMutValMsg").style.display = "block";
		}

		if (cmdArgs[3] == 'sim' && cmdArgs[22] == 'n/p') {
			userErrors = true;
			document.getElementById("simSeqValMsg").innerHTML = 'The input in this field is not correct.';
			document.getElementById("simSeqValMsg").style.display = "block";
		}

		// If user chose MbS simulated data, check simRecsel parameters
		if (cmdArgs[2] == 'snp' && cmdArgs[3] == 'sim') {
			if (verifySimrecselFieldA() == true && verifySimrecselFieldB() == true) {
				var InA = JSON.parse(document.getElementById("form1").simRecselA.value);
				var InB = JSON.parse(document.getElementById("form1").simRecselB.value);
				var str2 = Object.keys(InB).map(function(k){return InB[k]}).join("/");
				var str1 = InA.contigCausalMut + ',' + InA.posCausalMut + '~r~' + InA.numRecChrs;
				cmdArgs[21] = str2 + '~' + str1;
				//updateCmd();
			}
			if (verifySimrecselFieldA() == false) {
				userErrors = true;
				document.getElementById("simRecselAValMsg").innerHTML = 'The input in this field is not correct.';
				document.getElementById("simRecselAValMsg").style.display = "block";
			}
			if (verifySimrecselFieldB() == false) {
				userErrors = true;
				document.getElementById("simReqselBValMsg").innerHTML = 'The input in this field is not correct.';
				document.getElementById("simReqselBValMsg").style.display = "block";
			}
		}

		// In snp mode, check argument 23 (stringency). If user did not interact with the switch, select the default value 'high_stringency'
		if (cmdArgs[2] == 'snp' && cmdArgs[23] == 'n/p') {
			cmdArgs[23] = 'high_stringency';
		}

		if (userErrors == true) {
			document.getElementById("checkout-error").style.display = "block";
			document.getElementById("checkout-success").style.display = "none";
		} else {
			document.getElementById("checkout-error").style.display = "none";
			document.getElementById("checkout-success").style.display = "block";
		}		
	}

	function sleep(ms) {
		return new Promise(resolve => setTimeout(resolve, ms));
	}

	function goToManageProjects() {
		window.location.assign("manage-projects.htm");
	}

	// Function to trigger a new easymap execution and to redirect browser to manage-projects.htm
	async function runProject() {
		var http = new XMLHttpRequest();
		var url = "../cgi-bin/run-new-project-create-command.py";

		// Create POST string to send
		var argsStringToPost = "program=" + cmdArgs[0] +
							   "&project_name=" + cmdArgs[1] +
							   "&workflow=" + cmdArgs[2] +
							   "&data_source=" + cmdArgs[3] +
							   "&ref_seq=" + cmdArgs[4] +
							   "&ins_seq=" + cmdArgs[5] +
							   "&gff_file=" + cmdArgs[6] +
							   "&ann_file=" + cmdArgs[7] +
							   "&read_s=" + cmdArgs[8] +
							   "&read_f=" + cmdArgs[9] +
							   "&read_r=" + cmdArgs[10] +
							   "&lib_type_sample=" + cmdArgs[11] +
							   "&read_s_ctrl=" + cmdArgs[12] +
							   "&read_f_ctrl=" + cmdArgs[13] +
							   "&read_r_ctrl=" + cmdArgs[14] +
							   "&lib_type_ctrl=" + cmdArgs[15] +
							   "&is_ref_strain=" + cmdArgs[16] +
							   "&cross_type=" + cmdArgs[17] +
							   "&snp_analysis_type=" + cmdArgs[18] +
							   "&control_parental=" + cmdArgs[19] +
							   "&sim_mut=" + cmdArgs[20] +
							   "&sim_recsel=" + cmdArgs[21] +
							   "&sim_seq=" + cmdArgs[22] +
							   "&stringency=" + cmdArgs[23];

		//console.log('argsStringToPost: ' + argsStringToPost);

		http.open("POST", url, true);
		// Request header
		http.setRequestHeader("Content-type", "application/x-www-form-urlencoded");

		http.onreadystatechange = function() { // Call a function when the state changes
		    if(http.readyState == 4 && http.status == 200) {
		        //console.log(http.responseText);
		    }
		}
		// Send request
		http.send(argsStringToPost);

		await sleep(1000);
		// Redirect browser to projects page
		goToManageProjects()
	}

	/* OLD WAYS OF RUNNING NEW PROJECT
	// Function to trigger a new easymap execution and to redirect browser to manage-projects.htm 
	function runProject() {
		var cmdArgs = ['aaa', 'bbb', 'ccc'];
		var xhr = new XMLHttpRequest();
		xhr.open("POST", "run-new-project-create-command.py", true);
		xhr.setRequestHeader("Content-type", "application/json");
		xhr.send(JSON.stringify(cmdArgs));
		console.log('clicked');
		//Maybe desirable: A way to check if the request was processed correctly by php
		xhr.onreadystatechange = function () {
		    if (xhr.readyState === 4 && xhr.status === 200) {
		        if (this.responseText == 'success') {
		        	//
		        } else {
		        	alert('Error: the server that hosts easymap could not complete your request.');
		        }
		    }
		};
		//goToManageProjects();		
	}

	<a href="manage-projects.htm" class="button" onclick="runProject()">Run workflow</a>

	function runProject() {
		//alert('Button pressed');
		var xmlhttp = new XMLHttpRequest();
		xmlhttp.onreadystatechange = function() {
			if (this.readyState == 4 && this.status == 200) {
				document.getElementById("commandResponse").innerHTML = this.responseText;
			}
		};
		
		var arguments = "test-arg";
		
		xmlhttp.open("GET", "run-new-project.php?args=" + arguments, true);
		xmlhttp.send();
	}
	*/

	// End of functions ***************************************************************************************************
	
	
	// Define array with all the command arguments
/*	var cmdArgs = ['./easymap.sh','project_name','workflow','data_source','ref_seq','ins_seq','gff_file','ann_file',
					'read_s','read_f','read_r','lib_type_sample',
					'read_s_ctrl','read_f_ctrl','read_r_ctrl','lib_type_ctrl',
					'is_ref_strain','cross_type','snp_analysis_type','control_parental',
					'sim_mut','sim_recsel','sim_seq','stringency'];
*/	
	var cmdArgs = ['./easymap.sh','n/p','n/p','n/p','n/p','n/p','n/p','n/p',
					'n/p','n/p','n/p','n/p',
					'n/p','n/p','n/p','n/p',
					'n/p','n/p','n/p','n/p',
					'n/p','n/p','n/p','n/p'];

	// Create the command string for the first time (for development purposes only)
	//updateCmd();

	// React to interactions with text inputs
	// Reset default content when user clicks on input box
	document.getElementById("form1").projectName.onfocus = resetTextField;
	
	// Verify input of text fields
	document.getElementById("form1").projectName.onblur = verifyProjectName;
	
	// React to interactions with the main 2-way selectors
	document.getElementById("button1").onclick = buttons_analysisType;
	document.getElementById("button2").onclick = buttons_analysisType;
	document.getElementById("button3").onclick = buttons_dataSource;
	document.getElementById("button4").onclick = buttons_dataSource;
	
	// React to single selectors (refSeq, insFile, gffFile, annFile...)
	document.getElementById("form1").refFileSelector.onblur = processSingleSelectors;
	document.getElementById("form1").insFileSelector.onblur = processSingleSelectors;
	document.getElementById("form1").gffFileSelector.onblur = processSingleSelectors;
	document.getElementById("form1").annFileSelector.onblur = processSingleSelectors;
	
	// React to interactions with the MbS 2-way selectors
	document.getElementById("button11").onclick = buttons_mutBackground;
	document.getElementById("button12").onclick = buttons_mutBackground;
	document.getElementById("button13").onclick = buttons_crossType;
	document.getElementById("button14").onclick = buttons_crossType;
	document.getElementById("button15").onclick = buttons_contType;
	document.getElementById("button16").onclick = buttons_contType;
	document.getElementById("button17").onclick = buttons_contType;

	//React to interactions with reads selectors
	document.getElementById("form1").readsProblemSelector.onclick = checkProblemReads;
	document.getElementById("form1").readsControlSelector.onclick = checkControlReads;

	//React to interactions with simulation fields
	document.getElementById("form1").simMut.onblur = verifySimMut;
	document.getElementById("form1").simSeq.onblur = verifySimSeq;
	document.getElementById("form1").simRecselA.onblur = verifySimrecselFieldA;
	document.getElementById("form1").simRecselB.onblur = verifySimrecselFieldB;

	// React to interactions with backgroundCrossCtype buttons
	document.getElementById("backgroundCrossCtype").onclick = checkBackgroundCrossCtypeIntermediateCheck;

	// React to interactions with stringency button
	document.getElementById("stringency").onclick = checkStringency;

	// React to interactions with button to check input and go to the gateway to run a new project
	document.getElementById("checkFormButton").onclick = commandFinalCheck;

	// React to interactions with button to run project
	document.getElementById("runProjectButton").onclick = runProject;
}


