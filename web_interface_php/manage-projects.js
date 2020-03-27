/*
This .js file simply calls PHP files via AJAX
easymap.htm --> ajax.js --> fire-wf.php --> easymap.sh --> log.log --> read_log.php ... 

This file contains the function to update project info on screen, to stop running projects,
and to delete from disk finished/killed/error projects


*/

// Enclose all code inside a onload condition?????????????????????????


// Function to communicate html with php via AJAX to retrieve projects info 
function projectsInfo() {
	var xmlhttp = new XMLHttpRequest();
	xmlhttp.onreadystatechange = function() {
		if (this.readyState == 4 && this.status == 200) {
			document.getElementById("projectsInfo").innerHTML = this.responseText;
		}
	};
	
	xmlhttp.open("GET", "manage-projects-list-projects.php", true);
	xmlhttp.send();
}

// Call projectsInfo() 3 times separated by one second. This is done everytime page is
// loaded. Before implementing this, I observed that sometimes, when user was redirected
// from run-new-project.php just after clicking on 'Run new project', the new project
// did not appear in the list. With the current setting, this is solved.
var x = 0;
var intervalID = setInterval(function() {
	projectsInfo();
	if (++x === 3) {
		window.clearInterval(intervalID);
	}
}, 1000);

// Call function projectsInfo() every 50 seconds to update projects status regularly
setInterval(projectsInfo, 50000);

// Function to display <div id="stopPrj_PROJECTNAME">
function displayStopProject(projectName) {
	document.getElementById("stopPrj_" + projectName).style.display = 'block';
}

// Function to hide <div id="stopPrj_PROJECTNAME">
function hideStopProject(projectName) {
	document.getElementById("stopPrj_" + projectName).style.display = 'none';
}

// Function to communicate html with php via AJAX to stop a project 
function stopProject(projectName) {
	var xmlhttp = new XMLHttpRequest();
	xmlhttp.onreadystatechange = function() {
		if (this.readyState == 4 && this.status == 200) {
			// Update projects info on screen
			projectsInfo();
		}
	};
	xmlhttp.open("GET", "stop-project.php?p="+projectName, true);
	xmlhttp.send();
}

// Function to display <div id="rmPrj_PROJECTNAME">
function displayRemoveProject(projectName) {
	document.getElementById("rmPrj_" + projectName).style.display = 'block';
}

// Function to hide <div id="rmPrj_PROJECTNAME">
function hideRemoveProject(projectName) {
	document.getElementById("rmPrj_" + projectName).style.display = 'none';
}

// Function to communicate html with php via AJAX to remove a project from disk 
function removeProject(projectName) {
	var xmlhttp = new XMLHttpRequest();
	xmlhttp.onreadystatechange = function() {
		if (this.readyState == 4 && this.status == 200) {
			// Update projects info on screen
			projectsInfo();
		}
	};
	xmlhttp.open("GET", "remove-project.php?p="+projectName, true);
	xmlhttp.send();
}



