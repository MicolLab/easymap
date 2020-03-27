/*
This .js file simply calls Python CGI file via AJAX
*/

window.onload = function() {

	var projectName = parent.document.URL.substring(parent.document.URL.indexOf('?')+3, parent.document.URL.length);

	// Function to communicate html with php via AJAX to read a log file 
	function readReport(projectName) {
		var xmlhttp = new XMLHttpRequest();
		xmlhttp.onreadystatechange = function() {
			if (this.readyState == 4 && this.status == 200) {
				 document.getElementById("reportInfo").innerHTML = this.responseText;
			}
		};
		xmlhttp.open("GET", "../cgi-bin/view-report-read.py?p=" + projectName, true);
		xmlhttp.send();
	}

	readReport(projectName)
	
}
