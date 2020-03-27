/*
This .js file simply calls PHP files via AJAX
manage-data-files.htm --> manage-data-files.js --> manage-data-files.php --> BASH commands 

This file contains the function to update input file info on screen and to delete files

This files contains plupload js code to upload files
*/


// Function to communicate html with php via AJAX to retrieve files info 
function filesInfo() {
	var xmlhttp = new XMLHttpRequest();
	xmlhttp.onreadystatechange = function() {
		if (this.readyState == 4 && this.status == 200) {
			document.getElementById("filesInfo").innerHTML = this.responseText;
		}
	};
	xmlhttp.open("GET", "../cgi-bin/manage-input-files-list-files.py", true);
	xmlhttp.send();
}

function ShowRemoveFile(fileName) {
	document.getElementById("rmFile_" + fileName).style.display = "block";
}

function HideRemoveFile(fileName) {
	document.getElementById("rmFile_" + fileName).style.display = "none";
}

// Function to communicate html with php via AJAX to remove a project from disk 
function removeFile(fileName) {
	var xmlhttp = new XMLHttpRequest();
	xmlhttp.onreadystatechange = function() {
		if (this.readyState == 4 && this.status == 200) {
			// Update files info on screen
			filesInfo()
		}
	};
	xmlhttp.open("GET", "../cgi-bin/manage-input-files-removeFile.py?f="+ fileName, true);
	xmlhttp.send();
}

window.onload = function() {
	
	//////////////////////////////////////////////////
	// Here starts plupload js code to upload files //
	//////////////////////////////////////////////////
	
	var uploader = new plupload.Uploader({
		runtimes : 'html5,flash,silverlight,html4',
		browse_button : 'pickfiles', // you can pass an id...
		container: document.getElementById('container'), // ... or DOM Element itself
		url : '../cgi-bin/manage-input-files-upload.py',
		chunk_size: '50mb',
    	max_retries: 3,
		flash_swf_url : 'plupload-2.3.1/js/Moxie.swf',
		silverlight_xap_url : 'plupload-2.3.1/js/Moxie.xap',
	
		filters : {
			max_file_size : '200gb',
			mime_types: [
				{title : "GFF files", extensions : "gff,gff3"},
				{title : "Gene annotation files", extensions : "txt"},
				{title : "FASTA files", extensions : "fa"},
				{title : "FASTQ files", extensions : "fq,fastq"},
				{title : "gzipped files", extensions : "gz"}
			]
		},

		init: {
			PostInit: function() {
				document.getElementById('filelist').innerHTML = '';

				document.getElementById('uploadfiles').onclick = function() {
					uploader.start();
					return false;
				};
			},

			FilesAdded: function(up, files) {
				plupload.each(files, function(file) {
					document.getElementById('filelist').innerHTML += '<div id="' + file.id + '">' + file.name + ' (' + plupload.formatSize(file.size) + ') <b></b></div>';
				});
			},

			UploadProgress: function(up, file) {
				document.getElementById(file.id).getElementsByTagName('b')[0].innerHTML = '<span>' + file.percent + "%</span>";
			},

			Error: function(up, err) {
				document.getElementById('console').appendChild(document.createTextNode("\nError #" + err.code + ": " + err.message));
			}
		}
	});

	uploader.init();
	
	////////////////////////////////////////////////
	// Here ends plupload js code to upload files //
	////////////////////////////////////////////////
	
	
	// Call filesInfo() when page loads
	filesInfo();
	// Call function filesInfo() every 10 seconds to update files status regularly
	setInterval(filesInfo, 10000);
	
}
