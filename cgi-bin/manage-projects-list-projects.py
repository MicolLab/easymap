#!src/Python-2.7.18/.localpython/bin/python2

import os, cgi, cgitb, math
cgitb.enable()

#form = cgi.FieldStorage()
dirname = 'user_projects'

# Get list of files
dirpaths = []
for basename in os.listdir(dirname):
    filename = os.path.join(dirname, basename)
    if os.path.isdir(filename):
        dirpaths.append(filename)

# Funtion to get directory size
def get_size(start_path):
    total_size = 0
    for dirpath, dirnames, filenames in os.walk(start_path):
        for f in filenames:
            fp = os.path.join(dirpath, f)
            total_size += os.path.getsize(fp)
    return total_size

# Re-populate list with filename, size tuples
for i in xrange(len(dirpaths)):
    dirpaths[i] = (dirpaths[i], get_size(dirpaths[i]))

# Funtion to manage file sizes
def convert_size(size_bytes):
   if size_bytes == 0:
       return "0B"
   size_name = ("B", "KB", "MB", "GB", "TB", "PB", "EB", "ZB", "YB")
   i = int(math.floor(math.log(size_bytes, 1024)))
   p = math.pow(1024, i)
   s = round(size_bytes / p, 1)
   return "%s %s" % (s, size_name[i])

# Fill html file
print "Content-type:text/html\r\n\r\n" 


if len(dirpaths)==0:
		print '<div class="emptyResult">There are no projects currently on the server. To start a new one go to "Run new project".</div>'

else:
	for f in reversed(sorted(dirpaths, key=lambda x: x[0])):
		project = str(f[0].split('user_projects/')[1])
		folder_size = convert_size(int(f[1]))
		try:
			with open('./user_projects/'+project+'/2_logs/status') as status_file:
				for line in reversed(status_file.readlines()):
					if line.startswith('status:'): 
						status = str(line.split(':')[1]).strip()
						break
		except:
			status='Status file not available'

		if status != 'Status file not available':
			print '''
				<!-- <div style="background-color:#79e59d; border:solid green 1px; border-radius:5px; padding:10px;"> -->
				<!-- <div style="background-color:rgb(250,250,250); box-shadow: 0px 0px 5px 5px rgb(220,220,220); border-radius:5px; padding:13px 15px 15px 25px;"> -->
				<div style="background-color: #77d197; border-top: 4px solid rgb(130,130,130); border-radius: 0 0 5px 5px; padding:13px 15px 15px 25px;">

					<h4>Project: '''+ project +'''</h4>
					Project status: ''' + status +'''<br>
					Project size: ''' + folder_size +'''<br><br>
					<div class="managing-buttons">
						<a href="view-log.htm?p='''+ project +'''" class="button">View log file</a>		
			'''

			if status == 'finished':
				print '''<a href="view-report.htm?p=''' + project + '''" class="button">View report</a> '''
				
			if status == 'running':
				print '''<form><input type="button" class="button" onclick="displayStopProject(\''''+ project +'''\')" value="Stop execution" /></form>'''

			if status != 'running':
				print '''<form><input type="button" class="button" onclick="displayRemoveProject(\''''+ project +'''\')" value="Remove from disk" /></form>'''

			print '</div>'

			print '''
				<div id="stopPrj_'''+ project +'''" style="display:none">
					Are you sure you want to stop running this project?
					<form>
						<input type="button" onclick="stopProject(\''''+ project +'''\')" value="Confirm"/>
						<input type="button" onclick="hideStopProject(\''''+ project +'''\')" value="Cancel"/>
					</form>
				</div>

				<div id="rmPrj_'''+ project +'''" style="display:none;">
					Are you sure you want to permanently remove this project from disk? (Your input data will not be removed).
					<form>
						<input type="button" onclick="removeProject(\''''+ project +'''\')" value="Confirm"/>
						<input type="button" onclick="hideRemoveProject(\''''+ project +'''\')" value="Cancel"/>
					</form>
				</div>
			'''

		else:
			print '''
				<div style="background-color:#79e59d; border:solid green 1px; border-radius:5px; padding:10px;">

				<h4>Project: '''+ project +'''</h4>

				Project status: The status file of this project could not be found on the server. Please delete this project and rerun it.<br>
				<div class="managing-buttons">
				<form><input type="button" class="button" onclick="displayRemoveProject(\''''+ project +'''\')" value="Remove from disk" /></form>
				</div>

				<div id="rmPrj_'''+ project +'''" style="display:none;">
					Are you sure you want to permanently remove this project from disk? (Your input data will not be removed).
					<form>
						<input type="button" onclick="removeProject(\''''+ project +'''\')" value="Confirm"/>
						<input type="button" onclick="hideRemoveProject(\''''+ project +'''\')" value="Cancel"/>
					</form>
				</div>
			'''
		
		print '</div><br>'


