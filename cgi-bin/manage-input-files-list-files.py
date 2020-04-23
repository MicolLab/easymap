#!./src/Python-2.7.18/.localpython/bin/python2

import os, cgi, cgitb, math
cgitb.enable()

dirname = './user_data'

# Get list of files
filepaths = []
for basename in os.listdir(dirname):
    filename = os.path.join(dirname, basename)
    if os.path.isfile(filename):
        filepaths.append(filename)

# Re-populate list with filename, size tuples
for i in xrange(len(filepaths)):
    filepaths[i] = (filepaths[i], os.path.getsize(filepaths[i]))

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
print ""

if len(filepaths)==0:
		print '<div class="emptyResult">There are no input files currently on the server.</div>'

else:
	print '<hr class="files-separator"></hr>'
	for f in sorted(filepaths):
		user_input_file = str(f[0].split('user_data/')[1])
		file_size = convert_size(int(f[1]))
		if user_input_file[-3:] != '.gz':
			print '''
					<div id="file_'''+ user_input_file +'''">
						<div class="files-item left">
							<h4>'''+ user_input_file +'''</h4>
						</div>
						<div class="files-item center">
							<h4>'''+ file_size +'''</h4>
						</div>
						<div class="files-item right1">
							<a href="../cgi-bin/manage-input-files-retrieve-header.py?f='''+ user_input_file +'''" target="_blank" class="button" style="width:100px">Preview</a>
						</div>
						<div class="files-item right2">
							<form>
								<input style="margin: 3 0 7px 0" type="button" class="button" onclick="ShowRemoveFile(\'''' + user_input_file + '''\')" value="Remove from disk" />
							</form>
						</div>
						<div id="rmFile_'''+ user_input_file +'''" style="display:none; clear:both;">
							Are you sure you want to permanently remove this file from disk?
							<form>
								<input type="button" onclick="removeFile(\''''+ user_input_file +'''\')" value="Confirm" />
								<input type="button" onclick="HideRemoveFile(\''''+ user_input_file +'''\')" value="Cancel" />
							</form>
						</div>
						<hr class="files-separator"></hr>
					</div>
			'''
