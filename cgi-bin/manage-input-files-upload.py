#!src/Python-2.7.18/.localpython/bin/python2

import cgi, os, shutil, subprocess
import cgitb; cgitb.enable()


print 'Content-Type: text/html'
print ''
print '{"jsonrpc" : "2.0", "result" : null, "id" : "id"}'

'''
try: # Windows needs stdio set for binary mode.
    import msvcrt
    msvcrt.setmode (0, os.O_BINARY) # stdin  = 0
    msvcrt.setmode (1, os.O_BINARY) # stdout = 1
except ImportError:
    #pass
    with open('out.txt', 'w') as out:
		out.write('module missing')
'''

form = cgi.FieldStorage()

# Generator to buffer file chunks
def fbuffer(f, chunk_size=10000):
	while True:
		chunk = f.read(chunk_size)
		if not chunk:
			break
		yield chunk

# A nested FieldStorage instance holds the file
fileItem = form['file']
fileName = form['name'].value
currentChunk = form['chunk'].value
totalChunks = str(form['chunks'].value)

# For debugging only
#testOut = open('out.txt', 'a')

# Test if the chunk was uploaded and append it to server temporary file
if fileItem.filename:
	# Strip leading path from file name to avoid directory traversal attacks
	fn = os.path.basename(fileItem.filename)
	f = open('web_interface/tmp_upload_files/' + fn, 'ab', 10000)

	for chunk in fbuffer(fileItem.file):
		f.write(chunk)
	f.close()
	message = ''

else:
	message = ''

#testOut.write(currentChunk + ', ' + totalChunks + ', ')

# If chunk is the last to complete a file, process the server temporal file
if currentChunk == str(int(totalChunks) - 1):

	# If file is gzipped...
	if fileName[-3:] == '.gz':

		'''
		Files are received from the client. If files are smaller than the chunk size set in the client side
		(manage-input-files.js), they arive as a single chunk and they retain their original name. If they
		are over the chunk size, they are divided in chunks and the name of the reconstructed file in the
		server ib 'blob'. The following Try/Except structures make the script work in both situations, even
		if I change the size o the chunks in manage-input-files.js. The chunk size is originally 50 mb.
		'''
		
		# If file was chunked, its name will be 'blob', so rename it to original name
		if int(totalChunks) > 1:
			os.rename('web_interface/tmp_upload_files/blob','web_interface/tmp_upload_files/' + fileName)

		# Move file to user_data dir (.gz files are not displayed to user by manage-input-files-list-files.py)
		shutil.move('web_interface/tmp_upload_files/' + fileName, 'user_data/' + fileName)
		
		# Once in user_data dir, gunzip the file (user can see the gunzip progress)
		subprocess.Popen('gzip -d user_data/' + fileName, shell=True)
	
	# Equivalent but simpler for non .gz files...
	else:
		if int(totalChunks) > 1:
			os.rename('web_interface/tmp_upload_files/blob','web_interface/tmp_upload_files/' + fileName)
		
		shutil.move('web_interface/tmp_upload_files/' + fileName, 'user_data/' + fileName)


#testOut.close()