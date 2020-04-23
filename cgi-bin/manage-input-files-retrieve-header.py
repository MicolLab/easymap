#!src/Python-2.7.18/.localpython/bin/python2

import cgi, cgitb
cgitb.enable()

print "Content-type:text/html\r\n\r\n" 
arguments = cgi.FieldStorage()
nbrOfLines = 1000
fileName = str(arguments['f'].value).strip()
c=0


# Header
print(
	"""
	<p>File: """ + fileName + """</p>
	<p>(showing first """ +str(nbrOfLines)+ """ lines)</p>
	<p><pre>***********************************************************</pre></p>
	"""
	)

# File
with open('user_data/' + fileName) as fileName:
	print('<pre>')
	for line in fileName:
		#print line.strip()
		print cgi.escape(line.strip()).encode('ascii', 'xmlcharrefreplace')
		c+=1
		if c == nbrOfLines: break
	print('</pre>')
