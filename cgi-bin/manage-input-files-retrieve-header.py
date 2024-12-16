#!src/Python-3.12.3/.localpython/bin/python3

import cgi, cgitb, html
cgitb.enable()

print("Content-type:text/html\r\n\r\n") 
arguments = cgi.FieldStorage()
nbrOfLines = 1000
fileName = str(arguments['f'].value).strip()
c=0


# Header
print((
	"""
	<p>File: """ + fileName + """</p>
	<p>(showing first """ +str(nbrOfLines)+ """ lines)</p>
	<p><pre>***********************************************************</pre></p>
	"""
	))

# File
with open('user_data/' + fileName) as fileName:
	print('<pre>')
	for line in fileName:
		#print line.strip()
		print(html.escape(line.strip()))
		c+=1
		if c == nbrOfLines: break
	print('</pre>')
