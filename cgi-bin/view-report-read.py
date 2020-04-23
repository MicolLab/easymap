#!src/Python-2.7.18/.localpython/bin/python2

import cgi, cgitb
cgitb.enable()
 
arguments = cgi.FieldStorage()
projectName = str(arguments['p'].value).strip()

print "Content-type:text/html\r\n\r\n"

with open('user_projects/' + projectName + '/3_workflow_output/rfeed.html') as reportFile:
	for line in reportFile:
		print line
