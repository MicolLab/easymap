#!src/Python-3.12.3/.localpython/bin/python3

import cgi, cgitb, html
cgitb.enable()
 
arguments = cgi.FieldStorage()
projectName = str(arguments['p'].value).strip()

print("Content-type:text/html\r\n\r\n")

with open('user_projects/' + projectName + '/2_logs/log.log') as logFile:
	for line in logFile:
		print(html.escape(line.strip()) + '<br>')
