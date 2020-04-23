#!./src/Python-2.7.18/.localpython/bin/python2

import cgi, cgitb, shutil
cgitb.enable() # For debugging only

# Set Content-type header so XMLHttpRequest in JS can understand the response
print "Content-Type: text/html"
print ""

# Get the p argument from the URL
arguments = cgi.FieldStorage()
directory = str(arguments['p'].value).strip()

# Remove the directory stored in p
shutil.rmtree('./user_projects/' + directory)
#shutil.rmtree('../user_projects/' + directory, ignore_errors=True) # Usefull with read-only files?