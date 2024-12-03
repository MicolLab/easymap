#!./src/Python-3.12.3/.localpython/bin/python3


import cgi, cgitb, os
cgitb.enable() # For debugging only

# Set Content-type header so XMLHttpRequest in JS can understand the response
print("Content-Type: text/html")
print("")

# Get the f argument from the URL
arguments = cgi.FieldStorage()
fileName = str(arguments['f'].value).strip()

# Remove the file stored in f
os.remove('./user_data/' + fileName)