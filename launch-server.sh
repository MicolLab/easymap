#!/bin/bash
# A simple script to launch the easymap server, takes one argument: a port number between 8100 and 8200

# Start the server in the brackground
nohup ./src/Python-3.12.3/.localpython/bin/python3 -m http.server --cgi $1 &

# Check if the server started successfully
sleep 2

# Print the link to the user
echo "Easymap server is running at: http://localhost:$1"