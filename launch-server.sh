#!/bin/bash
# A simple script to launch the easymap server, takes one argument: a port number between 8100 and 8200

nohup ./src/Python-3.12.3/.localpython/bin/python3 -m CGIHTTPServer $1 &
