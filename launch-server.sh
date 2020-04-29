#!/bin/bash
# A simple script to launch the easymap server, takes one argument: a port number between 8100 and 8200

nohup ./src/Python-2.7.18/.localpython/bin/python2 -m CGIHTTPServer $1 &
