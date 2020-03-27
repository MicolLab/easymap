#!/bin/bash

################################################################################
#
# This script perform tasks required to unistall easymap
# 	-Remove starting HTTPCGIServer from system-wide crontab
#
#
################################################################################

# Remove from /etc/crontab the line that executes Easymap dedicated HTTPCGI server at bootup 
sudo sed -i '/easymap/c\' /etc/crontab
