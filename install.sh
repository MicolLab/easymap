#!/bin/bash

################################################################################
#
# This script automates some the steps required after cloning or downloading
# in order to make easymap ready for execution.
#
################################################################################

################################################################################
#
# REQUIREMENTS:
#	
#	- To use Easymap through the command line
#		- ...
#
#	- To use Easymap through the web interface
#		- Web server that runs PHP
#
################################################################################

# Deal with argument provided by user
if ! [ $1 ]; then
	echo 'Please provide an argument specifying the type of installation: "cli" or "server". Example: "./install.sh server"'
	exit
fi

if ! [ $1 == server ] && ! [ $1 == cli ]; then
	echo 'Please choose between "cli" and "server". Example: "./install.sh server"'
	exit
fi

if [ $1 == server ]; then
	if ! [ $2 ]; then
		port=8100
	elif [ "$2" -ge 8100 ] && [ "$2" -le 8200 ]; then
		port=$2
	else
		echo 'Please choose a port number between 8100 and 8200. Example: "./install.sh server 8100"'
		exit
	fi
fi

################################################################################

# Create some folders not present in GitHub repo (e.g. 'user_data' and 'user_projects')

[ -d user_data ] || mkdir user_data
[ -d user_projects ] || mkdir user_projects

[ -d web_interface/tmp_upload_files ] || mkdir web_interface/tmp_upload_files

################################################################################

# Compile bcftools, bowtie, hisat and samtools

cd ./htslib
make clean
make

cd ../bcftools-1.3.1 
make clean
make

cd ../bowtie2 
make clean
make

cd ../samtools1 
make clean
make

cd ../hisat2
make clean
make

cd ..

################################################################################

# Create src directory 
if [ -d src ]; then rm -rf src; fi
mkdir src
cd src

# Get Python-2.7.18
wget https://www.python.org/ftp/python/2.7.18/Python-2.7.18.tgz
tar -zxvf Python-2.7.18.tgz

# Install Python-2.7.18
cd Python-2.7.18
[ -d .localpython ] || mkdir .localpython
./configure --prefix=$PWD/.localpython
make
make install
cd ..

# Get virtualenv-15.1.0
wget https://pypi.python.org/packages/d4/0c/9840c08189e030873387a73b90ada981885010dd9aea134d6de30cd24cb8/virtualenv-15.1.0.tar.gz#md5=44e19f4134906fe2d75124427dc9b716
tar -zxvf virtualenv-15.1.0.tar.gz

# Install virtualenv-15.1.0
cd virtualenv-15.1.0/
../Python-2.7.18/.localpython/bin/python setup.py install

# Create virtual environment "easymap-env"
../Python-2.7.18/.localpython/bin/python virtualenv.py easymap-env -p ../Python-2.7.18/.localpython/bin/python

# Install Pillow with pip
[ -d cache ] || mkdir cache
easymap-env/bin/pip -qq install Pillow --cache-dir cache

cd ../..

################################################################################

# Change permissions to the easymap folder and subfolders so Easymap can be used both from the
# web interface (server user -- e.g. www-data) and the command line of any user
sudo chmod -R 777 .

# In file 'easymap', set absolute path to the Python binaries of the virtual environment
# Rest of Python scripts don't need this because are executed after easymap.sh activates the virtual environment
#sed -i -e "s~ABS_PATH_ENV_PYTHON~${PWD}/src/Python-2.7.18/.localpython/bin/python2~g" easymap

################################################################################

# Check if Easymap functions properly by running a small project: 
cp fonts/check.1.fa user_data/
cp fonts/check.gff user_data/
run_result=`./easymap -n setup -w snp -sim -r check -g check.gff -ed ref_bc_parmut`

# Cleanup
rm  user_data/check.gff
rm  user_data/check.1.fa 
rm -rf user_projects/*

if [ "$run_result" == "Easymap analysis properly completed." ]; then

	# Set easymap dedicated  HTTP CGI server to run always in the background
	if [ $1 == server ]; then
		
		# Run server in the background
		nohup ./src/Python-2.7.18/.localpython/bin/python2 -m CGIHTTPServer $port &
		
		# Modify/create the etc/crontab file to always start easymap server at bootup
		echo "@reboot   root    cd $PWD; ./src/Python-2.7.18/.localpython/bin/python2 -m CGIHTTPServer $port" >> /etc/crontab

		# Save port number to /config/port for future reference for the user
		echo $port > config/port

	fi

	echo " "
	echo " "
	echo "###################################################################################"
	echo "#                                                                                 #"
	echo "#                                                                                 #"
	echo "#                   Easymap installation successfully completed                   #"
	echo "#                                                                                 #"
	echo "#                                                                                 #"
	echo "###################################################################################"
	echo " "
	echo " "

else

	echo " "
	echo " "
	echo "###################################################################################"
	echo "#                                                                                 #"
	echo "#                                                                                 #"
	echo "#                          Easymap installation failed                            #"
	echo "#                                                                                 #"
	echo "#                                                                                 #"
	echo "###################################################################################"
	echo " "
	echo " "

fi
