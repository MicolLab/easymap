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
echo "EMPEZANDO CON EL PROCESO DE INSTALACION"

# Install necessary dependencies for Python
echo "Instalando dependencias necesarias para Python"
sudo apt-get update
sudo apt-get install -y build-essential libssl-dev zlib1g-dev libbz2-dev \
                        libreadline-dev libsqlite3-dev wget curl llvm \
                        libncurses5-dev libncursesw5-dev xz-utils tk-dev \
                        liblzma-dev python3-openssl git python3-dev

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
	echo "MODO SERVIDOR SELECCIONADO"
	if ! [ $2 ]; then
		port=8100
		echo "USANDO EL PUERTO POR DEFECTO 8100"
	elif [ "$2" -ge 8100 ] && [ "$2" -le 8200 ]; then
		port=$2
		echo "USANDO EL PUERTO $port"
	else
		echo 'Please choose a port number between 8100 and 8200. Example: "./install.sh server 8100"'
		exit
	fi
fi

################################################################################
echo "CREANDO LOS DIRECTORIOS REQUERIDOS"

# Create some folders not present in GitHub repo (e.g. 'user_data' and 'user_projects')
[ -d user_data ] || mkdir user_data
[ -d user_projects ] || mkdir user_projects
[ -d web_interface/tmp_upload_files ] || mkdir web_interface/tmp_upload_files

################################################################################
echo "COMPILANDO LOS DIRECTORIOS REQUERIDOS"

# Compile bcftools, bowtie, hisat and samtools

cd ./htslib
echo "COMPILANDO HTSLIB"
make clean
make

cd ../bcftools-1.3.1 
echo "COMPILANDO BCFTOOLS"
make clean
make

cd ../bowtie2 
echo "COMPILANDO BOWTIE2"
make clean
make

cd ../samtools1 
echo "COMPILANDO SAMTOOLS"
make clean
make

cd ../hisat2
echo "COMPILANDO HISAT2"
make clean
make

cd ..

################################################################################
echo "CREANDO EL DIRECTORIO SRC E INSTALANDO PYTHON"

# Create src directory 
if [ -d src ]; then rm -rf src; fi
mkdir src
cd src

# Get Python 3
echo "DESCARGANDO PYTHON"
wget https://www.python.org/ftp/python/3.12.3/Python-3.12.3.tgz
tar -zxvf Python-3.12.3.tgz

# Install Python-3.12.3
cd Python-3.12.3
[ -d .localpython ] || mkdir .localpython
echo "CONFIGURANDO LA INSTALACION DE PYTHON"
./configure --prefix=$PWD/.localpython --enable-optimizations
make
make install
cd ..

# Get virtualenv-20.28.0
echo "DESCARGANDO VIRTUALENV"
wget https://files.pythonhosted.org/packages/bf/75/53316a5a8050069228a2f6d11f32046cfa94fbb6cc3f08703f59b873de2e/virtualenv-20.28.0.tar.gz
tar -zxvf virtualenv-20.28.0.tar.gz

# Install virtualenv-20.28.0
cd virtualenv-20.28.0/
echo "INSTALANDO VIRTUALENV"
../Python-3.12.3/.localpython/bin/python3 -m pip install virtualenv

# Create virtual environment "easymap-env"
echo "CREANDO EL ENTORNO VIRTUAL"
../Python-3.12.3/.localpython/bin/python3 -m virtualenv easymap-env -p ../Python-3.12.3/.localpython/bin/python3

# Ensure pip is up to date
echo "ACTUALIZANDO pip A LA ÚLTIMA VERSION"
easymap-env/bin/python3 -m pip install --upgrade pip

# Ensure pip is installed in the virtual environment
if ! [ -x "$(command -v easymap-env/bin/pip)" ]; then
    echo "Error: pip no está instalado en el entorno virtual, intentando instalar pip..."
    easymap-env/bin/python3 -m ensurepip --upgrade
fi

# Install Pillow with pip
[ -d cache ] || mkdir cache
echo "INSTALANDO PILLOW"
easymap-env/bin/pip -qq install Pillow --cache-dir cache
echo "INSTALACION CORRECTA DE PILLOW"

cd ../..

################################################################################
echo "CAMBIANDO LOS PERMISOS"
# Change permissions to the easymap folder and subfolders so Easymap can be used both from the
# web interface (server user -- e.g. www-data) and the command line of any user
sudo chmod -R 777 .

# In file 'easymap', set absolute path to the Python binaries of the virtual environment
# Rest of Python scripts don't need this because are executed after easymap.sh activates the virtual environment
#sed -i -e "s~ABS_PATH_ENV_PYTHON~${PWD}/src/Python-3.12.3/.localpython/bin/python3~g" easymap

################################################################################
echo "EASYMAP CHECK"
# Check if Easymap functions properly by running a small project: 
# Verify and copy files if they exist
cp fonts/check.1.fa user_data/
cp fonts/check.gff user_data/

echo "COMANDO RUN_RESULT"
run_result=`./easymap -n setup -w snp -sim -r check -g check.gff -ed ref_bc_parmut 2>&1`
echo "COMANDO RUN_RESULT BIEN"

# Cleanup
rm user_data/check.gff
m user_data/check.1.fa
rm -rf user_projects/*

if [ "$run_result" == "Easymap analysis properly completed." ]; then
	echo "LOS ANALISIS DE EASYMAP SE HAN COMPLETADO CORRECTAMENTE"
	# Set easymap dedicated  HTTP CGI server to run always in the background
	if [ $1 == server ]; then
		echo "CORRIENDO EN MODO SERVIDOR"
		# Run server in the background
		nohup ./src/Python-3.12.3/.localpython/bin/python3 -m http.server --cgi $port &
		
		# Modify/create the etc/crontab file to always start easymap server at bootup
		echo "@reboot   root    cd $PWD; ./src/Python-3.12.3/.localpython/bin/python3 -m CGIHTTPServer $port" >> /etc/crontab

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