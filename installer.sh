#!/bin/bash

# === install mono, ThermoRawFileParser and percolator as root === #

sudo apt -y update && sudo apt -y install mono-devel unzip

wget https://github.com/percolator/percolator/releases/download/rel-3-06-01/percolator-v3-06-linux-amd64.deb -O /tmp/percolator.deb
sudo dpkg -i /tmp/percolator.deb
rm /tmp/percolator.deb

ZIP=ThermoRawFileParser1.4.2.zip
wget https://github.com/compomics/ThermoRawFileParser/releases/download/v1.4.2/$ZIP -O /tmp/$ZIP
sudo mkdir -p /opt/compomics
yes | sudo unzip /tmp/$ZIP -d /opt/compomics
rm /tmp/$ZIP

# ============== user level conda setup installation ============= #

INSTALLED_MINICONDA=0
# Check if Anaconda is installed
if command -v conda &> /dev/null; then
    CONDA_INSTALL_DIR=$(conda info --base)
else
    CONDA_INSTALL_DIR="$HOME/miniconda3/"
    echo "Miniconda not found. Installing Miniconda..."
    wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh
    bash /tmp/miniconda.sh -b -p "$CONDA_INSTALL_DIR"
    rm /tmp/miniconda.sh
    conda init
    INSTALLED_MINICONDA=1
fi
# Add Miniconda/Anaconda binary directory to PATH if not already present
if [[ ":$PATH:" != *":$CONDA_INSTALL_DIR/bin:"* ]]; then
    export PATH="$CONDA_INSTALL_DIR/bin:$PATH"
fi

# Configure conda in script
. "$CONDA_INSTALL_DIR/etc/profile.d/conda.sh"
# Activate base environment

conda activate base

# Check if the desired environment exists, create if not
if ! conda env list | grep -q "oktoberfest"; then
    echo "Creating conda environment..."
    conda create -y -n oktoberfest python==3.10
else
    echo "Conda environment already exists. Skipping creation..."
fi

conda activate oktoberfest

# Check if pip installation succeeded
if [ $? -ne 0 ]; then
    echo "Activating the envionment failed."
    exit 1
fi

pip install oktoberfest jupyterlab

# Check if pip installation succeeded
if [ $? -ne 0 ]; then
    echo "Pip installation failed."
    exit 1
fi
# =========================== summary ============================ #

perc_version=`percolator --help 2>&1 | head -n 1`
if [ $? -ne 0 ]; then
    echo "Percolator installation failed. Check logs."
    exit 1
fi

mono_version=`mono --version | head -n 1`
if [ $? -ne 0 ]; then
    echo "Mono installation failed. Check logs."
    exit 1
fi

thermo_version=`mono /opt/compomics/ThermoRawFileParser.exe --version`
if [ $? -ne 0 ]; then
    echo "ThermoRawFileParser installation failed. Check logs."
    exit 1
fi

oktoberfest_version=`pip show oktoberfest | grep Version`
if [ $? -ne 0 ]; then
    echo "Oktoberfest installation failed. Check logs."
    exit 1
fi

echo -e "\n"
echo "percolator version: $perc_version"
echo "mono version: $mono_version"
echo "ThermoRawFileParser version: $thermo_version"
echo "Oktoberfest version: $oktoberfest_version"

echo -e "\nInstallation complete."
if [ $INSTALLED_MINICONDA -eq 1 ]; then
    echo "Miniconda was installed. You need to restart your shell once before you can use oktoberfest."
fi
echo -e "\nTo use oktoberfest, perform the following steps:"
echo -e "\t 1. activate your environment, using \"conda activate oktoberfest\""
echo -e "\t 2. start jupyterlab using \"jupyter lab\" or run oktoberfest directly. Type \"python -m oktoberfest\" to get info."


exit 0
