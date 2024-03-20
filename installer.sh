#!/bin/bash

# === install mono, ThermoRawFileParser and percolator as root === #

sudo apt -y update && sudo apt -y install mono-devel

wget https://github.com/percolator/percolator/releases/download/rel-3-06-01/percolator-v3-06-linux-amd64.deb -O /tmp/percolator.deb
sudo dpkg -i /tmp/percolator.deb
rm /tmp/percolator.deb

ZIP=ThermoRawFileParser1.4.2.zip
wget https://github.com/compomics/ThermoRawFileParser/releases/download/v1.4.2/$ZIP -O /tmp/$ZIP
sudo mkdir -p /opt/compomics
yes | sudo unzip /tmp/$ZIP -d /opt/compomics
rm /tmp/$ZIP

# ============== user level conda setup installation ============= #

# Check if Anaconda is installed
if command -v conda &> /dev/null; then
    CONDA_INSTALL_DIR=$(conda info --base)
else
    CONDA_INSTALL_DIR="$HOME/miniconda3/"
    echo "Miniconda not found. Installing Miniconda..."
    wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh
    /tmp/miniconda.sh -b -p "$CONDA_INSTALL_DIR"
    rm /tmp/miniconda.sh
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


echo -e "\n"
echo "percolator version: `percolator --help 2>&1 | head -n 1`"
echo "mono version: `mono --version | head -n 1`"
echo "ThermoRawFileParser version: `mono /opt/compomics/ThermoRawFileParser.exe --version`"
echo "Oktoberfest version: `pip show oktoberfest | grep Version`"

echo -e "\nInstallation complete. To use oktoberfest, perform the following steps:"
echo -e "\t 1. activate your environment, using \"conda activate oktoberfest\""
echo -e "\t 2. start jupyterlab using \"jupyter lab\" or run oktoberfest directly. Type \"python -m oktoberfest\" to get info."


exit 0
