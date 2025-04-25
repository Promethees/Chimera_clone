#!/bin/bash

# setup.sh: Automate bioinformatics environment setup on Ubuntu 24.04 LTS

# Exit on error
set -e

# Check Ubuntu version
echo "Checking Ubuntu version..."
UBUNTU_VERSION=$(lsb_release -rs)
if [[ "$UBUNTU_VERSION" != "24.04" ]]; then
    echo "Warning: This script is optimized for Ubuntu 24.04. You have $UBUNTU_VERSION."
fi

# Update system
echo "Updating system packages..."
sudo apt update
sudo apt upgrade -y
sudo apt install -y build-essential git wget curl unzip libgl1-mesa-glx libegl1

# Install Miniconda
if ! command -v conda &> /dev/null; then
    echo "Installing Miniconda..."
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
    bash miniconda.sh -b -p $HOME/miniconda3
    rm miniconda.sh
    export PATH="$HOME/miniconda3/bin:$PATH"
    conda init bash
    source ~/.bashrc
    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge
    conda update -y conda
else
    echo "Miniconda already installed."
fi

# Create and activate Conda environment
if ! conda env list | grep -q "bioinformatics"; then
    echo "Creating Conda environment 'bioinformatics'..."
    conda create -y -n bioinformatics python=3.9
else
    echo "Conda environment 'bioinformatics' already exists."
fi
source $HOME/miniconda3/bin/activate bioinformatics

# Install Python libraries
echo "Installing Python libraries..."
conda install -y numpy pandas scipy matplotlib seaborn jupyter
pip install biopython transformers torch deepchem

# Install bioinformatics tools
echo "Installing bioinformatics tools..."
conda install -y hmmer muscle pymol-open-source snakemake
sudo apt install -y autodock-vina
conda install -y gromacs

# Install Rosetta (requires user-provided binary)
if [ ! -d "/opt/rosetta" ]; then
    echo "Rosetta installation requires an academic license."
    echo "Please download rosetta_bin_linux_2024.01_bundle.tgz from https://www.rosettacommons.org/software/license-and-download"
    read -p "Enter the full path to rosetta_bin_linux_2024.01_bundle.tgz: " ROSETTA_TAR
    if [ -f "$ROSETTA_TAR" ]; then
        tar -xvzf "$ROSETTA_TAR"
        sudo mv rosetta_bin_linux_2024.01 /opt/rosetta
        echo 'export PATH=/opt/rosetta/main/source/bin:$PATH' >> ~/.bashrc
        source ~/.bashrc
    else
        echo "Error: Rosetta tarball not found. Skipping Rosetta installation."
    fi
else
    echo "Rosetta already installed at /opt/rosetta."
fi

# Optional: Install Docker
if ! command -v docker &> /dev/null; then
    read -p "Install Docker for reproducible environments? (y/n): " INSTALL_DOCKER
    if [ "$INSTALL_DOCKER" = "y" ]; then
        echo "Installing Docker..."
        sudo apt install -y docker.io
        sudo usermod -aG docker $USER
        echo "Docker installed. Log out and back in to apply group changes."
    fi
else
    echo "Docker already installed."
fi

# Optional: GPU setup
if lspci | grep -qi nvidia; then
    read -p "NVIDIA GPU detected. Install CUDA for GPU support? (y/n): " INSTALL_CUDA
    if [ "$INSTALL_CUDA" = "y" ]; then
        echo "Installing CUDA..."
        sudo apt install -y nvidia-driver-550 nvidia-cuda-toolkit
        pip install torch --index-url https://download.pytorch.org/whl/cu121
    fi
else
    echo "No NVIDIA GPU detected. Skipping CUDA installation."
fi

# Test installation
echo "Testing installations..."
python -c "import Bio; print('Biopython:', Bio.__version__)"
if command -v pymol &> /dev/null; then
    echo "PyMOL installed."
else
    echo "PyMOL installation may have failed."
fi
if command -v rosetta_scripts &> /dev/null; then
    echo "Rosetta installed."
else
    echo "Rosetta installation may have failed."
fi

echo "Setup complete! Activate environment with 'conda activate bioinformatics'."
echo "See README.md for usage and troubleshooting."