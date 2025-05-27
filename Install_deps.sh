#!/bin/bash
# ProteinFoldPro Dependency Installer

# Verify Conda installation
if ! command -v conda &> /dev/null
then
    echo "Conda not found! Install Miniconda first."
    exit 1
fi

# Create environment
conda create -n proteinfold -y python=3.10
conda activate proteinfold

# Core packages
conda install -y -c conda-forge -c bioconda \
    colabfold=1.5.2 \
    hhsuite=3.3.0 \
    openmm=8.0 \
    pdbfixer=1.8.1

# Python dependencies
pip install \
    biopython==1.81 \
    fair-esm==2.0.0 \
    requests==2.31.0 \
    pandas==2.0.3 \
    matplotlib==3.7.1

echo "Installation complete! Activate with 'conda activate proteinfold'"
