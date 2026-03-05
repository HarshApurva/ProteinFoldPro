# ProteinFold Pro 🧬🔬
*Advanced Multi-Method Protein Structure Prediction Pipeline*

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)]()

Predict protein 3D structures using state-of-the-art AI methods with automated quality checks.

## Features
- 🚀 ColabFold (AlphaFold2) integration
- ⚡ ESMFold for fast predictions
- 🔍 HHblits & BLAST MSA generation
- 🛠️ OpenMM structure refinement
- 📊 Confidence metrics (pLDDT, PAE)
- 🐳 Docker container support

## Quick Install
```bash
conda create -n proteinfold python=3.10
conda activate proteinfold
bash scripts/install_deps.sh
