## Complete User Guide

### 1. Command Line Interface
```bash
# Basic prediction
python -m proteinfoldpro --sequence MKALIV... --output results/

# Advanced options
python -m proteinfoldpro \
  --input my_protein.fasta \
  --use-colabfold \
  --num-models 5 \
  --enable-openmm
```

### 2. Python API
```python
from proteinfoldpro import ProteinFoldPro

# Create predictor with custom options
predictor = ProteinFoldPro(
    output_dir="custom_results",
    enable_gpu=True
)

# Run full prediction workflow
result = predictor.fold_protein(
    sequence="MKALIVLGLVLL...",
    protein_id="my_protein",
    use_all_methods=True
)
```

### 3. Output Files
| File Type          | Description                          |
|--------------------|--------------------------------------|
| `*.pdb`            | Predicted 3D structure               |
| `*_scores.json`    | Confidence metrics                   |
| `*_msa.a3m`        | Multiple sequence alignment          |
| `*_report.md`      | Validation summary                   |
