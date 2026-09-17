# 🧪 Fiedler Spectral Physics Simulator for Molecular Solvation & Co-Crystal Design

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![Field: Computational Chemistry](https://img.shields.io/badge/Field-Spectral_Physics-green.svg)]()

> **Non-Iterative Graph Spectral Engine for Solvation Energy, Co-Solvency, and Co-Crystal Alignment**

---

## 🌟 Key Features

1. **Exact Spectral Physics ($\mathcal{O}(N \log N)$)**: Computes Normalized Graph Laplacian ($\mathbf{L}_{\text{norm}}$) spectra, replacing heavy DFT/MD calculations.
2. **3D Water Lone-Pair Orientation ($109.5^\circ$)**: Incorporates directional hydrogen bonding into graph edge weights.
3. **8 Major Organic Solvents**: Supports Water, MeOH, EtOH, DMSO, DMF, Acetone, THF, and Hexane.
4. **Co-Crystal Alignment Without Grid Search**: Uses Fiedler phase complementarity ($\mathbf{v}_i \cdot \mathbf{v}_j < 0$) + Kabsch SVD for instant 3D docking.
5. **Inverse Spectral Tomography**: Reconstructs molecular structures and functional groups from experimental 8-solvent solubility profiles.

---

## 🚀 Quick Start

```bash
git clone https://github.com/yourusername/fiedler-solvation-engine.git
cd fiedler-solvation-engine
pip install rdkit networkx numpy pandas
```

### Run Multi-Solvent Solubility Benchmark
```python
from multi_solvent_physical_properties_engine import calculate_multi_solvent_fiedler

# Calculate Fiedler solubility value for Ibuprofen in Water vs DMSO
f_water = calculate_multi_solvent_fiedler("CC(C)Cc1ccc(cc1)C(C)C(=O)O", "water")
f_dmso  = calculate_multi_solvent_fiedler("CC(C)Cc1ccc(cc1)C(C)C(=O)O", "dmso")

print(f"Ibuprofen Water Fiedler: {f_water:.6f}") # Low solubility
print(f"Ibuprofen DMSO  Fiedler: {f_dmso:.6f}")  # High solvation
```

---

## 📄 Citation & Preprint

If you use this simulator in your research, please cite our preprint:

```bibtex
@article{independent2026fiedler,
  title={Spectral Graph Physics for Molecular Solvation and Co-Crystal Design: A Non-Iterative Fiedler Laplacian Architecture},
  author={Independent Researcher},
  journal={ChemRxiv Preprint},
  year={2026}
}
```
