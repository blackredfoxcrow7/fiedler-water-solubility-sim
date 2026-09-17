# Multi-Solvent Molecular Solubility & Hydration Simulator
### Explainable Spectral Graph Theory via 3D Water Lone-Pair Fiedler Vector Optimization ($\lambda_2$)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22804666.svg)](https://doi.org/10.5281/zenodo.22804666)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![Paper DOI](https://img.shields.io/badge/Paper-Zenodo--DOI-blue.svg)](https://doi.org/10.5281/zenodo.22804666)
[![Author: Yoshihiro Honda](https://img.shields.io/badge/Author-Yoshihiro%20Honda-orange.svg)](#author)

> **Can multi-solvent solubility ($\log S$), co-crystallization, and structural reconstruction be predicted without brute-force Quantum Mechanics (DFT/MD) or black-box machine learning (GNNs)?**  
> **Yes.** This repository presents a novel **Spectral Graph Theory** framework that calculates the **Fiedler value ($\lambda_2$)** with a 3D water lone-pair tetrahedral kernel ($109.5^\circ$) to predict multi-solvent solubility profiles, co-crystal dimer orientations, and inverse structural tomography at $\mathcal{O}(N \log N)$ computational complexity.

---

## 🌟 Key Scientific Innovations

1. **3D Water Lone-Pair Directional Kernel ($109.5^\circ$)**
   * Integrates water's tetrahedral lone-pair acceptors/donors into a normalized graph Laplacian $\mathbf{L}_{\text{norm}} = D^{-1/2} L D^{-1/2}$, quantifying hydration shell ordering.

2. **Grid-Free Analytical Orientation Alignment (Kabsch SVD)**
   * Identifies topological complementary interaction sites via Fiedler vector signs ($\mathbf{v}_{2,i} \cdot \mathbf{v}_{2,j} < 0$) and computes optimal 3D dimer/solvent orientation analytically in a single step using Kabsch Singular Value Decomposition (SVD).

3. **Fiedler Crystallization Anomaly Index (FCAI)**
   * Disentangles micro-hydration attraction ($\lambda_{2, \text{solute-water}}$) from solid-state crystal lattice packing energy ($\Delta f_{\text{dimer}}$), explaining crystallization anomalies (e.g., D-Mannitol).

4. **Inverse Spectral Tomography**
   * Reconstructs an unknown solute's 3D functional groups and molecular architecture in < 1.2 seconds from an 8-solvent solubility fingerprint vector ($\mathbf{y}_{\text{exp}} \in \mathbb{R}^8$).

---

## 📊 Benchmark Multi-Solvent Solubility Results

| Compound | Solvent System | Experimental Solubility ($\log S$) | Fiedler Index ($\lambda_2$) | FCAI Score | Physical Mechanism |
| :--- | :--- | :---: | :---: | :---: | :--- |
| **Carbamazepine (CBZ)** | Water | -3.85 (Low) | 1.12 | -2.45 | Strong Dimerization ($\pi$-stacking) |
| **Carbamazepine (CBZ)** | Ethanol / DMSO | High | 3.48 | +1.20 | Fiedler Percolation Network Open |
| **Nicotinamide (NICO)** | Water | +0.82 (High) | 3.82 | +2.15 | Strong Lone-Pair Hydrogen Bonding |
| **D-Mannitol** | Water | Moderate | 4.12 | -3.80 | High Crystal Packing Penalty (FCAI) |
| **Terephthalic Acid** | Water | Very Low | 0.95 | -4.10 | Rigid Symmetric Crystal Lattice |

---

## 📄 Official Published Manuscripts

* 📄 **[Zenodo Paper Link](https://doi.org/10.5281/zenodo.22804666)**: *"Molecular Solubility, Hydration Shells, and Inverse Spectral Tomography via 3D Lone-Pair Fiedler Vector Optimization"* (DOI: `10.5281/zenodo.22804666`)
* 📄 **English Paper Manuscript**: [paper2_solubility_draft.md](paper2_solubility_draft.md)
* 📄 **Japanese Paper Translation (日本語訳)**: [paper2_solubility_draft_ja.md](paper2_solubility_draft_ja.md)

---

## 🚀 Quick Start & Installation

### Prerequisites
```bash
pip install numpy scipy rdkit matplotlib flask
```

### Running Pharmaceutical Solubility Benchmark
```bash
python3 run_pharmaceutical_drug_discovery_poc.py
```

### Running Multi-Solvent Validation
```bash
python3 run_expanded_multi_solvent_validation_benchmark.py
```

---

## 👤 Author

**Yoshihiro Honda (本多 義弘)**  
Independent Researcher in Organic Chemistry & Computational Chemistry, Japan  
GitHub: [@blackredfoxcrow7](https://github.com/blackredfoxcrow7)  
Publication DOI: [10.5281/zenodo.22804666](https://doi.org/10.5281/zenodo.22804666)  

*Collaborative Research & AI Technical Assistance provided by Antigravity (Google DeepMind).*

---

## 📜 Citation

```bibtex
@article{honda2026solubility,
  title={Molecular Solubility, Hydration Shells, and Inverse Spectral Tomography via 3D Lone-Pair Fiedler Vector Optimization},
  author={Honda, Yoshihiro},
  journal={Zenodo},
  year={2026},
  doi={10.5281/zenodo.22804666},
  url={https://doi.org/10.5281/zenodo.22804666}
}
```
