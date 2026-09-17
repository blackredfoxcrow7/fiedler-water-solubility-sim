# 🧬 Program Architecture & Usage Guide (プログラム全構成ガイド)

This document provides a comprehensive overview of the core engines, pharmaceutical benchmark scripts, and 3D WebGL viewers included in the `fiedler-water-solubility-sim` repository.

*Note: A full Japanese translation is provided in the second half of this document. (後半に日本語訳を併記しています。)*

---

## 🇺🇸 English Guide

### 1. Core Computational Engines

These three primary engines implement the theoretical framework published in our paper (Zenodo DOI: [`10.5281/zenodo.22804666`](https://doi.org/10.5281/zenodo.22804666)).

* **3D Water Lone-Pair Tetrahedral Geometry Engine**:
  * File: [`fiedler_tetrahedral_lonepair_water_engine.py`](fiedler_tetrahedral_lonepair_water_engine.py)
  * Paper Section: Section 2.1
  * Description: Integrates water's tetrahedral lone-pair acceptors/donors ($109.5^\circ$) into a normalized graph Laplacian $\mathbf{L}_{\text{norm}}$, quantifying hydration shell ordering and algebraic connectivity ($\lambda_2$).
* **Grid-Free Analytical Orientation & Dimer Alignment Engine**:
  * File: [`fiedler_dimer_interaction_optimizer.py`](fiedler_dimer_interaction_optimizer.py)
  * Paper Section: Section 2.3
  * Description: Identifies topological complementary interaction sites via Fiedler vector signs ($\mathbf{v}_{2,i} \cdot \mathbf{v}_{2,j} < 0$) and computes optimal 3D dimer/solvent orientation analytically in a single step using Kabsch Singular Value Decomposition (SVD).
* **Multi-Solvent 8-Solvent Profile Engine**:
  * File: [`multi_solvent_fiedler_engine.py`](multi_solvent_fiedler_engine.py)
  * Paper Section: Section 2.2 & 3.1
  * Description: Batch computes Fiedler indices ($\lambda_2$) and solubility profiles across 8 benchmark organic solvents (Water, Ethanol, DMSO, Acetone, Methanol, 1-Butanol, 1-Hexanol, Diethyl Ether).

---

### 2. Pharmaceutical Benchmark & PoC Scripts

* **Pharmaceutical Drug Discovery PoC Script**:
  * File: [`run_pharmaceutical_drug_discovery_poc.py`](run_pharmaceutical_drug_discovery_poc.py)
  * Paper Table: Table 1 (Section 3.1)
  * Benchmarks: Carbamazepine (CBZ), Nicotinamide (NICO), D-Mannitol, Terephthalic Acid.
  * Description: Disentangles micro-hydration attraction ($\lambda_{2, \text{solute-water}}$) from solid-state crystal lattice packing energy ($\Delta f_{\text{dimer}}$) to compute the Fiedler Crystallization Anomaly Index (FCAI).
* **Inverse Spectral Tomography: Molecular Reconstruction**:
  * File: [`fiedler_inverse_solubility_structure_reconstruction.py`](fiedler_inverse_solubility_structure_reconstruction.py)
  * Paper Section: Section 3.2
  * Description: Reconstructs an unknown solute's 3D functional groups and molecular architecture in < 1.2 seconds from an 8-solvent solubility fingerprint vector ($\mathbf{y}_{\text{exp}} \in \mathbb{R}^8$).
* **Cosolvency Phase Diagram Analysis**:
  * File: [`run_mixed_solvent_fractal_fiedler_analysis.py`](run_mixed_solvent_fractal_fiedler_analysis.py)
  * Paper Section: Section 3.3
  * Description: Analyzes non-linear solubility peaks and clathrate-like percolation transitions in mixed solvent systems (Water-DMSO, Water-Ethanol).

---

### 3. Interactive 3D WebGL Viewers

* [`active_fractal_3d_glycerin.html`](active_fractal_3d_glycerin.html): 3D hydration network viewer for Glycerin
* [`real_space_3d_1_butanol.html`](real_space_3d_1_butanol.html): 3D hydration viewer for 1-Butanol
* [`real_space_3d_terephthalic.html`](real_space_3d_terephthalic.html): 3D dimer crystallization viewer for Terephthalic Acid

---

## 🇯🇵 日本語ガイド (Japanese Guide)

### 1. 核心計算エンジン

* **3D水分子Lone-Pair四面体幾何エンジン**: [`fiedler_tetrahedral_lonepair_water_engine.py`](fiedler_tetrahedral_lonepair_water_engine.py)
  * 水分子の孤立電子対（$109.5^\circ$）の3D方向性を規格化ラプラシアン $\mathbf{L}_{\text{norm}}$ に組み込み、水和シェル秩序化と Fiedler値（$\lambda_2$）を算定（論文第2.1節）。
* **一発解析アライメント・二量体最適化エンジン**: [`fiedler_dimer_interaction_optimizer.py`](fiedler_dimer_interaction_optimizer.py)
  * Fiedler固有ベクトルの位相符号（$\mathbf{v}_{2,i} \cdot \mathbf{v}_{2,j} < 0$）から最良相互作用部位を特定し、Kabsch SVD法で最良配向を一発アナリティック算出（論文第2.3節）。
* **8大溶媒Fiedlerプロファイル一括計算エンジン**: [`multi_solvent_fiedler_engine.py`](multi_solvent_fiedler_engine.py)
  * 8溶媒中でのFiedler指標および溶解度プロファイルを一括評価（論文第2.2/3.1節）。

### 2. 論文ベンチマーク検証スクリプト

* **創薬医薬品PoC検証スクリプト**: [`run_pharmaceutical_drug_discovery_poc.py`](run_pharmaceutical_drug_discovery_poc.py)
  * カルバマゼピン（CBZ）、ニコチンアミド（NICO）、D-マンニトール、テレフタル酸の FCAI 指標を算出（論文表1）。
* **逆問題：8溶媒溶解度プロファイルからの構造逆推定**: [`fiedler_inverse_solubility_structure_reconstruction.py`](fiedler_inverse_solubility_structure_reconstruction.py)
  * 8溶媒指紋ベクトル $\mathbf{y}_{\text{exp}}$ から、1.2秒未満で未知分子構造を逆特定（論文第3.2節）。
* **混合溶媒（共溶媒）相図解析**: [`run_mixed_solvent_fractal_fiedler_analysis.py`](run_mixed_solvent_fractal_fiedler_analysis.py)
  * Water-DMSO, Water-Ethanol 混合溶媒における非線形パーコレーションピークを解析（論文第3.3節）。

---

## 💻 Execution Example (実行手順)

```bash
# Run Pharmaceutical Solubility Benchmark (Table 1)
python3 run_pharmaceutical_drug_discovery_poc.py

# Run Multi-Solvent Profile
python3 example_multi_solvent_usage.py
```
