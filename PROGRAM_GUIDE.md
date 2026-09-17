# 🧬 分子溶解度・水和・構造逆推定シミュレーター：プログラム全構成ガイド

本ドキュメントは、本リポジトリ（`fiedler-water-solubility-sim`）に含まれるコアエンジン、創薬PoC検証スクリプト、3D WebGLビューアの役割と使用方法を体系的にまとめたガイドです。

---

## 1. 核心計算エンジン（Core Engines）

論文（DOI: `10.5281/zenodo.22804666`）で発表された数理理論を直接実行する 3 大コアエンジンです。

### ① 3D水分子Lone-Pair四面体幾何エンジン
* **ファイル**: [`fiedler_tetrahedral_lonepair_water_engine.py`](fiedler_tetrahedral_lonepair_water_engine.py)
* **対応論文セクション**: 第2.1節
* **機能**: 水分子の孤立電子対（Lone Pair: $109.5^\circ$）の3D指向性を規格化グラフ・ラプラシアン $\mathbf{L}_{\text{norm}}$ に組み込み、水和シェルの秩序化および Fiedler値（$\lambda_2$）を算定。

### ② 一発解析アライメント・二量体最適化エンジン
* **ファイル**: [`fiedler_dimer_interaction_optimizer.py`](fiedler_dimer_interaction_optimizer.py)
* **対応論文セクション**: 第2.3節
* **機能**: Fiedler固有ベクトルの位相符号（$\mathbf{v}_{2,i} \cdot \mathbf{v}_{2,j} < 0$）から最良相互作用部位を特定し、Kabsch SVD（特異値分解）を用いて格子探索なしで最良配向を一発アナリティック算出。

### ③ 8大溶媒Fiedlerプロファイル一括計算エンジン
* **ファイル**: [`multi_solvent_fiedler_engine.py`](multi_solvent_fiedler_engine.py)
* **対応論文セクション**: 第2.2節・第3.1節
* **機能**: 水、エタノール、DMSO、アセトン、メタノール、1-ブタノール、1-ヘキサノール、ジエチルエーテルの8溶媒中でのFiedler指標および溶解度プロファイルを一括評価。

---

## 2. 論文ベンチマーク検証スクリプト

論文内の実験比較表や逆問題解析を再現・実行するためのスクリプト群です。

### ① 創薬医薬品PoC検証スクリプト
* **ファイル**: [`run_pharmaceutical_drug_discovery_poc.py`](run_pharmaceutical_drug_discovery_poc.py)
* **対応論文セクション**: 第3.1節（表1）
* **計算対象**: カルバマゼピン（CBZ）、ニコチンアミド（NICO）、D-マンニトール、テレフタル酸
* **内容**: 単分子水和 Fiedler値と二量体自己会合利得の解離から「結晶化アノマリー指標 (FCAI)」を算出し、難溶性メカニズムを判定。

### ② 逆問題：8溶媒溶解度プロファイルからの構造逆推定
* **ファイル**: [`fiedler_inverse_solubility_structure_reconstruction.py`](fiedler_inverse_solubility_structure_reconstruction.py)
* **対応論文セクション**: 第3.2節
* **内容**: 8溶媒の溶解度指紋ベクトル $\mathbf{y}_{\text{exp}} \in \mathbb{R}^8$ から、1.2秒未満で未知分子の骨格および官能基分布を逆特定。

### ③ 混合溶媒（共溶媒）相図解析
* **ファイル**: [`run_mixed_solvent_fractal_fiedler_analysis.py`](run_mixed_solvent_fractal_fiedler_analysis.py)
* **対応論文セクション**: 第3.3節
* **内容**: Water-DMSO, Water-Ethanol 混合溶媒における非線形溶解度ピーク（パーコレーション開通）を解析。

---

## 3. 3D WebGL インタラクティブビューア

ブラウザで開くだけで、分子と水和水素結合ネットワークを3D表示するツールです。

* [`active_fractal_3d_glycerin.html`](active_fractal_3d_glycerin.html): グリセリンの3D水和ネットワーク表示
* [`real_space_3d_1_butanol.html`](real_space_3d_1_butanol.html): 1-ブタノールの水和3D表示
* [`real_space_3d_terephthalic.html`](real_space_3d_terephthalic.html): テレフタル酸の二量体・晶析3D表示

---

## 💻 実行手順の例

### 例1: 創薬医薬品ベンチマーク（表1）の再計算
```bash
python3 run_pharmaceutical_drug_discovery_poc.py
```

### 例2: 8溶媒溶解度プロファイル一括計算
```bash
python3 example_multi_solvent_usage.py
```
