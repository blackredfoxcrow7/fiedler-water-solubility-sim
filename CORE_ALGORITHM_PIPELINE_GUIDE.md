# 📘 分子式 (SMILES) 入力から結果出力までの統合アルゴリズムパイプライン解説書

本ドキュメントは、**「SMILES（分子構造式）の入力」** から始まり、**「ラプラシアン Fiedler 固有値の演算」「3D 座標の幾何生成」「Kabsch SVD アライメント」「結晶・溶解度解析」** に至るコアエンジンの内部アルゴリズムとプログラム構造を詳細かつ分かりやすくまとめた仕様書です。

---

## ⚙️ 1. 全体パイプライン概要

プログラム全体は、以下の 7 つのステップで自動処理されます：

```
 [入力: SMILES 文字列] (例: "NC(Cc1ccc(O)cc1)C(=O)O")
          │
          ▼  【Step 1】分子グラフ構築 & 水素付加
     (RDKit: Chem.AddHs)
          │
          ▼  【Step 2】実空間 3D デカルト座標の力場生成 (Å 単位)
     (RDKit: EmbedMolecule + MMFF94)
          │
          ▼  【Step 3】分子間相互作用（H結合/塩橋/π-π/疎水）候補エッジの抽出
     (自動属性判別)
          │
          ▼  【Step 4】正規化ラプラシアン行列 L_norm の構築
     (L_norm = I - D^{-1/2} A D^{-1/2})
          │
          ▼  【Step 5】スペクトル分析（Fiedler値 λ2 & 固有ベクトル v の算出）
     (Eigen Solver: np.linalg.eigvalsh)
          │
          ▼  【Step 6】Fiedler 最大化 (max λ2) によるトポロジー選択
     (貪欲アルゴリズム)
          │
          ▼  【Step 7】Kabsch SVD 最適化による実空間 3D 最短アライメント
     (SVD 分解: H = U Σ V^T ➔ R, t)
          │
          ▼
 [出力: 実効溶解度 ΔΔf / 結晶化利得 Δf / 1:1:1 実空間 3D HTML 分子模型]
```

---

## 📄 2. 主要関数とその詳細解説

### ① `create_integrated_system(smiles_list)`
* **役割**: 入力された SMILES 文字列を解釈し、原子ノードと共有結合エッジを持つグラフ構造を構築する。
* **内部処理**:
  1. `Chem.MolFromSmiles(smiles)` で分子構造をパース。
  2. `Chem.AddHs(mol)` により、極性・水素結合に不可欠な明示的水素（`-H`）を全て付加。
  3. 各原子に元素記号、形式電荷、水素結合可能数（`max_h_bonds`）、芳香属性（`is_aromatic`）の属性タグを注記。

---

### ② `get_3d_conformer(smiles)`
* **役割**: 近似グラフレイアウトではなく、実際の物理化学法則に基づく**実空間 3D デカルト座標 $\mathbf{X} \in \mathbb{R}^{N \times 3}$（オングストローム $\text{Å}$ 単位）** を生成する。
* **内部処理**:
  1. `AllChem.EmbedMolecule(mol)`: 距離幾何学（Distance Geometry）により 3D 骨格を初期発生。
  2. `AllChem.MMFFOptimizeMolecule(mol)`: MMFF94 分子力場でエネルギー最小化を実行（C-C 結合長 $1.54\,\text{Å}$、C=O 結合長 $1.22\,\text{Å}$、結合角 $109.5^\circ / 120^\circ$ 等の真の3D空間幾何が確立）。

---

### ③ 分子間相互作用エッジの自動識別モジュール
* **役割**: 2つの分子（A と B）の全原子対 $(A_i, B_j)$ に対し、化学的相補性に基づいて分子間結合エッジと重み $w_{ij}$ を自動割り当てする。
* **判定ルール**:
  * **塩橋 / イオン結合**（カルボキシル O $\leftrightarrow$ アミノ N）: **重み $w = 2.0$**
  * **水素結合**（ドナー O-H, N-H $\leftrightarrow$ アクセプター O, N）: **重み $w = 1.5$**
  * **芳香環 $\pi-\pi$ スタッキング**（芳香族 C $\leftrightarrow$ 芳香族 C）: **重み $w = 1.2$**
  * **疎水結合**（非極性 C $\leftrightarrow$ 非極性 C）: **重み $w = 0.7$**

---

### ④ `compute_fiedler_normalized(G)`
* **役割**: 構築された分子・水グラフ $G$ の**正規化ラプラシアン行列（Normalized Laplacian）**を作成し、第2固有値（Fiedler値 $\lambda_2$）を解出する。
* **数理演算**:
  1. 隣接行列 $\mathbf{A}$（要素 $A_{ij} = w_{ij}$）と次数行列 $\mathbf{D}$（対角要素 $D_{ii} = \sum_j A_{ij}$）を作成。
  2. 正規化ラプラシアン行列を計算：
     $$\mathbf{L}_{\text{norm}} = \mathbf{I} - \mathbf{D}^{-1/2} \mathbf{A} \mathbf{D}^{-1/2}$$
  3. エルミート行列の固有値ソルバー（`np.linalg.eigvalsh`）を実行し、昇順ソートした固有値群 $0 = \lambda_0 \le \lambda_1 \le \lambda_2 \dots$ から、グラフの代数的接続度を表す **$\lambda_1$ または $\lambda_2$（Fiedler値）** を抽出。

---

### ⑤ Fiedler 最大化アルゴリズム（トポロジー最適化）
* **役割**: 分子間候補エッジの中から、トポロジー全体の Fiedler 値 $\lambda_2$ が**最も高くなる（構造が最も安定化する）エッジ網を貪欲選択（Greedy Selection）**する。
* **内部処理**:
  1. 候補エッジを重みの高い順にソート。
  2. 1原子あたりの最大分子間エッジ数（例: 最大2本）の物理制約を守りつつ、エッジを1本ずつ追加。
  3. 追加後の Fiedler 値 $\lambda_{2,\text{test}}$ が従来の最高値を超える場合のみエッジを確定採用。

---

### ⑥ `kabsch_rigid_transform(P, Q, weights)`
* **役割**: Fiedler 最大化で選ばれた分子間結合ペアに対し、**分子 B を実空間 3D 空間内で最適に回転・平行移動させ、分子間結合距離を物理的最小に収束**させる。
* **数理構造（Kabsch SVD アルゴリズム）**:
  1. 重心移動: $P_c = P - \bar{P}, \quad Q_c = Q - \bar{Q}$
  2. 加重クロス共分散行列の計算:
     $$\mathbf{H} = \mathbf{Q}_c^T \mathbf{W} \mathbf{P}_c$$
  3. 特異値分解（SVD）:
     $$\mathbf{H} = \mathbf{U} \mathbf{\Sigma} \mathbf{V}^T$$
  4. 最適回転行列 $\mathbf{R}$ と平行移動ベクトル $\mathbf{t}$ の確定:
     $$\mathbf{R} = \mathbf{V} \mathbf{U}^T, \quad \mathbf{t} = \bar{P} - \bar{Q} \mathbf{R}^T$$
  5. 座標変換: $\mathbf{X}_B' = \mathbf{X}_B \mathbf{R}^T + \mathbf{t}$（**平均分子間距離が $1.39\,\text{Å} \sim 3.70\,\text{Å}$ の物理極小へ一発収束**）。

---

### ⑦ 代数的位相ベクトル（$\mathbf{d}_{\text{phase}}$）モジュール
* **役割**: Fiedler 固有ベクトルのスカラー値 $v(i)$ と 3D 座標 $\mathbf{x}_i$ から、**「分子の親水/疎水・電位の向きを示す空間軸ベクトル」** を計算する。
* **数理演算**:
  $$\mathbf{d}_{\text{phase}} = \sum_{i=1}^N v(i) \cdot (\mathbf{x}_i - \bar{\mathbf{x}})$$
  * 異種分子やドッキング計算において、双方の $\mathbf{d}_{\text{phase}}$ を**逆平行（180°対向）**に向き合わせることで、6次元空間探索なしでベストな 3D 配向を一発指定。

---

## 📊 3. コアエンジンの統一性まとめ

本パイプラインは、入力される SMILES の種類や構成を変えるだけで、同一のコードベース上で以下の全ての物理量を算出します：

* **単分子水和入力**: 単分子水和力 $\Delta f_{\text{solvated}}$ ➔ 二段階実効溶解度 $\Delta \Delta f$
* **同種2分子入力 (`[A, A]`)**: ホモダイマー Fiedler $\lambda_{2,\text{homo}}$ ➔ 結晶化自己会合利得 $\Delta f_{\text{homo}}$
* **異種2分子入力 (`[A, B]`)**: ヘテロダイマー Fiedler $\lambda_{2,\text{hetero}}$ ➔ 共結晶・ドッキング利得 $\Delta f_{\text{hetero}}$

共通のプログラムファイル: [`/home/eldenring/waterMain/fiedler_dimer_interaction_optimizer.py`](file:///home/eldenring/waterMain/fiedler_dimer_interaction_optimizer.py)
