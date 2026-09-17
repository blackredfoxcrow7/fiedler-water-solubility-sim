# 📘 Fiedler トポロジーシミュレータ 最新版プログラム使用説明書 (Usage Guide)

本ドキュメントは、これまでに開発・完成された**「水素結合官能基動的設定」**、**「無次元空間適応型 衝突回避」**、**「背景構造オフセット差分算出」**、および**「積極的フラクタル水網成長」**を完備した最新版シミュレータ群の構成および使用方法をまとめたマニュアルです。

---

## 📁 1. ファイル構成一覧

すべての主要プログラムは `/home/eldenring/waterMain/` ディレクトリ内に配置されています。

| 役割 | ファイル名 | 主な機能・補正仕様 |
| :--- | :--- | :--- |
| **🏆 決定版メインシミュレータ** | [`fiedler_sim_hbond_offset_clash_free.py`](file:///home/eldenring/waterMain/fiedler_sim_hbond_offset_clash_free.py) | **・衝突回避（重なり防止）内蔵**<br>**・H数に応じた次数動的判定**<br>**・背景オフセット差分（$\Delta f$）自動算出** |
| **🌀 強フラクタル水網シミュレータ** | [`fiedler_fractal_water_sim.py`](file:///home/eldenring/waterMain/fiedler_fractal_water_sim.py) | **・パワーロー結合 ($r^{-\alpha}$)**<br>**・優先的選択付着 (カイレイ木構造)**<br>**・自己相似的フラクタル水網の成長** |
| **🎨 3モード3D可視化ツール** | [`visualize_toggle_water.py`](file:///home/eldenring/waterMain/visualize_toggle_water.py) | **・「🌐全表示」「🧪溶質のみ」「💧水網のみ」のワンクリック切替ボタン付き3D HTML生成** |
| **🧪 衝突ゼロ検証スクリプト** | [`verify_clash_free.py`](file:///home/eldenring/waterMain/verify_clash_free.py) | **・原子めり込みゼロの自動テスト** |
| **📊 実験・物性検証スクリプト** | [`run_proposal_step1.py`](file:///home/eldenring/waterMain/run_proposal_step1.py) | **・ポリオール類・環状vs直鎖分子の物性検証** |
| **⚖️ モデル比較スクリプト** | [`compare_fractal_sim.py`](file:///home/eldenring/waterMain/compare_fractal_sim.py) | **・標準モデル vs フラクタルモデル比較** |

---

## 💻 2. 決定版メインシミュレータの使い方

`fiedler_sim_hbond_offset_clash_free.py` をPythonスクリプトから呼び出し、分子のSMILES表記から水溶性評価値（実効変化量 $\Delta f_{\text{final}}$）を算出する基本コードです。

```python
import sys
import os

# 1. パスの設定とモジュールのインポート
sys.path.append("/home/eldenring/waterMain")
import fiedler_sim_hbond_offset_clash_free as fsim

# 2. 評価したい分子の SMILES を指定して初期システム（溶質＋水分子網）を構築
smiles = "CC(C)(C)O"  # 例: tert-ブタノール
G_init, w_os = fsim.create_integrated_system([smiles], water_per_atom=3, is_loop=True)

# 3. 衝突回避＆次数動的制限付きシミュレーションを実行（100ステップ）
f_history, frames_data, delta_f_final, delta_f_history = fsim.run_fiedler_opt_sim(
    G_init, w_os, steps=100
)

# 4. 評価結果の取得
print(f"初期 Fiedler 値 (f_init)   : {f_history[0]:.6f}")
print(f"最終 Fiedler 値 (f_final)  : {f_history[-1]:.6f}")
print(f"実効水和変化量 (Δf_final)  : {delta_f_final:.6f}")  # f_final - f_init
```

---

## 🎨 3. 3モード切替 3D 可視化ツールの使い方

`visualize_toggle_water.py` を使用して、画面上のボタンで「水分子あり」「溶質のみ」「水分子のみ」をワンクリックで切り替えられる 3D Plotly HTML を作成します。

```python
import sys
sys.path.append("/home/eldenring/waterMain")
import visualize_toggle_water as vis

# 任意の分子の SMILES と分子名を入力して 3D HTML を出力
html_path = vis.generate_toggle_3d_visualization("OCC(O)CO", "glycerin", steps=100)

print(f"生成された 3D HTML ファイル: {html_path}")
```

### 💡 ブラウザ画面での操作機能
生成された `.html` ファイルをウェブブラウザで開くと、左上に以下の**3切替ボタン**が配置されています：
1. **`🌐 すべて表示 (Show All)`**: 溶質・水分子・共有結合・水素結合・水網を全表示。
2. **`🧪 溶質のみ (Solute Only)`**: 水分子を非表示にし、水網から押し出された**溶質自体の立体構造**のみを視察。
3. **`💧 水分子のみ (Water Only)`**: 溶質を非表示にし、**水分子同士が形成する水網（シアン色のクラスタ）**のみを観察。

---

## 🌀 4. 強フラクタル水網シミュレータの使い方

水分子同士を積極的に自己相似な樹状フラクタル構造（パワーロー結合）へと成長させたい場合は `fiedler_fractal_water_sim.py` を使用します。

```python
import sys
sys.path.append("/home/eldenring/waterMain")
import fiedler_fractal_water_sim as fsim_frac

# 1. 初期システム構築
G_init, w_os = fsim_frac.fsim_base.create_integrated_system(["CCCO"], water_per_atom=3, is_loop=True)

# 2. パワーロー指数 alpha_powerlow を指定して強フラクタル化シミュレーション実行
f_history, frames_data, delta_f_final = fsim_frac.run_fiedler_fractal_opt_sim(
    G_init, w_os, steps=100, alpha_powerlow=2.0
)

print(f"強フラクタル水網 Δf_final: {delta_f_final:.6f}")
```

---

## 🚀 5. コマンドライン（ターミナル）からの直接実行方法

各種検証テストスクリプトは、Linux ターミナルから以下のコマンドで即座に実行可能です。

```bash
# A. 原子の衝突・めり込みがゼロ（0）であることの検証
python3 /home/eldenring/waterMain/verify_clash_free.py

# B. ポリオール類（グリセリン等）や環状分子（THF, シクロヘキサノール）の物性シミュレーション実行
python3 /home/eldenring/waterMain/run_proposal_step1.py

# C. 標準シミュレータ vs 強フラクタルシミュレータの比較実行
python3 /home/eldenring/waterMain/compare_fractal_sim.py
```

---

## ⚙️ 6. 主要な物理・数理パラメーター仕様まとめ

* **衝突回避閾値 (Clash Avoidance Thresholds)**:
  * 重原子同士: `0.065` 未満で斥力作用
  * 水素を含むペア: `0.045` 未満で斥力作用
* **水素結合能の動的判定ルール ($n_H$)**:
  * 酸素 (O): $n_H \ge 1 \implies \text{max\_h\_bonds}=2, \text{num\_waters}=2$ (ドナー＋アクセプター)
  * 酸素 (O): $n_H = 0 \implies \text{max\_h\_bonds}=1, \text{num\_waters}=1$ (アクセプターのみ)
  * 窒素 (N): $n_H \ge 1 \implies \text{max\_h\_bonds}=n_H+1, \text{num\_waters}=n_H+1$
  * 窒素 (N): $n_H = 0 \implies \text{max\_h\_bonds}=1, \text{num\_waters}=1$
* **実効水和変化量 (Effective Solvation Delta)**:
  * $\Delta f_{\text{final}} = f(T) - f(0)$ （構造背景ノイズ・中心性バイアスの完全消去）
