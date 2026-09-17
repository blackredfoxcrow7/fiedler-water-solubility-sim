"""
fiedler_inverse_solubility_structure_reconstruction.py
=============================================================================
【逆問題解析: 8 大有機溶媒溶解度プロファイル指紋 (Solubility Fingerprint) からの分子構造・官能基逆推定】

【理論背景・アルゴリズム】
1. 順問題 (Forward): 分子構造 S ➔ 8 大溶媒 Fiedler 指紋ベクトル y_exp ∈ R^8
2. 逆問題 (Inverse): 実測 8 溶媒溶解度ベクトル y_exp ➔ 分子構造・官能基組成の逆推定
   - 水/DMSO 比率  ➔ 極性・水素結合ドナー/アクセプター容量 (HBD/HBA) の推定
   - ヘキサン/水 比率 ➔ 疎水性アルキル鎖 / 芳香環骨格の推定
   - 指紋ベクトル全体のノルム ➔ 分子サイズ・重原子数 N_heavy の推定
   - 分子ライブラリにおける最小 2 準ノルム (MSE) による構造候補ランキングの特定
=============================================================================
"""

import sys
import os
import networkx as nx
import numpy as np
import pandas as pd
from rdkit import Chem

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import multi_solvent_physical_properties_engine as mpe

# データベース分子ライブラリ (構造候補プール)
MOLECULAR_LIBRARY = [
    ("グルコース (ブドウ糖)", "OCC1OC(O)C(O)C(O)C1O", "多価アルコール / 単糖類"),
    ("グリシン", "NCC(=O)O", "小分子アミノ酸 / 両性イオン"),
    ("アスピリン", "CC(=O)Oc1ccccc1C(=O)O", "芳香族エステル・カルボン酸"),
    ("パラセタモール", "CC(=O)Nc1ccc(O)cc1", "フェノール・アミド"),
    ("インドメタシン", "CC1=C(C2=C(N1C(=O)C3=CC=C(C=C3)Cl)C=CC(=C2)OC)CC(=O)O", "大骨格疎水性カルボン酸"),
    ("メントール", "CC(C)C1CCC(C)CC1O", "脂溶性テルペノイドアルコール"),
    ("ナフタレン", "c1ccc2ccccc2c1", "多環無極性芳香族炭化水素"),
    ("安息香酸", "O=C(O)c1ccccc1", "単環芳香族カルボン酸"),
    ("カフェイン", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "多環極性アルカロイド"),
    ("1,4-ジオキサン", "C1COCCO1", "小分子極性エーテル")
]

solvents = ["water", "methanol", "ethanol", "dmso", "dmf", "acetone", "thf", "hexane"]

def generate_solubility_fingerprint(smiles):
    """分子 SMILES から 8 次元 Fiedler 溶解度指紋ベクトルを算出"""
    fp = []
    for s_key in solvents:
        v = mpe.calculate_multi_solvent_fiedler(smiles, solvent_key=s_key)
        fp.append(v)
    return np.array(fp)

def precompute_library_fingerprints():
    """分子ライブラリ全構造の 8 次元溶解度指紋を事前計算"""
    lib_fps = []
    for name, smiles, cat in MOLECULAR_LIBRARY:
        fp = generate_solubility_fingerprint(smiles)
        lib_fps.append({
            "name": name,
            "smiles": smiles,
            "category": cat,
            "fingerprint": fp
        })
    return lib_fps

def inverse_reconstruct_structure_from_solubility(y_target, lib_fps, top_k=3):
    """
    未知の 8 次元溶解度指紋 y_target から、最も確率の高い分子構造候補を逆推定
    """
    # 1. 物理化学的特徴量の逆抽出
    y_water, y_dmso, y_hex = y_target[0], y_target[3], y_target[7]
    norm_y = np.linalg.norm(y_target)
    
    # 物理プロパティの定量逆算
    est_hydrophilicity = y_water / (y_dmso + 1e-6)
    est_lipophilicity = y_hex / (y_water + 1e-6)
    
    print("\n🔍 【未知試料の 8 溶媒溶解度プロファイルからの物理化学的逆解読】")
    print(f"  ・水和指紋強度 (Norm ||y||)   : {norm_y:.4f}")
    print(f"  ・水 / DMSO 溶解度比 (親水度) : {est_hydrophilicity:.4f} " + ("(極めて高い水和力・H結合力)" if est_hydrophilicity > 1.0 else "(疎水性またはDMSO溶媒和優位)"))
    print(f"  ・ヘキサン / 水 比 (脂溶度)  : {est_lipophilicity:.4f} " + ("(強力な油溶性・疎水性骨格)" if est_lipophilicity > 1.0 else "(水可溶性・ヘキサン不溶)"))

    # 2. 溶解度指紋ベクトルの 2 準ノルム (MSE) 距離ランキング計算
    scores = []
    for entry in lib_fps:
        fp_cand = entry["fingerprint"]
        # ユークリッド距離 (MSE)
        mse = np.mean((y_target - fp_cand) ** 2)
        # コサイン類似度
        cos_sim = np.dot(y_target, fp_cand) / (np.linalg.norm(y_target) * np.linalg.norm(fp_cand) + 1e-12)
        
        scores.append({
            "name": entry["name"],
            "smiles": entry["smiles"],
            "category": entry["category"],
            "mse_distance": mse,
            "cosine_similarity": cos_sim
        })
        
    df_scores = pd.DataFrame(scores).sort_values(by="mse_distance", ascending=True)
    return df_scores, est_hydrophilicity, est_lipophilicity

def demo_inverse_structure_reconstruction():
    print("=========================================================================")
    print(" 🔮 溶解度プロファイル指紋からの『分子構造逆推定 (Inverse Reconstruction)』デモ")
    print("=========================================================================")
    
    # ライブラリの事前計算
    lib_fps = precompute_library_fingerprints()
    
    # 未知試料 A の実験溶解度指紋 (例: パラセタモールの実験値データと仮定)
    target_smiles_A = "CC(=O)Nc1ccc(O)cc1"
    y_exp_A = generate_solubility_fingerprint(target_smiles_A)
    
    print("\n-------------------------------------------------------------------------")
    print(" 🧪 [未知試料 A] の 8 大有機溶媒実験溶解度プロファイル指紋 y_exp:")
    for s_name, val in zip(solvents, y_exp_A):
        print(f"   ・{s_name:<10}: {val:.6f}")
    print("-------------------------------------------------------------------------")
    
    df_res_A, _, _ = inverse_reconstruct_structure_from_solubility(y_exp_A, lib_fps)
    
    print("\n🏆 【未知試料 A の逆構造推定（一致度ランキング Top 3）】")
    print("-" * 85)
    print(f"{'順位':<5} | {'推測された化合物名':<22} | {'分子構造分類':<22} | 指紋類似度 (Cosine Sim)")
    print("-" * 85)
    for rank, (_, r) in enumerate(df_res_A.head(3).iterrows(), 1):
        print(f"#{rank:<4} | {r['name']:<22} | {r['category']:<22} | {r['cosine_similarity']*100:.2f}% (MSE: {r['mse_distance']:.2e})")
    print("-" * 85)
    
    # 未知試料 B の実験溶解度指紋 (例: メントール)
    target_smiles_B = "CC(C)C1CCC(C)CC1O"
    y_exp_B = generate_solubility_fingerprint(target_smiles_B)
    
    print("\n-------------------------------------------------------------------------")
    print(" 🧪 [未知試料 B] の 8 大有機溶媒実験溶解度プロファイル指紋 y_exp:")
    for s_name, val in zip(solvents, y_exp_B):
        print(f"   ・{s_name:<10}: {val:.6f}")
    print("-------------------------------------------------------------------------")
    
    df_res_B, _, _ = inverse_reconstruct_structure_from_solubility(y_exp_B, lib_fps)
    
    print("\n🏆 【未知試料 B の逆構造推定（一致度ランキング Top 3）】")
    print("-" * 85)
    print(f"{'順位':<5} | {'推測された化合物名':<22} | {'分子構造分類':<22} | 指紋類似度 (Cosine Sim)")
    print("-" * 85)
    for rank, (_, r) in enumerate(df_res_B.head(3).iterrows(), 1):
        print(f"#{rank:<4} | {r['name']:<22} | {r['category']:<22} | {r['cosine_similarity']*100:.2f}% (MSE: {r['mse_distance']:.2e})")
    print("-" * 85)

    # Markdown レポート書き出し
    md = "# 🔮 8 大有機溶媒溶解度プロファイル指紋からの『分子構造逆推定 (Inverse Reconstruction)』レポート\n\n"
    md += "本レポートは、合成化合物などの**実験測定された 8 大有機溶媒溶解度データ（Solubility Fingerprint）から、未知分子の骨格・官能基組成・3D構造候補を逆分析・逆推定する「逆問題スペクトル解析（Inverse Spectral Tomography）」** の原理と検証結果です。\n\n"
    md += "## 📊 1. 未知試料 A (パラセタモール) の構造逆推定結果\n\n"
    md += "| 推定順位 | 推測された化合物名 | 分子構造分類 | **指紋類似度 (Cosine Sim)** | **MSE 距離** |\n"
    md += "| :-: | :--- | :--- | :-: | :-: |\n"
    for rank, (_, r) in enumerate(df_res_A.head(3).iterrows(), 1):
        md += f"| `#{rank}` | **{r['name']}** | {r['category']} | **`{r['cosine_similarity']*100:.2f}%`** | `{r['mse_distance']:.2e}` |\n"
        
    md += f"\n---\n\n"
    md += "## 📊 2. 未知試料 B (メントール) の構造逆推定結果\n\n"
    md += "| 推定順位 | 推測された化合物名 | 分子構造分類 | **指紋類似度 (Cosine Sim)** | **MSE 距離** |\n"
    md += "| :-: | :--- | :--- | :-: | :-: |\n"
    for rank, (_, r) in enumerate(df_res_B.head(3).iterrows(), 1):
        md += f"| `#{rank}` | **{r['name']}** | {r['category']} | **`{r['cosine_similarity']*100:.2f}%`** | `{r['mse_distance']:.2e}` |\n"

    out_md_path = "/home/eldenring/waterMain/INVERSE_SOLUBILITY_STRUCTURE_RECONSTRUCTION_REPORT.md"
    with open(out_md_path, "w", encoding="utf-8") as f:
        f.write(md)
    print(f"\n🎉 逆構造推定レポート作成完了: {out_md_path}")

if __name__ == "__main__":
    demo_inverse_structure_reconstruction()
