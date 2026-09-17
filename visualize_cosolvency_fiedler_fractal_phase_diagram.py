"""
visualize_cosolvency_fiedler_fractal_phase_diagram.py
=============================================================================
【二元混合溶媒 (Water + Cosolvent) 混合比率 χ に対する (λ2, Df) 相判定スキャン】

【研究目的】
共溶媒モル分率 χ (0.0 ~ 1.0) を変化させた際、溶質-水-共溶媒 ネットワークの 3D 空間グラフから
1) Fiedler 代数的接続度 λ2 (結合統合度)
2) 3D グラフ空間幾何 フラクタル次元 Df (パーコレーション次元)
を同時追跡し、構造相転移を定量描画する。
=============================================================================
"""

import sys
import os
import networkx as nx
import numpy as np
import pandas as pd
from rdkit import Chem

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import fiedler_binary_mixed_solvent_engine as bse

def compute_graph_spatial_fractal_dimension(G):
    if G.number_of_nodes() < 4 or not nx.is_connected(G):
        return 1.0
        
    spl = dict(nx.all_pairs_shortest_path_length(G))
    nodes = list(G.nodes())
    
    r_vals = [1, 2, 3, 4, 5]
    n_avg = []
    
    for r in r_vals:
        counts = []
        for u in nodes:
            c = sum(1 for v, d in spl[u].items() if d <= r)
            counts.append(c)
        n_avg.append(np.mean(counts))
        
    log_r = np.log(r_vals)
    log_n = np.log(n_avg)
    
    slope, _ = np.polyfit(log_r, log_n, 1)
    return float(slope)

def scan_full_cosolvency_phase_trajectory(name, smiles, cosolv_key="ethanol"):
    cosolv_name = bse.COSOLVENT_DATABASE[cosolv_key]["name"]
    print(f"\n=========================================================================================")
    print(f" 🧪 【{name}】における {cosolv_name} 混合比率 χ (0% ➔ 100%) (λ2, Df) 構造相転移プロット")
    print("=========================================================================================")
    
    chis = np.linspace(0.0, 1.0, 11)
    records = []
    
    for chi in chis:
        G = bse.build_binary_mixed_solvent_graph(smiles, cosolv_key=cosolv_key, chi_cosolv=chi, total_solvent_mols=24)
        
        if nx.is_connected(G):
            L_norm = nx.normalized_laplacian_matrix(G).toarray()
            evals = np.linalg.eigvalsh(L_norm)
            f_val = float(evals[1])
        else:
            f_val = 0.0
            
        df_val = compute_graph_spatial_fractal_dimension(G)
        
        # ネットワーク構造相状態の物理的判定
        if f_val >= 0.05 and df_val >= 1.92:
            phase_state = "🌟 [最適高溶解相] 3D バルク均一コソルベント溶媒和 (Homogeneous Cosolvent Phase)"
        elif f_val > 0.03 and df_val >= 1.85:
            phase_state = "📦 [可溶溶媒和網] 2D/3D 水和ケージ・ネットワーク相 (Solvated Network)"
        else:
            phase_state = "❄️ [難溶・相分離相] 疎水反発・1D パーコレーション不連続 (Phase Segregated)"
            
        records.append({
            "chi": chi,
            "pct": int(chi * 100),
            "dielectric": bse.compute_mixed_dielectric(chi, bse.COSOLVENT_DATABASE[cosolv_key]["dielectric"]),
            "fiedler_lambda2": f_val,
            "fractal_df": df_val,
            "phase": phase_state
        })
        
    df_traj = pd.DataFrame(records)
    
    print("-" * 125)
    print(f"{'混合比 χ':<10} | {'共溶媒%':<8} | {'混合誘電率 ε':<12} | {'Fiedler λ2 (結合力)':<22} | {'フラクタル次元 Df (空間幾何)':<26} | ネットワーク構造相")
    print("-" * 125)
    for _, r in df_traj.iterrows():
        print(f"χ = {r['chi']:<6.2f} | {r['pct']:<4}%   | ε = {r['dielectric']:<8.1f} | λ2 = {r['fiedler_lambda2']:<18.6f} | Df = {r['fractal_df']:<20.4f} | {r['phase']}")
    print("-" * 125)
    
    peak_row = df_traj.loc[df_traj["fiedler_lambda2"].idxmax()]
    print(f"\n🏆 【物理的発見】: {name} の最高溶解度達成ポイント")
    print(f"  ・共溶媒最適混合比率 χ_opt : {peak_row['chi']:.2f} ({peak_row['pct']}%)")
    print(f"  ・最大 Fiedler 接続度 λ2   : {peak_row['fiedler_lambda2']:.6f}")
    print(f"  ・到達フラクタル次元 Df    : {peak_row['fractal_df']:.4f} ({peak_row['phase']})")

    out_md_path = f"/home/eldenring/waterMain/COSOLVENCY_PHASE_DIAGRAM_{cosolv_key.upper()}_REPORT.md"
    with open(out_md_path, "w", encoding="utf-8") as f:
        f.write(f"# 🧪 {name} における {cosolv_name} 混合比率 $\\chi$ による $(\\lambda_2, D_f)$ 構造相転移レポート\n\n")
        f.write(df_traj.to_markdown(index=False))
    print(f"🎉 レポート作成完了: {out_md_path}")
    return df_traj

if __name__ == "__main__":
    scan_full_cosolvency_phase_trajectory("安息香酸 (Benzoic Acid)", "O=C(O)c1ccccc1", cosolv_key="ethanol")
    scan_full_cosolvency_phase_trajectory("安息香酸 (Benzoic Acid)", "O=C(O)c1ccccc1", cosolv_key="dmso")
