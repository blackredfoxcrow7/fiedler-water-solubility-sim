"""
run_percolation_phase_interface_fiedler_analysis.py
=============================================================================
【パーコレーション グラフ理論 ➔ マクロ相転移・界面生成 (Interface & Phase Boundary) 解析】

【理論背景】
1. パーコレーション転移 (p >= p_c): 巨大連結成分 C_giant の出現 == マクロな相 (結晶/ゲル) の誕生
2. Fiedler 値 λ2: 新相内部の代数的統合度・結晶格子凝集力 (Cohesive Energy)
3. チーガーカット h(G) = min |∂S| / vol(S): グラフ境界カット == 物理的な「新しい界面 (Interface Area)」
=============================================================================
"""

import sys
import os
import networkx as nx
import numpy as np
import pandas as pd

def simulate_percolation_phase_transition(num_molecules=120, p_steps=12):
    print("=========================================================================")
    print(" 🔬 溶液からの相出現 (結晶化/ゲル化) ＆ グラフパーコレーション・界面解析")
    print("=========================================================================")
    
    p_vals = np.linspace(0.05, 0.45, p_steps)
    records = []
    
    np.random.seed(42)
    coords = np.random.randn(num_molecules, 3) * 4.0
    
    for p in p_vals:
        G = nx.Graph()
        for i in range(num_molecules):
            G.add_node(i, pos=coords[i])
            
        for i in range(num_molecules):
            for j in range(i + 1, num_molecules):
                dist = np.linalg.norm(coords[i] - coords[j])
                if dist < 4.5:
                    if np.random.rand() < p:
                        G.add_edge(i, j, weight=1.0 / dist)
                        
        components = list(nx.connected_components(G))
        largest_comp_size = max(len(c) for c in components) if components else 0
        ratio_giant = largest_comp_size / num_molecules
        
        if largest_comp_size > 3:
            G_giant = G.subgraph(max(components, key=len)).copy()
            L_norm = nx.normalized_laplacian_matrix(G_giant).toarray()
            evals = np.linalg.eigvalsh(L_norm)
            fiedler_val = float(evals[1]) if len(evals) > 1 else 0.0
            
            # 界面積を表すグラフカット |∂S|
            try:
                cut_val, _ = nx.stoer_wagner(G_giant)
                interface_area = float(cut_val)
            except:
                interface_area = 0.0
        else:
            fiedler_val = 0.0
            interface_area = 0.0
            
        # 相状態の判定
        if ratio_giant >= 0.70:
            phase_state = "💎 [マクロ相誕生] 完全結晶化 / 全系架橋ゲル相 (Solid / Gel Phase)"
        elif ratio_giant >= 0.30:
            phase_state = "🧪 [核生成・臨界相] ナノ結晶核 / 界面形成 (Nucleation / Interface)"
        else:
            phase_state = "💧 [均一溶液相] 独立モノマー / 自由溶媒和 (Homogeneous Solution)"
            
        records.append({
            "prob_p": p,
            "giant_ratio": ratio_giant,
            "fiedler_lambda2": fiedler_val,
            "interface_cut": interface_area,
            "phase": phase_state
        })
        
    df_res = pd.DataFrame(records)
    
    print("-" * 115)
    print(f"{'結合確率 p':<10} | {'巨大成分割合 %':<14} | {'Fiedler λ2 (凝集力)':<22} | {'界面カット数 |∂S| (界面積)':<24} | 相状態判定")
    print("-" * 115)
    for _, r in df_res.iterrows():
        print(f"p = {r['prob_p']:<7.3f} | {r['giant_ratio']*100:<12.1f}% | λ2 = {r['fiedler_lambda2']:<18.6f} | |∂S| = {r['interface_cut']:<18.2f} | {r['phase']}")
    print("-" * 115)
    
    out_md_path = "/home/eldenring/waterMain/PERCOLATION_PHASE_INTERFACE_REPORT.md"
    with open(out_md_path, "w", encoding="utf-8") as f:
        f.write("# 🔬 溶液からのマクロ相出現 (結晶化/ゲル化) ＆ グラフパーコレーション・界面論証レポート\n\n")
        f.write(df_res.to_markdown(index=False))
    print(f"\n🎉 相・界面理論レポート作成完了: {out_md_path}")
    return df_res

if __name__ == "__main__":
    simulate_percolation_phase_transition()
