"""
run_water_microcrystal_graph_size_analysis.py
=============================================================================
【水分子微結晶 (Water Micro-Crystal) のグラフサイズ・空間的拡がりの定量評価モジュール】

【評価する 4 大サイズ指標】
1. ノード数 N (水分子数) ＆ 水素結合エッジ数 |E|
2. グラフ直径 D_graph ＆ 平均パス長 L (結晶のネットワークスパン)
3. 代数的体積 vol(G) ＆ 3D 慣性半径 Rg (空間的広がり)
4. Fiedler 接続度 λ2 ＆ 構造密度の収束
=============================================================================
"""

import sys
import os
import networkx as nx
import numpy as np
import pandas as pd

def generate_water_microcrystal_graph(N_water):
    """
    氷 Ih 格子ライクな 4 配向四面体水素結合網を持つ水分子微結晶グラフを生成
    """
    np.random.seed(42 + N_water)
    # 3D 氷格子空間に水分子を配置
    side = int(np.ceil(N_water ** (1/3)))
    coords = []
    for x in range(side):
        for y in range(side):
            for z in range(side):
                if len(coords) < N_water:
                    # 氷の四面体ジグザグ配置の揺らぎ
                    coords.append([x * 2.75 + (y%2)*0.5, y * 2.75 + (z%2)*0.5, z * 2.75])
    coords = np.array(coords)
    
    G = nx.Graph()
    for i in range(N_water):
        G.add_node(i, pos=coords[i])
        
    # 水素結合 (距離 3.2 Å 以内の分子間にエッジ)
    for i in range(N_water):
        for j in range(i + 1, N_water):
            d = np.linalg.norm(coords[i] - coords[j])
            if d <= 3.2:
                G.add_edge(i, j, weight=1.0 / d)
                
    # 最大連結成分の抽出
    if not nx.is_connected(G):
        largest_cc = max(nx.connected_components(G), key=len)
        G = G.subgraph(largest_cc).copy()
        
    return G, coords

def analyze_microcrystal_graph_size():
    print("=========================================================================")
    print(" 🔬 水分子微結晶 (Water Micro-Crystal) のグラフサイズ ＆ 幾何広がり解析")
    print("=========================================================================")
    
    sizes = [4, 8, 16, 32, 64, 128]
    records = []
    
    for N in sizes:
        G, coords = generate_water_microcrystal_graph(N)
        actual_N = G.number_of_nodes()
        actual_E = G.number_of_edges()
        
        # 1. グラフ直径 (Diameter) ＆ 平均パス長
        try:
            diameter = nx.diameter(G)
            avg_path_length = nx.average_shortest_path_length(G)
        except:
            diameter = 0
            avg_path_length = 0.0
            
        # 2. 代数的体積 (Volume)
        vol_G = sum(dict(G.degree()).values())
        
        # 3. 3D 物理的慣性半径 Rg (Radius of Gyration)
        sub_coords = np.array([coords[n] for n in G.nodes()])
        centroid = np.mean(sub_coords, axis=0)
        rg = np.sqrt(np.mean(np.sum((sub_coords - centroid)**2, axis=1)))
        
        # 4. Fiedler 代数接続度 λ2
        L_norm = nx.normalized_laplacian_matrix(G).toarray()
        evals = np.linalg.eigvalsh(L_norm)
        fiedler_val = float(evals[1]) if len(evals) > 1 else 0.0
        
        records.append({
            "size_N": actual_N,
            "edges_E": actual_E,
            "diameter_D": diameter,
            "avg_path_L": avg_path_length,
            "volume_vol": vol_G,
            "radius_Rg": rg,
            "fiedler_lambda2": fiedler_val
        })
        
    df_res = pd.DataFrame(records)
    
    print("-" * 110)
    print(f"{'水分子数 N':<10} | {'水素結合数 |E|':<14} | {'グラフ直径 D':<14} | {'平均パス長 L':<14} | {'慣性半径 Rg (Å)':<16} | Fiedler λ2")
    print("-" * 110)
    for _, r in df_res.iterrows():
        print(f"N = {r['size_N']:<8} | |E| = {r['edges_E']:<10} | D = {r['diameter_D']:<10} | L = {r['avg_path_L']:<10.2f} | Rg = {r['radius_Rg']:<12.2f} Å | λ2 = {r['fiedler_lambda2']:.6f}")
    print("-" * 110)
    
    out_md_path = "/home/eldenring/waterMain/WATER_MICROCRYSTAL_GRAPH_SIZE_REPORT.md"
    with open(out_md_path, "w", encoding="utf-8") as f:
        f.write("# 🔬 水分子微結晶 (Water Micro-Crystal) のグラフサイズ・空間評価レポート\n\n")
        f.write(df_res.to_markdown(index=False))
    print(f"\n🎉 微結晶グラフサイズレポート作成完了: {out_md_path}")
    return df_res

if __name__ == "__main__":
    analyze_microcrystal_graph_size()
