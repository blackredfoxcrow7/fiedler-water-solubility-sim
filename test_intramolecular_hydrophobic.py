import sys
import os
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import fiedler_sim_hbond_offset as fsim

# 小分子（1-ブタノール）に対して分子内疎水結合を許可した場合と許可しない場合を比較するテスト

def run_test_with_intramolecular(smiles, allow_intra=False):
    G_init, w_os = fsim.create_integrated_system([smiles], water_per_atom=3, is_loop=True)
    
    # run_fiedler_opt_sim のロジックを少し拡張してテスト
    G = G_init.copy()
    pos = nx.spring_layout(G, dim=3, k=0.1, seed=42)
    s_nodes = [n for n, d in G.nodes(data=True) if d.get('type') == 'solute_atom']
    
    f_history = []
    
    for t in range(50):
        # 分子内疎水結合の試行
        if t >= 20 and allow_intra:
            for i, u in enumerate(s_nodes):
                for v in s_nodes[i+1:]:
                    # 同一分子内かつ両方疎水アトム
                    if G.nodes[u].get('category') == 'hydrophobic' and G.nodes[v].get('category') == 'hydrophobic':
                        if not G.has_edge(u, v):
                            dist = np.linalg.norm(pos[u] - pos[v])
                            if dist < 0.20:
                                G.add_edge(u, v, edge_type='agg', color='orange')
                                
        pos = nx.spring_layout(G, dim=3, pos=pos, iterations=4, seed=42)
        evals, f_val, _ = fsim.analyze_laplacian_with_vector(G)
        f_history.append(f_val)
        
    return f_history[0], f_history[-1], f_history[-1] - f_history[0], G.number_of_edges()

import networkx as nx

init_no, final_no, delta_no, edges_no = run_test_with_intramolecular("CCCCO", allow_intra=False)
init_yes, final_yes, delta_yes, edges_yes = run_test_with_intramolecular("CCCCO", allow_intra=True)

print("=========================================================================")
print(" 🧪 Intra-molecular Hydrophobic Bonding Test (1-Butanol)")
print("=========================================================================")
print(f"・分子内疎水結合なし : f_init={init_no:.6f}, f_final={final_no:.6f}, Δf={delta_no:.6f}, エッジ数={edges_no}")
print(f"・分子内疎水結合あり : f_init={init_yes:.6f}, f_final={final_yes:.6f}, Δf={delta_yes:.6f}, エッジ数={edges_yes}")
print("=========================================================================")
