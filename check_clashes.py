import sys
import os
import networkx as nx
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import fiedler_sim_hbond_offset as fsim

# fiedler_sim_hbond_offset.py の最終フレームにおける原子間めり込み（衝突）の判定
def check_atom_clashes(smiles, mol_name):
    G_init, w_os = fsim.create_integrated_system([smiles], water_per_atom=3, is_loop=True)
    _, frames_data, _, _ = fsim.run_fiedler_opt_sim(G_init, w_os, steps=100)
    
    G_final, pos_final, _ = frames_data[-1]
    nodelist = list(G_final.nodes())
    
    # 判定しきい値 (無次元スケール)
    # 重原子同士: 0.065 未満はめり込み
    # 水素含む: 0.045 未満はめり込み
    clashes = []
    min_dist_found = float('inf')
    
    for i in range(len(nodelist)):
        for j in range(i+1, len(nodelist)):
            u, v = nodelist[i], nodelist[j]
            sym_u = G_final.nodes[u].get('symbol', '')
            sym_v = G_final.nodes[v].get('symbol', '')
            
            # 共有結合しているペアは除外
            if G_final.has_edge(u, v) and G_final.edges[u, v].get('edge_type') in ['covalent', 'covalent_w']:
                continue
                
            dist = np.linalg.norm(pos_final[u] - pos_final[v])
            if dist < min_dist_found:
                min_dist_found = dist
                
            limit = 0.045 if (sym_u == 'H' or sym_v == 'H') else 0.065
            if dist < limit:
                clashes.append((u, v, sym_u, sym_v, dist, limit))
                
    print(f"\n--- Clash Test for {mol_name} ---")
    print(f"  ・最小非結合間距離  : {min_dist_found:.6f}")
    print(f"  ・検出された衝突数  : {len(clashes)}")
    if clashes:
        for u, v, su, sv, d, lim in clashes[:5]:
            print(f"    - 衝突: {u}({su}) - {v}({sv}) | 距離={d:.6f} < 限界値={lim:.3f}")
    return len(clashes)

check_atom_clashes("CC(C)(C)O", "tert-Butanol")
check_atom_clashes("CCOCC", "Diethyl Ether")
