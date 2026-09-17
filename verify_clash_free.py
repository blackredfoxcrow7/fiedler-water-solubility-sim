import sys
import os
import networkx as nx
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import fiedler_sim_hbond_offset_clash_free as fsim

def verify_clash_free(smiles, mol_name):
    G_init, w_os = fsim.create_integrated_system([smiles], water_per_atom=3, is_loop=True)
    _, frames_data, delta_f_final, _ = fsim.run_fiedler_opt_sim(G_init, w_os, steps=100)
    
    G_final, pos_final, _ = frames_data[-1]
    nodelist = list(G_final.nodes())
    
    clashes = []
    min_dist_found = float('inf')
    
    for i in range(len(nodelist)):
        for j in range(i+1, len(nodelist)):
            u, v = nodelist[i], nodelist[j]
            sym_u = G_final.nodes[u].get('symbol', '')
            sym_v = G_final.nodes[v].get('symbol', '')
            
            if G_final.has_edge(u, v) and G_final.edges[u, v].get('edge_type') in ['covalent', 'covalent_w']:
                continue
                
            dist = np.linalg.norm(pos_final[u] - pos_final[v])
            if dist < min_dist_found:
                min_dist_found = dist
                
            limit = 0.045 if (sym_u == 'H' or sym_v == 'H') else 0.065
            if dist < limit:
                clashes.append((u, v, sym_u, sv, dist, limit))
                
    print(f"\n==================================================")
    print(f" 🛡️ Clash Free Verification: {mol_name} ({smiles})")
    print(f"==================================================")
    print(f"  ・最小非結合間距離  : {min_dist_found:.6f}")
    print(f"  ・検出された衝突数  : {len(clashes)} (0であれば完璧)")
    print(f"  ・実効変化量 Δf_final: {delta_f_final:.6f}")

verify_clash_free("CC(C)(C)O", "tert-Butanol")
verify_clash_free("CCOCC", "Diethyl Ether")
