import networkx as nx
import plotly.graph_objects as go
from rdkit import Chem
import numpy as np
import random

import fiedler_sim_hbond_offset_clash_free as fsim_base

# --- 積極的に水網をフラクタル化する拡張シミュレータ ---
def run_fiedler_fractal_opt_sim(G_init, w_os, steps=100, alpha_powerlow=2.0):
    G = G_init.copy()
    pos = nx.spring_layout(G, dim=3, k=0.1, seed=42)
    s_nodes = [n for n, d in G.nodes(data=True) if d.get('type') == 'solute_atom']
    polar_s_nodes = [n for n in s_nodes if G.nodes[n].get('category') in ['polar', 'ion']]
    
    candidate_nodes = w_os + polar_s_nodes
    f_history = []
    frames_data = []

    for t in range(steps):
        nodelist = list(G.nodes())
        node_to_idx = {n: i for i, n in enumerate(nodelist)}
        evals, f_val, f_vec = fsim_base.analyze_laplacian_with_vector(G)
        
        if t > 10 and nx.is_connected(G):
            # ■ 1. エッジ切断
            to_remove = []
            for u, v, d in G.edges(data=True):
                etype = d.get('edge_type')
                if etype in ['opt_water', 'dummy_solvation']:
                    u_idx, v_idx = node_to_idx.get(u), node_to_idx.get(v)
                    if u_idx is not None and v_idx is not None:
                        dist = np.linalg.norm(pos[u] - pos[v])
                        if dist > 0.8 and G.degree(u) > 1 and G.degree(v) > 1:
                            to_remove.append((u, v))
            G.remove_edges_from(to_remove)
            
            nodelist = list(G.nodes())
            node_to_idx = {n: i for i, n in enumerate(nodelist)}
            evals, f_val, f_vec = fsim_base.analyze_laplacian_with_vector(G)

            # ■ 2. 【フラクタル活性化】パワーロー距離 ＆ 優先的付着ルールによる水網形成
            if nx.is_connected(G):
                degrees = {n: G.degree(n) for n in candidate_nodes}
                potentials = []
                
                for i in range(len(candidate_nodes)):
                    for j in range(i+1, len(candidate_nodes)):
                        u, v = candidate_nodes[i], candidate_nodes[j]
                        if not G.has_edge(u, v):
                            max_u = G.nodes[u].get('max_h_bonds', 4)
                            max_v = G.nodes[v].get('max_h_bonds', 4)
                            
                            if degrees.get(u, 0) < max_u and degrees.get(v, 0) < max_v:
                                dist = np.linalg.norm(pos[u] - pos[v]) + 1e-5
                                
                                # フラクタルアトラクション確率 (P ~ r^(-alpha) * (deg_u * deg_v)^gamma)
                                pref_factor = (degrees.get(u, 1) * degrees.get(v, 1)) ** 0.5
                                fractal_score = (dist ** (-alpha_powerlow)) * pref_factor
                                
                                u_idx, v_idx = node_to_idx[u], node_to_idx[v]
                                diff = abs(f_vec[u_idx] - f_vec[v_idx])
                                
                                # トポロジー差とフラクタル優先度の合体
                                total_score = diff * 0.5 + fractal_score * 0.5
                                potentials.append((total_score, u, v))
                                    
                potentials.sort(reverse=True, key=lambda x: x[0])
                added_count = 0
                for score, u, v in potentials:
                    G.add_edge(u, v, edge_type='opt_water', color='cyan')
                    degrees[u] = degrees.get(u, 0) + 1
                    degrees[v] = degrees.get(v, 0) + 1
                    added_count += 1
                    if added_count >= 3:
                        break

        # フェーズ2: 疎水性排除
        if t >= 20:
            to_remove_hydro = []
            for u, v, d in G.edges(data=True):
                if d.get('edge_type') == 'dummy_solvation':
                    s_node = u if 'S_' in str(u) else v
                    if G.nodes[s_node].get('category') == 'hydrophobic' and random.random() < 0.5:
                        to_remove_hydro.append((u, v))
            G.remove_edges_from(to_remove_hydro)

        pos = nx.spring_layout(G, dim=3, pos=pos, iterations=4, seed=42)
        pos = fsim_base.apply_clash_avoidance(G, pos)
        
        evals, f_val, _ = fsim_base.analyze_laplacian_with_vector(G)
        f_history.append(f_val)
        frames_data.append((G.copy(), pos.copy(), evals))

    f_offset = f_history[0]
    delta_f_final = f_history[-1] - f_offset
    return f_history, frames_data, delta_f_final
