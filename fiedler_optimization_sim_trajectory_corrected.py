import networkx as nx
import plotly.graph_objects as go
from rdkit import Chem
import numpy as np
import random
import matplotlib.pyplot as plt

# --- 1. システム構築 (SMILES + 原子レベル水 + 数珠つなぎ) ---
def create_integrated_system(smiles_list, water_per_atom=2, is_loop=True):
    G = nx.Graph()
    all_water_os = []

    atom_props = {
        'H': {'c': 0.05, 'color': 'white', 'size': 4, 'cat': 'hydrophobic'},
        'C': {'c': 0.1, 'color': 'gray', 'size': 10, 'cat': 'hydrophobic'},
        'O': {'c': -0.4, 'color': 'red', 'size': 10, 'cat': 'polar'},
        'N': {'c': -0.3, 'color': 'blue', 'size': 10, 'cat': 'polar'},
        'P': {'c': 0.5, 'color': 'orange', 'size': 13, 'cat': 'polar'}, 
        'S': {'c': 0.4, 'color': 'yellow', 'size': 13, 'cat': 'polar'}, 
        'F': {'c': -0.2, 'color': 'lightgreen', 'size': 8, 'cat': 'polar'},
        'Cl': {'c': -0.5, 'color': 'green', 'size': 11, 'cat': 'polar'},
        'Br': {'c': -0.4, 'color': 'darkgreen', 'size': 13, 'cat': 'polar'},
        'Na': {'c': 0.9, 'color': 'purple', 'size': 12, 'cat': 'ion'},  
        'K':  {'c': 0.8, 'color': 'violet', 'size': 14, 'cat': 'ion'}
    }

    for mol_idx, smiles in enumerate(smiles_list):
        mol_name = f"MOL_{mol_idx}"
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        for atom in mol.GetAtoms():
            idx, sym = atom.GetIdx(), atom.GetSymbol()
            node_id = f"S_{mol_idx}_{idx}"
            
            prop = atom_props.get(sym, {'c': 0.0, 'color': 'pink', 'size': 10, 'cat': 'polar'})
            charge = prop['c'] + atom.GetFormalCharge()
            category = 'polar' if (prop['cat'] == 'polar' or atom.GetFormalCharge() != 0) else prop['cat']

            # 水素結合の最大許容数（次数制限）を設定
            if sym == 'O':
                max_h_bonds = 2
                num_waters = 2
            elif sym == 'N':
                max_h_bonds = 3
                num_waters = 3
            elif category in ['polar', 'ion']:
                max_h_bonds = 2
                num_waters = 2
            else:
                max_h_bonds = 1
                num_waters = 1

            G.add_node(node_id, symbol=sym, type='solute_atom', molecule_id=mol_name,
                       category=category, charge=charge, color=prop['color'], size=prop['size'],
                       max_h_bonds=max_h_bonds)
            
            # 初期ダミー水分子を次数制限（num_waters）に合わせて割り当てる
            for _ in range(num_waters):
                w_idx = len(all_water_os)
                o, h1, h2 = f"W{w_idx}_O", f"W{w_idx}_H1", f"W{w_idx}_H2"
                G.add_node(o, symbol='O', type='water_atom', color='blue', size=7, max_h_bonds=4)
                G.add_node(h1, symbol='H', type='water_atom', color='white', size=3)
                G.add_node(h2, symbol='H', type='water_atom', color='white', size=3)
                G.add_edge(o, h1, edge_type='covalent_w')
                G.add_edge(o, h2, edge_type='covalent_w')
                G.add_edge(node_id, o, edge_type='dummy_solvation')
                all_water_os.append(o)
                
        for bond in mol.GetBonds():
            G.add_edge(f"S_{mol_idx}_{bond.GetBeginAtomIdx()}", f"S_{mol_idx}_{bond.GetEndAtomIdx()}", edge_type='covalent')

    if all_water_os:
        for i in range(len(all_water_os) - 1):
            G.add_edge(all_water_os[i], all_water_os[i+1], edge_type='water_net', color='cyan')
        if is_loop and len(all_water_os) > 2:
            G.add_edge(all_water_os[-1], all_water_os[0], edge_type='water_net', color='cyan') 
            
    return G, all_water_os

# --- 2. スペクトル・Fiedlerベクトル解析 ---
def analyze_laplacian_with_vector(G):
    if not nx.is_connected(G): 
        return np.zeros(len(G)), 0, np.zeros(len(G))
    L = nx.laplacian_matrix(G).toarray()
    evals, evecs = np.linalg.eigh(L)
    if len(evals) > 1:
        fiedler = evals[1]
        fiedler_vec = evecs[:, 1]
    else:
        fiedler = 0
        fiedler_vec = np.zeros(len(G))
    return evals, fiedler, fiedler_vec

# --- 3. 軌跡に基づくFiedler値の数理補正式 (Trajectory-based Correction) ---
def calculate_corrected_fiedler_score(init, peak, final):
    r_peak_init = peak / init if init > 0 else 1.0
    r_final_init = final / init if init > 0 else 1.0
    r_peak_final = peak / final if final > 0 else 1.0
    r_final_peak = final / peak if peak > 0 else 1.0
    
    # BFGS局所最適化によって得られた重み係数
    w = [0.0083, -0.2211, -0.1934, 0.0740, -0.5618, -0.1102, -0.5859]
    score = (w[0] * init
             + w[1] * peak
             + w[2] * final
             + w[3] * r_peak_init
             + w[4] * r_final_init
             + w[5] * r_peak_final
             + w[6] * r_final_peak)
    return score

# --- 4. 最適化志向AIシミュレーター ---
def run_fiedler_opt_sim(G_init, w_os, steps=100, allow_long_distance_water=False):
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
        evals, f_val, f_vec = analyze_laplacian_with_vector(G)
        
        # --- 新規: Fiedlerベクトルを利用したトポロジー最適化フェーズ ---
        if t > 10 and nx.is_connected(G):
            # ■ 1. エッジの切断
            to_remove = []
            for u, v, d in G.edges(data=True):
                etype = d.get('edge_type')
                
                if etype in ['opt_water', 'dummy_solvation']:
                    is_polar_dummy = False
                    if etype == 'dummy_solvation':
                        s_node = u if 'S_' in str(u) else v
                        if G.nodes[s_node].get('category') in ['polar', 'ion']:
                            is_polar_dummy = True
                    
                    u_idx, v_idx = node_to_idx.get(u), node_to_idx.get(v)
                    if u_idx is not None and v_idx is not None:
                        diff = abs(f_vec[u_idx] - f_vec[v_idx])
                        dist = np.linalg.norm(pos[u] - pos[v])
                        is_ww = (u in w_os and v in w_os)
                        
                        over_dist = dist > 0.8
                        if allow_long_distance_water and is_ww:
                            over_dist = False
                            
                        should_cut = over_dist
                        if not is_polar_dummy and (diff < 0.05 and random.random() < 0.2):
                            should_cut = True
                            
                        if should_cut and G.degree(u) > 1 and G.degree(v) > 1:
                            to_remove.append((u, v))
                            
            G.remove_edges_from(to_remove)
            
            nodelist = list(G.nodes())
            node_to_idx = {n: i for i, n in enumerate(nodelist)}
            evals, f_val, f_vec = analyze_laplacian_with_vector(G)

            # ■ 2. エッジの追加
            if nx.is_connected(G):
                degrees = {n: sum(1 for neighbor in G.neighbors(n) if G.edges[n, neighbor].get('edge_type') in ['opt_water', 'dummy_solvation', 'water_net']) for n in candidate_nodes}
                
                potentials = []
                for i in range(len(candidate_nodes)):
                    for j in range(i+1, len(candidate_nodes)):
                        u, v = candidate_nodes[i], candidate_nodes[j]
                        if not G.has_edge(u, v):
                            max_u = G.nodes[u].get('max_h_bonds', 4)
                            max_v = G.nodes[v].get('max_h_bonds', 4)
                            
                            if degrees.get(u, 0) < max_u and degrees.get(v, 0) < max_v:
                                dist = np.linalg.norm(pos[u] - pos[v])
                                is_ww = (u in w_os and v in w_os)
                                
                                valid_dist = dist < 0.6
                                if allow_long_distance_water and is_ww:
                                    valid_dist = True
                                
                                if valid_dist:
                                    u_idx, v_idx = node_to_idx[u], node_to_idx[v]
                                    diff = abs(f_vec[u_idx] - f_vec[v_idx])
                                    potentials.append((diff, u, v))
                                    
                potentials.sort(reverse=True, key=lambda x: x[0])
                added_count = 0
                for diff, u, v in potentials:
                    G.add_edge(u, v, edge_type='opt_water', color='lime')
                    degrees[u] = degrees.get(u, 0) + 1
                    degrees[v] = degrees.get(v, 0) + 1
                    added_count += 1
                    if added_count >= 3:
                        break

        # --- フェーズ2: 疎水性排除 ---
        if t >= 20:
            to_remove_hydro = []
            for u, v, d in G.edges(data=True):
                if d.get('edge_type') == 'dummy_solvation':
                    s_node = u if 'S_' in str(u) else v
                    if G.nodes[s_node].get('category') == 'hydrophobic':
                        if random.random() < 0.5:
                            to_remove_hydro.append((u, v))
            G.remove_edges_from(to_remove_hydro)

        # --- フェーズ3: 疎水性凝集 ---
        if t >= 30:
            for i, u in enumerate(s_nodes):
                for v in s_nodes[i+1:]:
                    if G.nodes[u].get('molecule_id') != G.nodes[v].get('molecule_id') and \
                       G.nodes[u].get('category') == 'hydrophobic' and np.linalg.norm(pos[u] - pos[v]) < 0.2:
                        G.add_edge(u, v, edge_type='agg', color='orange')

        pos = nx.spring_layout(G, dim=3, pos=pos, iterations=4, seed=42)
        
        evals, f_val, _ = analyze_laplacian_with_vector(G)
        f_history.append(f_val)
        frames_data.append((G.copy(), pos.copy(), evals))

    # 最終的な軌跡ベースの補正スコアを算出
    corrected_score = calculate_corrected_fiedler_score(f_history[0], np.max(f_history), f_history[-1])

    return f_history, frames_data, corrected_score
