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

            G.add_node(node_id, symbol=sym, type='solute_atom', molecule_id=mol_name,
                       category=category, charge=charge, color=prop['color'], size=prop['size'])
            
            for _ in range(water_per_atom):
                w_idx = len(all_water_os)
                o, h1, h2 = f"W{w_idx}_O", f"W{w_idx}_H1", f"W{w_idx}_H2"
                G.add_node(o, symbol='O', type='water_atom', color='blue', size=7)
                G.add_node(h1, symbol='H', type='water_atom', color='white', size=3)
                G.add_node(h2, symbol='H', type='water_atom', color='white', size=3)
                G.add_edge(o, h1, edge_type='covalent_w', weight=5.0)
                G.add_edge(o, h2, edge_type='covalent_w', weight=5.0)
                G.add_edge(node_id, o, edge_type='dummy_solvation')
                all_water_os.append(o)
                
        for bond in mol.GetBonds():
            G.add_edge(f"S_{mol_idx}_{bond.GetBeginAtomIdx()}", f"S_{mol_idx}_{bond.GetEndAtomIdx()}", edge_type='covalent', weight=5.0)

    if all_water_os:
        for i in range(len(all_water_os) - 1):
            G.add_edge(all_water_os[i], all_water_os[i+1], edge_type='water_net', color='cyan')
        if is_loop and len(all_water_os) > 2:
            G.add_edge(all_water_os[-1], all_water_os[0], edge_type='water_net', color='cyan') 
            
    return G, all_water_os

# --- 2. スペクトル・Fiedlerベクトル解析 ---
def analyze_laplacian_with_vector(G):
    """
    Fiedler値だけでなく、Fiedlerベクトル（最適化の指針）も同時に計算する拡張版
    """
    if not nx.is_connected(G): 
        return np.zeros(len(G)), 0, np.zeros(len(G))
    L = nx.laplacian_matrix(G).toarray()
    
    # np.linalg.eigh は実対称行列の固有値（昇順）と固有ベクトルを返す
    evals, evecs = np.linalg.eigh(L)
    if len(evals) > 1:
        fiedler = evals[1]
        fiedler_vec = evecs[:, 1] # 第2固有バクトル
    else:
        fiedler = 0
        fiedler_vec = np.zeros(len(G))
    return evals, fiedler, fiedler_vec

# --- 3. 最適化志向AIシミュレーター ---
def run_fiedler_opt_sim(G_init, w_os, steps=100, allow_long_distance_water=False):
    G = G_init.copy()
    pos = nx.spring_layout(G, dim=3, k=0.1, seed=42)
    s_nodes = [n for n, d in G.nodes(data=True) if d.get('type') == 'solute_atom']
    polar_s_nodes = [n for n in s_nodes if G.nodes[n].get('category') in ['polar', 'ion']]
    
    # 水の最適化探索の対象グループ
    candidate_nodes = w_os + polar_s_nodes
    
    f_history = []
    frames_data = []

    mode_text = "[距離無視]" if allow_long_distance_water else "[距離制約あり]"
    print(f"Fiedler駆動最適化シミュレーション中 {mode_text} ... (全 {steps} ステップ)")

    for t in range(steps):
        # 現状のFiedlerベクトルを計算
        nodelist = list(G.nodes())
        node_to_idx = {n: i for i, n in enumerate(nodelist)}
        evals, f_val, f_vec = analyze_laplacian_with_vector(G)
        
        # --- 新規: Fiedlerベクトルを利用したトポロジー最適化フェーズ ---
        # 水ネットワークが分断されていない前提で（Step 10以降本格化）
        if t > 10 and nx.is_connected(G):
            # ■ 1. エッジの切断（Fiedler値の上昇に寄与していない不要な結合の排除）
            to_remove = []
            for u, v, d in G.edges(data=True):
                etype = d.get('edge_type')
                
                if etype in ['opt_water', 'dummy_solvation']:
                    
                    # 絶対保護：極性基（アルコールのOなど）と繋がっているエッジは距離が離れても絶対に切断しない
                    is_polar_dummy = False
                    if etype == 'dummy_solvation':
                        s_node = u if 'S_' in str(u) else v
                        if G.nodes[s_node].get('category') in ['polar', 'ion']:
                            is_polar_dummy = True
                    
                    if not is_polar_dummy:
                        # どちらかのノードが孤立する（つながりがゼロになる）切断は許可しない
                        if G.degree(u) > 1 and G.degree(v) > 1:
                            u_idx, v_idx = node_to_idx.get(u), node_to_idx.get(v)
                            if u_idx is not None and v_idx is not None:
                                diff = abs(f_vec[u_idx] - f_vec[v_idx])
                                dist = np.linalg.norm(pos[u] - pos[v])
                                is_ww = (u in w_os and v in w_os)
                                
                                # 削除条件の判定
                                over_dist = dist > 0.8
                                if allow_long_distance_water and is_ww:
                                    over_dist = False # フラグONかつ水同士なら距離理由では切られない
                                    
                                # ①距離が遠い、または ②ベクトル差が小さい（無価値）かつ揺らぎで切断
                                if over_dist or (diff < 0.05 and random.random() < 0.2):
                                    to_remove.append((u, v))
            G.remove_edges_from(to_remove)
            
            # トポロジー変化後、ベクトルを再計算
            nodelist = list(G.nodes())
            node_to_idx = {n: i for i, n in enumerate(nodelist)}
            evals, f_val, f_vec = analyze_laplacian_with_vector(G)

            # ■ 2. エッジの追加（Fiedler値が最も上昇する特効薬的な結合の探索）
            if nx.is_connected(G):
                # 水や極性基の結合数（最大配位数=4まで）を計測
                degrees = {n: sum(1 for neighbor in G.neighbors(n) if G.edges[n, neighbor].get('edge_type') in ['opt_water', 'dummy_solvation', 'water_net']) for n in candidate_nodes}
                
                potentials = []
                for i in range(len(candidate_nodes)):
                    for j in range(i+1, len(candidate_nodes)):
                        u, v = candidate_nodes[i], candidate_nodes[j]
                        if not G.has_edge(u, v):
                            if degrees.get(u, 0) < 4 and degrees.get(v, 0) < 4:
                                dist = np.linalg.norm(pos[u] - pos[v])
                                is_ww = (u in w_os and v in w_os)
                                
                                valid_dist = dist < 0.6
                                if allow_long_distance_water and is_ww:
                                    valid_dist = True # 距離無視フラグON時は空間を飛び越える
                                
                                if valid_dist:
                                    u_idx, v_idx = node_to_idx[u], node_to_idx[v]
                                    diff = abs(f_vec[u_idx] - f_vec[v_idx])
                                    potentials.append((diff, u, v))
                                    
                # 最もFiedler値上昇に貢献する上位3ペアを結ぶ
                potentials.sort(reverse=True, key=lambda x: x[0])
                added_count = 0
                for diff, u, v in potentials:
                    G.add_edge(u, v, edge_type='opt_water', color='lime', weight=1.0)
                    degrees[u] = degrees.get(u, 0) + 1
                    degrees[v] = degrees.get(v, 0) + 1
                    added_count += 1
                    if added_count >= 3:
                        break

        # --- フェーズ2: 疎水性排除（水と疎水基の初期結合の完全切断） ---
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
        # 潰れを防ぐ上限は残しつつ、各部位がしっかり密着するように「上限本数」と「引力」を調整
        if t >= 30:
            agg_degrees = {n: 0 for n in s_nodes}
            for u, v, d in G.edges(data=True):
                if d.get('edge_type') == 'agg':
                    if u in agg_degrees: agg_degrees[u] += 1
                    if v in agg_degrees: agg_degrees[v] += 1

            for i, u in enumerate(s_nodes):
                for v in s_nodes[i+1:]:
                    if G.nodes[u].get('molecule_id') != G.nodes[v].get('molecule_id') and \
                       G.nodes[u].get('category') == 'hydrophobic' and np.linalg.norm(pos[u] - pos[v]) < 0.3:
                        # 凝集上限を4本へ緩和し、もう少し密に絡み合えるようにする
                        if not G.has_edge(u, v) and agg_degrees[u] < 4 and agg_degrees[v] < 4:
                            # バネの重みを引き上げ（1.5 -> 2.5）、適度な距離までしっかりと引き寄せる
                            G.add_edge(u, v, edge_type='agg', color='orange', weight=2.5)
                            agg_degrees[u] += 1
                            agg_degrees[v] += 1



        # ばねモデルによる形状再配置
        pos = nx.spring_layout(G, dim=3, pos=pos, iterations=4, seed=42)
        
        # 記録用
        evals, f_val, _ = analyze_laplacian_with_vector(G)
        f_history.append(f_val)
        frames_data.append((G.copy(), pos.copy(), evals))

    return f_history, frames_data

# --- 4. 可視化モジュール ---
def visualize_results(f_history, frames_data, allow_long_distance_water=False):
    frames = []
    for t, (G_t, pos_t, _) in enumerate(frames_data):
        edge_x, edge_y, edge_z = [], [], []
        agg_x, agg_y, agg_z = [], [], []
        opt_x, opt_y, opt_z = [], [], [] # 最適化された新しい水素結合（黄緑色）
        
        for u, v, d in G_t.edges(data=True):
            p1, p2 = pos_t[u], pos_t[v]
            etype = d.get('edge_type')
            if etype == 'agg':
                agg_x.extend([p1[0], p2[0], None])
                agg_y.extend([p1[1], p2[1], None])
                agg_z.extend([p1[2], p2[2], None])
            elif etype == 'opt_water':
                opt_x.extend([p1[0], p2[0], None])
                opt_y.extend([p1[1], p2[1], None])
                opt_z.extend([p1[2], p2[2], None])
            else:
                edge_x.extend([p1[0], p2[0], None])
                edge_y.extend([p1[1], p2[1], None])
                edge_z.extend([p1[2], p2[2], None])
        
        node_x = [pos_t[n][0] for n in G_t.nodes()]
        node_y = [pos_t[n][1] for n in G_t.nodes()]
        node_z = [pos_t[n][2] for n in G_t.nodes()]
        node_colors = [G_t.nodes[n].get('color', 'green') for n in G_t.nodes()]
        node_sizes = [G_t.nodes[n].get('size', 5) for n in G_t.nodes()]
        
        plot_data = [
            go.Scatter3d(x=edge_x, y=edge_y, z=edge_z, mode='lines', line=dict(color='cyan', width=1), opacity=0.4),
            go.Scatter3d(x=agg_x, y=agg_y, z=agg_z, mode='lines', line=dict(color='orange', width=2), opacity=0.9),
            go.Scatter3d(x=opt_x, y=opt_y, z=opt_z, mode='lines', line=dict(color='lime', width=2), opacity=0.9), # AI生成エッジ
            go.Scatter3d(x=node_x, y=node_y, z=node_z, mode='markers', marker=dict(size=node_sizes, color=node_colors))
        ]

        frames.append(go.Frame(data=plot_data, name=f'step{t}'))

    mode_text = " (Long-Distance OFF)" if not allow_long_distance_water else " (Long-Distance ON)"

    sliders = [dict(
        active=0, yanchor="top", xanchor="left",
        currentvalue=dict(font=dict(size=14), prefix="Simulation Step: ", visible=True, xanchor="right"),
        transition=dict(duration=100, easing="cubic-in-out"),
        pad=dict(b=10, t=50), len=0.9, x=0.1, y=0,
        steps=[dict(
            args=[[f.name], dict(frame=dict(duration=100, redraw=True), mode="immediate", transition=dict(duration=0))],
            label=str(k), method="animate"
        ) for k, f in enumerate(frames)]
    )]

    updatemenus = [dict(
        type="buttons", showactive=False, y=0, x=0, xanchor="right", yanchor="top", pad=dict(t=50, r=10),
        buttons=[
            dict(label="Play", method="animate", args=[None, dict(frame=dict(duration=100, redraw=True), fromcurrent=True, transition=dict(duration=0))]),
            dict(label="Pause", method="animate", args=[[None], dict(frame=dict(duration=0, redraw=False), mode="immediate", transition=dict(duration=0))])
        ]
    )]

    fig = go.Figure(
        data=frames[0].data if frames else [],
        layout=go.Layout(
            title=f"Topology Optimization via Fiedler Vector{mode_text}",
            scene=dict(xaxis_visible=False, yaxis_visible=False, zaxis_visible=False),
            updatemenus=updatemenus, sliders=sliders
        ),
        frames=frames
    )
    fig.show()

    best_step = np.argmax(f_history)
    _, _, best_evals = frames_data[best_step]

    plt.figure(figsize=(10, 4))
    plt.subplot(1, 2, 1)
    plt.plot(f_history, color='blue')
    plt.axvline(x=10, color='gray', linestyle=':', label='Optimization Start')
    plt.axvline(x=best_step, color='red', linestyle='--', label=f'Peak Stability')
    plt.title("Evolution by Fiedler Maximization")
    plt.xlabel("Step")
    plt.legend()
    plt.grid(True)
    
    plt.subplot(1, 2, 2)
    plt.hist(best_evals, bins=80, color='purple', alpha=0.7)
    plt.title(f"Eigenvalue Spikes (Best State at Step {best_step})")
    plt.xlabel("Eigenvalue")
    plt.grid(True)
    
    plt.tight_layout()
    plt.show()

    # --- 最終結果のコンソール出力 ---
    final_fiedler = f_history[-1]
    peak_fiedler = np.max(f_history)
    print("\n" + "="*40)
    print(" 🏁 Simulation Complete")
    print("="*40)
    print(f"   Peak Fiedler Value  : {peak_fiedler:.6f} (at Step {best_step})")
    print(f"   Final Fiedler Value : {final_fiedler:.6f}")
    print("="*40 + "\n")

# --- 実行用コード（フラグ例） ---
if __name__ == "__main__":
    smiles_list = ["CC(=O)N[C@@H](C)C(=O)N[C@@H](C)C(=O)N[C@@H](C)C(=O)N[C@@H](C)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](C)C(=O)N[C@@H](C)C(=O)N[C@@H](C)C(=O)N[C@@H](C)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](C)C(=O)N[C@@H](C)C(=O)N[C@@H](C)C(=O)N[C@@H](C)C(=O)N[C@@H](CCCCN)C(N)=O"]#* 6#, "P(=O)(O)(O)O", "C1=CC=CC=C1"]"CC(O)C""CC(C)(C)O"
    print(smiles_list)
    ALLOW_LONG_DISTANCE = False # Trueで水同士が空間を越えて繋がるワープエッジを許可["CCCCCCO"] * 6

    G_init, w_os = create_integrated_system(smiles_list, water_per_atom=4, is_loop=True)
    f_history, frames_data = run_fiedler_opt_sim(G_init, w_os, steps=100, allow_long_distance_water=ALLOW_LONG_DISTANCE)
    visualize_results(f_history, frames_data, allow_long_distance_water=ALLOW_LONG_DISTANCE)
