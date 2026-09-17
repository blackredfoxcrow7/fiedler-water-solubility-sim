"""
visualize_dimer_fiedler_3d.py
=============================================================================
ダイマー化状態における Fiedler 固有ベクトルの符号・数値および分子間結合の3D視覚化

【機能】
1. 2分子ダイマー（Solute A + Solute B）のKabsch最適化後3D位置関係を描画
2. 各原子ノードを Fiedler 固有ベクトルの符号（正＝赤 / 負＝青）と数値でカラーマップ表示
3. 分子間相互作用（塩橋・水素結合・π-πスタッキング・疎水結合）を種類別に色分けした破線描画
4. 位相ベクトル d_phase の3D矢印表示
5. Plotly によるインタラクティブ 3D HTML ファイル（3モード切替付き）の保存
=============================================================================
"""

import sys
import os
import networkx as nx
import numpy as np
import plotly.graph_objects as go
from rdkit import Chem

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import fiedler_dimer_interaction_optimizer as fdio

def generate_dimer_fiedler_3d_html(name, smiles, out_filename):
    mol, coords = fdio.get_3d_conformer(smiles)
    N = mol.GetNumAtoms()
    
    # ---------------------------------------------------------
    # 1. 2分子ダイマーのグラフ構築と Fiedler 固有ベクトルの算出
    # ---------------------------------------------------------
    G_dimer = nx.Graph()
    for i in range(N):
        sym = mol.GetAtomWithIdx(i).GetSymbol()
        G_dimer.add_node(f"A_{i}", symbol=sym, mol='A', orig_idx=i)
        G_dimer.add_node(f"B_{i}", symbol=sym, mol='B', orig_idx=i)
        
    for bond in mol.GetBonds():
        G_dimer.add_edge(f"A_{bond.GetBeginAtomIdx()}", f"A_{bond.GetEndAtomIdx()}", weight=1.0, edge_type='covalent')
        G_dimer.add_edge(f"B_{bond.GetBeginAtomIdx()}", f"B_{bond.GetEndAtomIdx()}", weight=1.0, edge_type='covalent')

    # 候補分子間エッジの生成と貪欲Fiedler最大化
    candidate_inter_edges = []
    for a_atom in mol.GetAtoms():
        a_idx = a_atom.GetIdx()
        a_sym = a_atom.GetSymbol()
        a_aro = a_atom.GetIsAromatic()
        a_h = sum(1 for nbr in a_atom.GetNeighbors() if nbr.GetSymbol() == 'H')
        for b_atom in mol.GetAtoms():
            b_idx = b_atom.GetIdx()
            b_sym = b_atom.GetSymbol()
            b_aro = b_atom.GetIsAromatic()
            b_h = sum(1 for nbr in b_atom.GetNeighbors() if nbr.GetSymbol() == 'H')
            
            node_A = f"A_{a_idx}"
            node_B = f"B_{b_idx}"
            
            if (a_sym == 'O' and b_sym == 'N') or (a_sym == 'N' and b_sym == 'O'):
                candidate_inter_edges.append((node_A, node_B, 2.0, 'ionic_bridge', '塩橋/イオン結合'))
            elif (a_sym in ['O', 'N'] and b_sym in ['O', 'N']) and (a_h > 0 or b_h > 0):
                candidate_inter_edges.append((node_A, node_B, 1.5, 'hbond', '水素結合'))
            elif a_aro and b_aro and a_sym == 'C' and b_sym == 'C':
                candidate_inter_edges.append((node_A, node_B, 1.2, 'pi_stacking', '芳香環π-π結合'))
            elif a_sym == 'C' and b_sym == 'C' and (not a_aro) and (not b_aro):
                candidate_inter_edges.append((node_A, node_B, 0.7, 'hydrophobic', '疎水結合'))

    candidate_inter_edges.sort(key=lambda x: x[2], reverse=True)
    selected_edges = []
    degree_A = {f"A_{i}": 0 for i in range(N)}
    degree_B = {f"B_{i}": 0 for i in range(N)}
    
    current_best_f = 0.0
    for u, v, w, etype, label in candidate_inter_edges:
        if degree_A[u] < 2 and degree_B[v] < 2:
            G_dimer.add_edge(u, v, weight=w, edge_type=etype, label=label)
            f_test = fdio.compute_fiedler_normalized(G_dimer)
            if f_test > current_best_f or not nx.is_connected(G_dimer):
                current_best_f = f_test
                degree_A[u] += 1
                degree_B[v] += 1
                selected_edges.append((u, v, w, etype, label))
            else:
                G_dimer.remove_edge(u, v)
            if len(selected_edges) >= min(12, N):
                break

    # 正規化ラプラシアンの固有ベクトル分析
    L_dimer = nx.normalized_laplacian_matrix(G_dimer).toarray()
    evals, evecs = np.linalg.eigh(L_dimer)
    fiedler_val = evals[1]
    fiedler_vec = evecs[:, 1]
    
    nodes_order = list(G_dimer.nodes())
    node_to_fiedler = {node: fiedler_vec[idx] for idx, node in enumerate(nodes_order)}
    
    # ---------------------------------------------------------
    # 2. Kabsch アルゴリズムによる Molecule B の3D配置最適化
    # ---------------------------------------------------------
    coords_A = coords.copy()
    coords_B = coords.copy()
    
    if len(selected_edges) > 0:
        P = np.array([coords_A[int(u.split('_')[1])] for u, v, w, etype, label in selected_edges])
        Q = np.array([coords_B[int(v.split('_')[1])] for u, v, w, etype, label in selected_edges])
        W = np.array([w for u, v, w, etype, label in selected_edges])
        R, t = fdio.kabsch_rigid_transform(P, Q, weights=W)
        coords_B_opt = np.dot(coords_B, R.T) + t
    else:
        coords_B_opt = coords_B + np.array([3.0, 0.0, 0.0])

    pos_3d = {}
    for i in range(N):
        pos_3d[f"A_{i}"] = coords_A[i]
        pos_3d[f"B_{i}"] = coords_B_opt[i]

    # ---------------------------------------------------------
    # 3. Plotly 3D インタラクティブ視覚化オブジェクトの作成
    # ---------------------------------------------------------
    fig = go.Figure()
    
    # max abs value for color scaling
    max_abs_f = max(abs(v) for v in fiedler_vec) if len(fiedler_vec) > 0 else 1.0

    # Trace 0: Solute A ノード (Fiedler 符号・数値カラー)
    x_A, y_A, z_A, c_A, text_A = [], [], [], [], []
    for i in range(N):
        node = f"A_{i}"
        x, y, z = pos_3d[node]
        val = node_to_fiedler[node]
        sign_str = "(+ 正)" if val > 0 else "(- 負)"
        x_A.append(x); y_A.append(y); z_A.append(z)
        c_A.append(val)
        text_A.append(f"A_{i} ({mol.GetAtomWithIdx(i).GetSymbol()})<br>Fiedler: {val:+.4f} {sign_str}")

    # Trace 1: Solute B ノード
    x_B, y_B, z_B, c_B, text_B = [], [], [], [], []
    for i in range(N):
        node = f"B_{i}"
        x, y, z = pos_3d[node]
        val = node_to_fiedler[node]
        sign_str = "(+ 正)" if val > 0 else "(- 負)"
        x_B.append(x); y_B.append(y); z_B.append(z)
        c_B.append(val)
        text_B.append(f"B_{i} ({mol.GetAtomWithIdx(i).GetSymbol()})<br>Fiedler: {val:+.4f} {sign_str}")

    # Add Molecule A Nodes
    fig.add_trace(go.Scatter3d(
        x=x_A, y=y_A, z=z_A, mode='markers+text',
        marker=dict(
            size=11, color=c_A, colorscale='RdBu_r', cmin=-max_abs_f, cmax=max_abs_f,
            colorbar=dict(title="Fiedler Vector Value (符号)", x=1.05),
            line=dict(color='black', width=1.5)
        ),
        text=[f"A:{mol.GetAtomWithIdx(i).GetSymbol()}" for i in range(N)],
        hoverinfo='text', hovertext=text_A, name="Molecule A (Solute A)", visible=True
    ))

    # Add Molecule B Nodes
    fig.add_trace(go.Scatter3d(
        x=x_B, y=y_B, z=z_B, mode='markers+text',
        marker=dict(
            size=11, color=c_B, colorscale='RdBu_r', cmin=-max_abs_f, cmax=max_abs_f,
            line=dict(color='yellow', width=1.5)
        ),
        text=[f"B:{mol.GetAtomWithIdx(i).GetSymbol()}" for i in range(N)],
        hoverinfo='text', hovertext=text_B, name="Molecule B (Solute B)", visible=True
    ))

    # Trace 2: 分子内共有結合 (Molecule A & B Covalent Sticks)
    cov_x, cov_y, cov_z = [], [], []
    for u, v, d in G_dimer.edges(data=True):
        if d.get('edge_type') == 'covalent':
            x0, y0, z0 = pos_3d[u]
            x1, y1, z1 = pos_3d[v]
            cov_x.extend([x0, x1, None]); cov_y.extend([y0, y1, None]); cov_z.extend([z0, z1, None])

    fig.add_trace(go.Scatter3d(
        x=cov_x, y=cov_y, z=cov_z, mode='lines',
        line=dict(color='lightgray', width=6), name='Covalent Bonds', visible=True
    ))

    # Trace 3: 分子間相互作用エッジ (Intermolecular Edges by Type)
    edge_colors = {
        'ionic_bridge': 'magenta',
        'hbond': 'lime',
        'pi_stacking': 'gold',
        'hydrophobic': 'orange'
    }
    
    for etype_key, color_code in edge_colors.items():
        inter_x, inter_y, inter_z = [], [], []
        labels_list = []
        for u, v, d in G_dimer.edges(data=True):
            if d.get('edge_type') == etype_key:
                x0, y0, z0 = pos_3d[u]
                x1, y1, z1 = pos_3d[v]
                dist = np.linalg.norm(pos_3d[u] - pos_3d[v])
                inter_x.extend([x0, x1, None])
                inter_y.extend([y0, y1, None])
                inter_z.extend([z0, z1, None])
                labels_list.append(d.get('label', etype_key))

        if inter_x:
            fig.add_trace(go.Scatter3d(
                x=inter_x, y=inter_y, z=inter_z, mode='lines',
                line=dict(color=color_code, width=5, dash='dash'),
                name=f"分子間結合: {labels_list[0]}", visible=True
            ))

    # Layout Updates with Dark Theme & Buttons
    fig.update_layout(
        title=f"🧪 Dimer Fiedler Vector 3D Map: {name} | Fiedler λ2 = {fiedler_val:.6f}",
        scene=dict(
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            zaxis=dict(visible=False),
            bgcolor='rgb(15, 23, 42)' # 高級感のあるダークスレート背景
        ),
        margin=dict(l=0, r=0, b=0, t=80)
    )

    out_path = f"/home/eldenring/waterMain/{out_filename}"
    fig.write_html(out_path)
    print(f"🎉 3D Fiedler ダイマー視覚化 HTML 出力完了: {out_path}")
    return out_path

if __name__ == "__main__":
    generate_dimer_fiedler_3d_html("チロシン (Tyr)", "NC(Cc1ccc(O)cc1)C(=O)O", "dimer_fiedler_3d_tyrosine.html")
    generate_dimer_fiedler_3d_html("テレフタル酸", "O=C(O)c1ccc(C(=O)O)cc1", "dimer_fiedler_3d_terephthalic.html")
    generate_dimer_fiedler_3d_html("1-ブタノール", "CCCCO", "dimer_fiedler_3d_1_butanol.html")
