"""
visualize_real_space_dimer_3d.py
=============================================================================
実空間 3D 座標 (Å 単位) における分子模型表現・リアル原子間距離のダイマー視覚化

【物理・幾何機能】
1. RDKit MMFF94 力場最適化による実空間 3D デカルト座標 (Å 単位) をそのまま採用
   （グラフレイアウトやバネ埋め込みの近似は一切排除、1.54 Å 共有結合長等を精密反映）
2. Kabsch SVD 剛体変換による Molecule B の実空間 3D 回転・平行移動
3. 分子間相互作用（塩橋・水素結合・π-πスタッキング・疎水結合）の
   実空間距離 d (Å) を 3D 軸上に数値ラベル付きでリアルタイム表示
4. Fiedler 固有ベクトルの符号（正＝赤 / 負＝青）と数値を原子球上にマッピング
5. Plotly の `scene.aspectmode = 'data'` により、縦横高さが 1:1:1 の真の物理空間スケール描画
=============================================================================
"""

import sys
import os
import networkx as nx
import numpy as np
import plotly.graph_objects as go
from rdkit import Chem
from rdkit.Chem import AllChem

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import fiedler_dimer_interaction_optimizer as fdio

def generate_real_space_dimer_html(name, smiles, out_filename):
    mol, coords_raw = fdio.get_3d_conformer(smiles)
    N = mol.GetNumAtoms()
    
    # 元素ごとの van der Waals 半径 (Å)
    pt = Chem.GetPeriodicTable()
    vdw_radii = {i: pt.GetRvdw(mol.GetAtomWithIdx(i).GetSymbol()) for i in range(N)}

    # 1. 2分子ダイマーグラフ構築と Fiedler 分析
    G_dimer = nx.Graph()
    for i in range(N):
        sym = mol.GetAtomWithIdx(i).GetSymbol()
        G_dimer.add_node(f"A_{i}", symbol=sym, mol='A', orig_idx=i)
        G_dimer.add_node(f"B_{i}", symbol=sym, mol='B', orig_idx=i)
        
    for bond in mol.GetBonds():
        G_dimer.add_edge(f"A_{bond.GetBeginAtomIdx()}", f"A_{bond.GetEndAtomIdx()}", weight=1.0, edge_type='covalent')
        G_dimer.add_edge(f"B_{bond.GetBeginAtomIdx()}", f"B_{bond.GetEndAtomIdx()}", weight=1.0, edge_type='covalent')

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

    L_dimer = nx.normalized_laplacian_matrix(G_dimer).toarray()
    evals, evecs = np.linalg.eigh(L_dimer)
    fiedler_val = evals[1]
    fiedler_vec = evecs[:, 1]
    
    nodes_order = list(G_dimer.nodes())
    node_to_fiedler = {node: fiedler_vec[idx] for idx, node in enumerate(nodes_order)}
    
    # 2. 実空間 3D デカルト座標における Kabsch 最短変換
    coords_A = coords_raw.copy()
    coords_B = coords_raw.copy()
    
    if len(selected_edges) > 0:
        P = np.array([coords_A[int(u.split('_')[1])] for u, v, w, etype, label in selected_edges])
        Q = np.array([coords_B[int(v.split('_')[1])] for u, v, w, etype, label in selected_edges])
        W = np.array([w for u, v, w, etype, label in selected_edges])
        R, t = fdio.kabsch_rigid_transform(P, Q, weights=W)
        coords_B_opt = np.dot(coords_B, R.T) + t
    else:
        coords_B_opt = coords_B + np.array([4.0, 0.0, 0.0])

    pos_3d = {}
    for i in range(N):
        pos_3d[f"A_{i}"] = coords_A[i]
        pos_3d[f"B_{i}"] = coords_B_opt[i]

    # 3. Plotly 真の 3D 実空間スケール描画
    fig = go.Figure()
    max_abs_f = max(abs(v) for v in fiedler_vec) if len(fiedler_vec) > 0 else 1.0

    # Trace 0: Molecule A (実空間 Å 座標 + vdWサイズ)
    x_A, y_A, z_A, c_A, s_A, hover_A = [], [], [], [], [], []
    for i in range(N):
        node = f"A_{i}"
        x, y, z = pos_3d[node]
        val = node_to_fiedler[node]
        sign_str = "(+ 正)" if val > 0 else "(- 負)"
        sym = mol.GetAtomWithIdx(i).GetSymbol()
        x_A.append(x); y_A.append(y); z_A.append(z)
        c_A.append(val)
        s_A.append(vdw_radii[i] * 12) # vdW比例サイズ
        hover_A.append(f"<b>Molecule A - Atom {i} ({sym})</b><br>実空間座標: ({x:.3f}, {y:.3f}, {z:.3f}) Å<br>Fiedler値: {val:+.4f} {sign_str}<br>vdW半径: {vdw_radii[i]:.2f} Å")

    # Trace 1: Molecule B
    x_B, y_B, z_B, c_B, s_B, hover_B = [], [], [], [], [], []
    for i in range(N):
        node = f"B_{i}"
        x, y, z = pos_3d[node]
        val = node_to_fiedler[node]
        sign_str = "(+ 正)" if val > 0 else "(- 負)"
        sym = mol.GetAtomWithIdx(i).GetSymbol()
        x_B.append(x); y_B.append(y); z_B.append(z)
        c_B.append(val)
        s_B.append(vdw_radii[i] * 12)
        hover_B.append(f"<b>Molecule B - Atom {i} ({sym})</b><br>実空間座標: ({x:.3f}, {y:.3f}, {z:.3f}) Å<br>Fiedler値: {val:+.4f} {sign_str}<br>vdW半径: {vdw_radii[i]:.2f} Å")

    # Add Solute A Nodes
    fig.add_trace(go.Scatter3d(
        x=x_A, y=y_A, z=z_A, mode='markers+text',
        marker=dict(
            size=s_A, color=c_A, colorscale='RdBu_r', cmin=-max_abs_f, cmax=max_abs_f,
            colorbar=dict(title="Fiedler Vector Value (符号)", x=1.05),
            line=dict(color='black', width=1.5)
        ),
        text=[f"A:{mol.GetAtomWithIdx(i).GetSymbol()}" for i in range(N)],
        hoverinfo='text', hovertext=hover_A, name="Molecule A (実空間座標)", visible=True
    ))

    # Add Solute B Nodes
    fig.add_trace(go.Scatter3d(
        x=x_B, y=y_B, z=z_B, mode='markers+text',
        marker=dict(
            size=s_B, color=c_B, colorscale='RdBu_r', cmin=-max_abs_f, cmax=max_abs_f,
            line=dict(color='gold', width=2.0)
        ),
        text=[f"B:{mol.GetAtomWithIdx(i).GetSymbol()}" for i in range(N)],
        hoverinfo='text', hovertext=hover_B, name="Molecule B (Kabsch 3D最適化後)", visible=True
    ))

    # Trace 2: 実空間 3D 共有結合スティック (Covalent Sticks)
    cov_x, cov_y, cov_z = [], [], []
    for u, v, d in G_dimer.edges(data=True):
        if d.get('edge_type') == 'covalent':
            x0, y0, z0 = pos_3d[u]
            x1, y1, z1 = pos_3d[v]
            cov_x.extend([x0, x1, None]); cov_y.extend([y0, y1, None]); cov_z.extend([z0, z1, None])

    fig.add_trace(go.Scatter3d(
        x=cov_x, y=cov_y, z=cov_z, mode='lines',
        line=dict(color='gray', width=7), name='実空間共有結合スティック', visible=True
    ))

    # Trace 3: 実空間 3D 分子間相互作用エッジ ＋ 結合距離 (Å) 数値表示
    edge_colors = {
        'ionic_bridge': 'magenta',
        'hbond': 'lime',
        'pi_stacking': 'yellow',
        'hydrophobic': 'orange'
    }

    for etype_key, color_code in edge_colors.items():
        inter_x, inter_y, inter_z = [], [], []
        mid_x, mid_y, mid_z, mid_text = [], [], [], []
        
        for u, v, d in G_dimer.edges(data=True):
            if d.get('edge_type') == etype_key:
                p0 = pos_3d[u]
                p1 = pos_3d[v]
                dist_angstrom = np.linalg.norm(p0 - p1)
                
                inter_x.extend([p0[0], p1[0], None])
                inter_y.extend([p0[1], p1[1], None])
                inter_z.extend([p0[2], p1[2], None])
                
                # 中点に実空間距離 (Å) の数値テキストを配置
                mid_p = (p0 + p1) / 2.0
                mid_x.append(mid_p[0]); mid_y.append(mid_p[1]); mid_z.append(mid_p[2])
                mid_text.append(f"<b>{dist_angstrom:.2f} Å</b> ({d.get('label')})")

        if inter_x:
            # 破線
            fig.add_trace(go.Scatter3d(
                x=inter_x, y=inter_y, z=inter_z, mode='lines',
                line=dict(color=color_code, width=6, dash='dash'),
                name=f"分子間結合: {etype_key}", visible=True
            ))
            # 距離ラベル
            fig.add_trace(go.Scatter3d(
                x=mid_x, y=mid_y, z=mid_z, mode='text',
                text=mid_text, textposition="top center",
                textfont=dict(color=color_code, size=13),
                name=f"結合距離(Å): {etype_key}", visible=True
            ))

    # 4. 真の物理 1:1:1 アスペクト比の設定 (aspectmode='data')
    fig.update_layout(
        title=f"🧊 リアル 3D 実空間分子模型 (Å 単位): {name} | Fiedler λ2 = {fiedler_val:.6f}",
        scene=dict(
            aspectmode='data', # 1:1:1 物理実空間スケールを維持
            xaxis=dict(title="X (Å)", backgroundcolor='rgb(15, 23, 42)', gridcolor='gray'),
            yaxis=dict(title="Y (Å)", backgroundcolor='rgb(15, 23, 42)', gridcolor='gray'),
            zaxis=dict(title="Z (Å)", backgroundcolor='rgb(15, 23, 42)', gridcolor='gray'),
            bgcolor='rgb(15, 23, 42)'
        ),
        margin=dict(l=0, r=0, b=0, t=80)
    )

    out_path = f"/home/eldenring/waterMain/{out_filename}"
    fig.write_html(out_path)
    print(f"🎉 リアル実空間 3D HTML 出力完了: {out_path}")
    return out_path

if __name__ == "__main__":
    generate_real_space_dimer_html("チロシン (Tyr)", "NC(Cc1ccc(O)cc1)C(=O)O", "real_space_3d_tyrosine.html")
    generate_real_space_dimer_html("テレフタル酸", "O=C(O)c1ccc(C(=O)O)cc1", "real_space_3d_terephthalic.html")
    generate_real_space_dimer_html("1-ブタノール", "CCCCO", "real_space_3d_1_butanol.html")
