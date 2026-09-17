"""
visualize_terephthalic_pi_stacking.py
=============================================================================
テレフタル酸のベンゼン環重なり（芳香環 π-π スタッキング）実空間 3D 視覚化

【結晶物理的特徴】
1. テレフタル酸（1,4-ベンゼンジカルボン酸）の実測結晶格子（Form I / Form II）では、
   1,4-カルボキシル基間の直線水素結合鎖に加え、
   隣り合うベンゼン環同士が面対向で平行重なり（π-πスタッキング: 面間距離 3.40 Å）を形成。
2. 本スクリプトは、実空間 3D デカルト座標 (Å 単位) において、
   重なる2つのベンゼン環の面間距離 (3.40 Å) および環中心距離 (3.61 Å) を
   リアルタイム 3D 数値表示付きで CPK / Ball-and-Stick 模型化。
3. 原子球上に Fiedler 固有ベクトルの符号（赤＝正 / 青＝負）を連続マッピング。
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

def generate_terephthalic_pi_stacking_3d(out_filename="terephthalic_pi_stacking_dimer.html"):
    smiles = "O=C(O)c1ccc(C(=O)O)cc1" # テレフタル酸
    mol, coords_raw = fdio.get_3d_conformer(smiles)
    N = mol.GetNumAtoms()
    
    # 元素の vdW 半径
    pt = Chem.GetPeriodicTable()
    vdw_radii = [pt.GetRvdw(mol.GetAtomWithIdx(i).GetSymbol()) for i in range(N)]

    # ベンゼン環炭素の原子インデックスを抽出
    ring_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIsAromatic()]
    
    # Molecule A: ベンゼン環中心を原点に配置
    coords_A = coords_raw - np.mean(coords_raw[ring_atoms], axis=0)
    
    # ベンゼン環 A の法線ベクトル（面直方向）を SVD から計算
    _, _, Vt = np.linalg.svd(coords_A[ring_atoms])
    normal_A = Vt[2, :] # 環平面の法線ベクトル
    inplane_u = Vt[0, :] # 環平面内の主軸ベクトル

    # Molecule B: 面間距離 3.40 Å (π-πスタッキングの標準距離) かつ ズレ重なり (Displaced Stacking) で配置
    interplanar_dist = 3.40 # 面直離隔距離 (Å)
    displaced_shift = 1.20  # 面内ずらし距離 (Å)
    
    shift_vec = normal_A * interplanar_dist + inplane_u * displaced_shift
    coords_B = coords_A + shift_vec

    # 3D 位置辞書
    pos_3d = {}
    for i in range(N):
        pos_3d[f"A_{i}"] = coords_A[i]
        pos_3d[f"B_{i}"] = coords_B[i]

    # Fiedler 固有ベクトルの算出 (ベンゼン環重なりダイマー)
    G_dimer = nx.Graph()
    for i in range(N):
        sym = mol.GetAtomWithIdx(i).GetSymbol()
        G_dimer.add_node(f"A_{i}", symbol=sym, mol='A')
        G_dimer.add_node(f"B_{i}", symbol=sym, mol='B')
    for bond in mol.GetBonds():
        G_dimer.add_edge(f"A_{bond.GetBeginAtomIdx()}", f"A_{bond.GetEndAtomIdx()}", weight=1.0, edge_type='covalent')
        G_dimer.add_edge(f"B_{bond.GetBeginAtomIdx()}", f"B_{bond.GetEndAtomIdx()}", weight=1.0, edge_type='covalent')

    # ベンゼン環同士の π-π スタッキングエッジ (重み 1.2)
    for r_a in ring_atoms:
        for r_b in ring_atoms:
            dist_ab = np.linalg.norm(pos_3d[f"A_{r_a}"] - pos_3d[f"B_{r_b}"])
            if dist_ab < 4.2:
                G_dimer.add_edge(f"A_{r_a}", f"B_{r_b}", weight=1.2, edge_type='pi_stacking')

    # 端部カルボキシル基間の水素結合エッジ (重み 1.5)
    for i_a in range(N):
        for i_b in range(N):
            sym_a = mol.GetAtomWithIdx(i_a).GetSymbol()
            sym_b = mol.GetAtomWithIdx(i_b).GetSymbol()
            if sym_a == 'O' and sym_b == 'O':
                dist_o = np.linalg.norm(pos_3d[f"A_{i_a}"] - pos_3d[f"B_{i_b}"])
                if dist_o < 3.5:
                    G_dimer.add_edge(f"A_{i_a}", f"B_{i_b}", weight=1.5, edge_type='hbond')

    L_dimer = nx.normalized_laplacian_matrix(G_dimer).toarray()
    evals, evecs = np.linalg.eigh(L_dimer)
    fiedler_val = evals[1]
    fiedler_vec = evecs[:, 1]
    
    nodes_order = list(G_dimer.nodes())
    node_to_fiedler = {node: fiedler_vec[idx] for idx, node in enumerate(nodes_order)}
    max_abs_f = max(abs(v) for v in fiedler_vec) if len(fiedler_vec) > 0 else 1.0

    # Plotly 真の実空間 3D 描画
    fig = go.Figure()

    # Trace 0: Molecule A (ベンゼン環1)
    x_A, y_A, z_A, c_A, s_A, text_A = [], [], [], [], [], []
    for i in range(N):
        node = f"A_{i}"
        x, y, z = pos_3d[node]
        val = node_to_fiedler[node]
        sign_str = "(+ 正)" if val > 0 else "(- 負)"
        sym = mol.GetAtomWithIdx(i).GetSymbol()
        x_A.append(x); y_A.append(y); z_A.append(z)
        c_A.append(val)
        s_A.append(vdw_radii[i] * 12)
        text_A.append(f"<b>テレフタル酸 A - Atom {i} ({sym})</b><br>実空間座標: ({x:.3f}, {y:.3f}, {z:.3f}) Å<br>Fiedler値: {val:+.4f} {sign_str}")

    # Trace 1: Molecule B (重なるベンゼン環2)
    x_B, y_B, z_B, c_B, s_B, text_B = [], [], [], [], [], []
    for i in range(N):
        node = f"B_{i}"
        x, y, z = pos_3d[node]
        val = node_to_fiedler[node]
        sign_str = "(+ 正)" if val > 0 else "(- 負)"
        sym = mol.GetAtomWithIdx(i).GetSymbol()
        x_B.append(x); y_B.append(y); z_B.append(z)
        c_B.append(val)
        s_B.append(vdw_radii[i] * 12)
        text_B.append(f"<b>テレフタル酸 B - Atom {i} ({sym})</b><br>実空間座標: ({x:.3f}, {y:.3f}, {z:.3f}) Å<br>Fiedler値: {val:+.4f} {sign_str}")

    # Add Molecule A Nodes
    fig.add_trace(go.Scatter3d(
        x=x_A, y=y_A, z=z_A, mode='markers+text',
        marker=dict(
            size=s_A, color=c_A, colorscale='RdBu_r', cmin=-max_abs_f, cmax=max_abs_f,
            colorbar=dict(title="Fiedler Vector Value (符号)", x=1.05),
            line=dict(color='black', width=1.5)
        ),
        text=[f"A:{mol.GetAtomWithIdx(i).GetSymbol()}" for i in range(N)],
        hoverinfo='text', hovertext=text_A, name="テレフタル酸 A (ベンゼン環1)", visible=True
    ))

    # Add Molecule B Nodes
    fig.add_trace(go.Scatter3d(
        x=x_B, y=y_B, z=z_B, mode='markers+text',
        marker=dict(
            size=s_B, color=c_B, colorscale='RdBu_r', cmin=-max_abs_f, cmax=max_abs_f,
            line=dict(color='gold', width=2.5)
        ),
        text=[f"B:{mol.GetAtomWithIdx(i).GetSymbol()}" for i in range(N)],
        hoverinfo='text', hovertext=text_B, name="テレフタル酸 B (重なるベンゼン環2)", visible=True
    ))

    # Trace 2: 共有結合スティック A & B
    cov_x, cov_y, cov_z = [], [], []
    for u, v, d in G_dimer.edges(data=True):
        if d.get('edge_type') == 'covalent':
            x0, y0, z0 = pos_3d[u]
            x1, y1, z1 = pos_3d[v]
            cov_x.extend([x0, x1, None]); cov_y.extend([y0, y1, None]); cov_z.extend([z0, z1, None])

    fig.add_trace(go.Scatter3d(
        x=cov_x, y=cov_y, z=cov_z, mode='lines',
        line=dict(color='lightgray', width=7), name='実空間共有結合骨格', visible=True
    ))

    # Trace 3: ベンゼン環中心同士を結ぶ π-π スタッキング軸 ＋ 面間距離ラベル
    center_A = np.mean(coords_A[ring_atoms], axis=0)
    center_B = np.mean(coords_B[ring_atoms], axis=0)
    center_dist = np.linalg.norm(center_A - center_B)

    fig.add_trace(go.Scatter3d(
        x=[center_A[0], center_B[0]], y=[center_A[1], center_B[1]], z=[center_A[2], center_B[2]],
        mode='lines+text',
        line=dict(color='gold', width=9, dash='dash'),
        text=["", f"<b>芳香環 π-π 重なり面間距離: {interplanar_dist:.2f} Å (中心間: {center_dist:.2f} Å)</b>"],
        textposition="top center",
        textfont=dict(color='gold', size=14),
        name="ベンゼン環 π-π スタッキング軸 (3.40 Å)"
    ))

    fig.update_layout(
        title=f"🏛️ テレフタル酸 ベンゼン環重なり（π-πスタッキング）実空間3Dモデル | Fiedler λ2 = {fiedler_val:.6f}",
        scene=dict(
            aspectmode='data', # 1:1:1 実空間スケール
            xaxis=dict(title="X (Å)", backgroundcolor='rgb(15, 23, 42)', gridcolor='gray'),
            yaxis=dict(title="Y (Å)", backgroundcolor='rgb(15, 23, 42)', gridcolor='gray'),
            zaxis=dict(title="Z (Å)", backgroundcolor='rgb(15, 23, 42)', gridcolor='gray'),
            bgcolor='rgb(15, 23, 42)'
        ),
        margin=dict(l=0, r=0, b=0, t=80)
    )

    out_path = f"/home/eldenring/waterMain/{out_filename}"
    fig.write_html(out_path)
    print(f"🎉 テレフタル酸ベンゼン環重なり 3D HTML 出力完了: {out_path}")
    return out_path

if __name__ == "__main__":
    generate_terephthalic_pi_stacking_3d()
