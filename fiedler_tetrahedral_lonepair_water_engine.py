"""
fiedler_tetrahedral_lonepair_water_engine.py
=============================================================================
【水分子の孤立電子対 (Lone Pair) 四面体配向ベクトル組み込み型 Fiedler 水和エンジン】

【理論・アルゴリズム概要】
1. 水分子の酸素原子（O）が持つ 2 対の孤立電子対（Lone Pairs: LP1, LP2）を、
   O-H 結合（104.5°）と直交する tetrahedral（109.5°）四面体幾何ベクトルとして 3D 空間内に明示的に構築。
2. 溶質分子の -OH / -NH ドナー結合の 3D 方向ベクトル v_donor と、水分子の Lone Pair ベクトル v_acceptor
   との立体交差角 θ_align を計算。
3. グラフエッジ選択規則 (Edge Selection Rule) に cos^2(θ_align) 方向適合関数を組み込み：
   w_edge = w_base * max(0, cos(θ_align))^2
   完璧に四面体配向適合（0°）のときエッジ重みが最大化され、立体衝突（90°直交）ではエッジ重みが 0 となる。
=============================================================================
"""

import numpy as np
import networkx as nx
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

def compute_tetrahedral_lonepairs(o_pos, h1_pos, h2_pos):
    """
    水分子の O-H1, O-H2 3D 座標から、四面体幾何 (109.5°) に沿った
    2 つの孤立電子対 (LP1, LP2) の 3D 座標を幾何学的に自動算出する
    """
    v_oh1 = h1_pos - o_pos
    v_oh2 = h2_pos - o_pos
    
    v_oh1 /= np.linalg.norm(v_oh1)
    v_oh2 /= np.linalg.norm(v_oh2)
    
    # H-O-H 二等分線ベクトル (Bisector)
    v_bisect = v_oh1 + v_oh2
    v_bisect /= np.linalg.norm(v_bisect)
    
    # H-O-H 平面の法線ベクトル (Normal)
    v_normal = np.cross(v_oh1, v_oh2)
    v_normal /= np.linalg.norm(v_normal)
    
    # 四面体アングル: 二等分線の反転方向 + 法線方向への合成
    r_lp = 0.96 # LP までの実効距離 (Å)
    angle_half = np.radians(109.5 / 2.0)
    
    lp1_dir = -np.cos(angle_half) * v_bisect + np.sin(angle_half) * v_normal
    lp2_dir = -np.cos(angle_half) * v_bisect - np.sin(angle_half) * v_normal
    
    lp1_pos = o_pos + r_lp * lp1_dir
    lp2_pos = o_pos + r_lp * lp2_dir
    
    return lp1_pos, lp2_pos

def build_tetrahedral_lonepair_water_system(solute_smiles, num_waters=12):
    """
    孤立電子対四面体ベクトル配向則をグラフエッジ選択に適用した統合ネットワークを構築
    """
    mol = Chem.MolFromSmiles(solute_smiles)
    if mol is None:
        raise ValueError(f"SMILESパースエラー: {solute_smiles}")
        
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    AllChem.MMFFOptimizeMolecule(mol)
    conf = mol.GetConformer()
    
    N_solute = mol.GetNumAtoms()
    coords_solute = conf.GetPositions()
    
    G = nx.Graph()
    
    # 1. 溶質ノード追加
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        sym = atom.GetSymbol()
        G.add_node(f"S_{idx}", symbol=sym, is_polar=(sym in ['O','N']), pos=coords_solute[idx])
        
    for b in mol.GetBonds():
        G.add_edge(f"S_{b.GetBeginAtomIdx()}", f"S_{b.GetEndAtomIdx()}", weight=1.0, type='solute_covalent')

    # 2. 水分子群の構築（四面体 Lone Pair 付与）
    centroid = np.mean(coords_solute, axis=0)
    radius = 4.5
    
    for w_i in range(num_waters):
        phi = np.arccos(1 - 2 * (w_i + 0.5) / num_waters)
        theta = np.pi * (1 + 5**0.5) * w_i
        
        o_pos = centroid + radius * np.array([np.sin(phi)*np.cos(theta), np.sin(phi)*np.sin(theta), np.cos(phi)])
        
        v_in = centroid - o_pos
        v_in /= np.linalg.norm(v_in)
        
        v_perpendicular = np.array([-v_in[1], v_in[0], 0.0])
        if np.linalg.norm(v_perpendicular) < 1e-3:
            v_perpendicular = np.array([1.0, 0.0, 0.0])
        v_perpendicular /= np.linalg.norm(v_perpendicular)
        
        h1_pos = o_pos + 0.96 * (v_in * np.cos(np.radians(52.25)) + v_perpendicular * np.sin(np.radians(52.25)))
        h2_pos = o_pos + 0.96 * (v_in * np.cos(np.radians(52.25)) - v_perpendicular * np.sin(np.radians(52.25)))
        
        lp1_pos, lp2_pos = compute_tetrahedral_lonepairs(o_pos, h1_pos, h2_pos)
        
        w_o_id = f"W_{w_i}_O"
        w_lp1_id = f"W_{w_i}_LP1"
        w_lp2_id = f"W_{w_i}_LP2"
        
        G.add_node(w_o_id, symbol='O', pos=o_pos, type='water')
        G.add_node(w_lp1_id, symbol='LP', pos=lp1_pos, type='lonepair')
        G.add_node(w_lp2_id, symbol='LP', pos=lp2_pos, type='lonepair')
        
        G.add_edge(w_o_id, w_lp1_id, weight=0.8, type='water_internal')
        G.add_edge(w_o_id, w_lp2_id, weight=0.8, type='water_internal')
        
        # 3. エッジ選択規則: 孤立電子対ベクトル配向適合度 cos^2(θ_align) による重み付け
        for s_idx in range(N_solute):
            s_atom = mol.GetAtomWithIdx(s_idx)
            if s_atom.GetSymbol() == 'O':
                h_neighbors = [nbr for nbr in s_atom.GetNeighbors() if nbr.GetSymbol() == 'H']
                if h_neighbors:
                    h_idx = h_neighbors[0].GetIdx()
                    v_donor = coords_solute[h_idx] - coords_solute[s_idx]
                    v_donor_norm = v_donor / np.linalg.norm(v_donor)
                    
                    v_acceptor1 = lp1_pos - o_pos
                    v_acceptor1_norm = v_acceptor1 / np.linalg.norm(v_acceptor1)
                    
                    cos_theta1 = -np.dot(v_donor_norm, v_acceptor1_norm)
                    
                    if cos_theta1 > 0:
                        edge_w = 1.5 * (cos_theta1 ** 2) # 方向適合規則
                        G.add_edge(f"S_{s_idx}", w_lp1_id, weight=edge_w, type='hbond_lonepair_aligned')
                        
    return G

def calculate_tetrahedral_fiedler(solute_smiles):
    """孤立電子対四面体配向則を適用した Fiedler 値 λ2 を算出"""
    G = build_tetrahedral_lonepair_water_system(solute_smiles)
    if not nx.is_connected(G):
        return 0.0
    L_norm = nx.normalized_laplacian_matrix(G).toarray()
    evals = np.linalg.eigvalsh(L_norm)
    return float(evals[1])

if __name__ == "__main__":
    print("=========================================================================")
    print(" 🧪 水分子孤立電子対 (Lone Pair) 四面体配向エンジン デモ検証")
    print("=========================================================================")
    dataset = [
        ("D-マンニトール (難溶性)", "OCC(O)C(O)C(O)C(O)CO"),
        ("グルコース (ブドウ糖)", "OCC1OC(O)C(O)C(O)C1O"),
        ("フルクトース (果糖)", "OCC1OC(O)(CO)C(O)C1O"),
        ("スクロース (ショ糖)", "C(C1C(C(C(C(O1)OC2(C(C(C(O2)CO)O)O)CO)O)O)O)O")
    ]
    for name, smiles in dataset:
        f_val = calculate_tetrahedral_fiedler(smiles)
        print(f"  ・{name:<25} : 孤立電子対四面体 Fiedler λ2 = {f_val:.6f}")
