"""
multi_solvent_fiedler_engine.py
=============================================================================
【汎用多種溶媒・物理ネットワーク Fiedler シミュレーションエンジン】

【目的】
・水 (H2O) だけでなく、有機溶媒（エタノール, メタノール, DMSO, DMF, アセトン, THF, ヘキサン等）
  における溶質分子のミクロ溶媒和親和性および溶解性を統一的に計算・予測する。

【溶媒物理パラメータモデル】
・誘電率 (Dielectric Constant, εr)
・分子量 (Mw) および モル体積 (Vm = Mw / ρ)
・水素結合ドナー数 (H-bond Donor Capacity)
・水素結合アクセプター数 (H-bond Acceptor Capacity)
・両親媒性 / 疎水性アルキル鎖パラメータ
=============================================================================
"""

import sys
import os
import networkx as nx
import numpy as np
import pandas as pd
from rdkit import Chem

# 汎用溶媒物理データベース
SOLVENT_DATABASE = {
    "water": {
        "name": "水 (Water)",
        "smiles": "O",
        "dielectric": 80.1,
        "molar_volume": 18.07,
        "donor_cap": 2,
        "acceptor_cap": 2,
        "type": "polar_protic"
    },
    "methanol": {
        "name": "メタノール (MeOH)",
        "smiles": "CO",
        "dielectric": 32.7,
        "molar_volume": 40.45,
        "donor_cap": 1,
        "acceptor_cap": 2,
        "type": "polar_protic"
    },
    "ethanol": {
        "name": "エタノール (EtOH)",
        "smiles": "CCO",
        "dielectric": 24.5,
        "molar_volume": 58.39,
        "donor_cap": 1,
        "acceptor_cap": 2,
        "type": "amphiphilic"
    },
    "dmso": {
        "name": "DMSO (ジメチルスルホキシド)",
        "smiles": "CS(=O)C",
        "dielectric": 46.7,
        "molar_volume": 71.03,
        "donor_cap": 0,
        "acceptor_cap": 2,
        "type": "polar_aprotic"
    },
    "dmf": {
        "name": "DMF (ジメチルホルムアミド)",
        "smiles": "CN(C)C=O",
        "dielectric": 36.7,
        "molar_volume": 77.43,
        "donor_cap": 0,
        "acceptor_cap": 1,
        "type": "polar_aprotic"
    },
    "acetone": {
        "name": "アセトン (Acetone)",
        "smiles": "CC(=O)C",
        "dielectric": 20.7,
        "molar_volume": 74.08,
        "donor_cap": 0,
        "acceptor_cap": 1,
        "type": "polar_aprotic"
    },
    "hexane": {
        "name": "n-ヘキサン (Hexane)",
        "smiles": "CCCCCC",
        "dielectric": 1.89,
        "molar_volume": 131.57,
        "donor_cap": 0,
        "acceptor_cap": 0,
        "type": "non_polar"
    }
}

def create_multi_solvent_integrated_system(solute_smiles, solvent_key="water", solvent_count=12):
    """任意溶媒中での溶質分子統合ネットワークグラフを構築"""
    if solvent_key not in SOLVENT_DATABASE:
        raise ValueError(f"未対応の溶媒キーです: {solvent_key}")
        
    solv_info = SOLVENT_DATABASE[solvent_key]
    G = nx.Graph()
    
    # 1. 溶質分子ノードと結合の構築
    mol_solute = Chem.AddHs(Chem.MolFromSmiles(solute_smiles))
    N_solute = mol_solute.GetNumAtoms()
    
    for atom in mol_solute.GetAtoms():
        idx, sym = atom.GetIdx(), atom.GetSymbol()
        node_id = f"S_{idx}"
        h_count = sum(1 for nbr in atom.GetNeighbors() if nbr.GetSymbol() == 'H')
        is_polar = sym in ['O', 'N', 'S', 'P', 'F', 'Cl', 'Br'] or atom.GetFormalCharge() != 0
        G.add_node(node_id, symbol=sym, type='solute', is_polar=is_polar, h_count=h_count)
        
    for bond in mol_solute.GetBonds():
        G.add_edge(f"S_{bond.GetBeginAtomIdx()}", f"S_{bond.GetEndAtomIdx()}", weight=1.0, edge_type='covalent_solute')
        
    # 2. 溶媒分子ノードとネットワークの構築
    solvent_smiles = solv_info["smiles"]
    mol_solv = Chem.AddHs(Chem.MolFromSmiles(solvent_smiles))
    N_solv_atom = mol_solv.GetNumAtoms()
    
    solv_centers = []
    
    for s_idx in range(solvent_count):
        solv_prefix = f"V_{s_idx}"
        # 溶媒分子内部の描画
        for atom in mol_solv.GetAtoms():
            a_idx, sym = atom.GetIdx(), atom.GetSymbol()
            v_node = f"{solv_prefix}_{a_idx}"
            is_polar = sym in ['O', 'N', 'S']
            G.add_node(v_node, symbol=sym, type='solvent', solvent_id=s_idx, is_polar=is_polar)
            if is_polar and (f"{solv_prefix}_center" not in solv_centers):
                solv_centers.append(v_node)
                
        for bond in mol_solv.GetBonds():
            G.add_edge(f"{solv_prefix}_{bond.GetBeginAtomIdx()}", f"{solv_prefix}_{bond.GetEndAtomIdx()}", weight=1.0, edge_type='covalent_solvent')

    # 3. 溶媒同士の物理ネットワーク (Solvent-Solvent Network)
    # 誘電率 εr および モル体積 Vm に依存するエッジ架橋
    solv_net_weight = 0.5 * (solv_info["dielectric"] / 80.1) # 水を 0.5 とする規格化重み
    
    if len(solv_centers) > 1:
        for i in range(len(solv_centers) - 1):
            G.add_edge(solv_centers[i], solv_centers[i+1], weight=solv_net_weight, edge_type='solvent_net')
        if len(solv_centers) > 2:
            G.add_edge(solv_centers[-1], solv_centers[0], weight=solv_net_weight, edge_type='solvent_net')

    # 4. 溶質-溶媒間の溶媒和相互作用エッジ (Solute-Solvent Solvation Edges)
    for s_idx in range(N_solute):
        s_node = f"S_{s_idx}"
        s_sym = G.nodes[s_node]["symbol"]
        s_polar = G.nodes[s_node]["is_polar"]
        
        for v_center in solv_centers:
            v_sym = G.nodes[v_center]["symbol"]
            
            # H結合親和性
            if s_polar and solv_info["type"] in ["polar_protic", "amphiphilic", "polar_aprotic"]:
                # 誘電率と親水性を反映したエッジ重み
                edge_w = 1.2 * (solv_info["dielectric"] / 80.1) ** 0.5
                G.add_edge(s_node, v_center, weight=edge_w, edge_type='solvation_hbond')
            elif not s_polar and solv_info["type"] == "non_polar":
                # 無極性・疎水性溶媒和 (Hexane等)
                edge_w = 0.8
                G.add_edge(s_node, v_center, weight=edge_w, edge_type='solvation_vdw')

    return G

def compute_multi_solvent_fiedler(solute_smiles, solvent_key="water"):
    """指定された溶媒中での溶質の Fiedler 溶媒和値を算出"""
    G = create_multi_solvent_integrated_system(solute_smiles, solvent_key=solvent_key)
    
    if not nx.is_connected(G):
        return 0.0
        
    L_norm = nx.normalized_laplacian_matrix(G).toarray()
    evals = np.linalg.eigvalsh(L_norm)
    fiedler_val = float(evals[1])
    return fiedler_val

def compare_solute_across_solvents(name, smiles):
    print(f"\n🧪 【{name}】の各種有機溶媒中での Fiedler 溶媒和親和性比較")
    print("-" * 75)
    print(f"{'溶媒名':<25} | {'溶媒タイプ':<16} | {'誘電率 εr':<10} | {'Fiedler 溶媒和親和性 λ2'}")
    print("-" * 75)
    
    res = []
    for key, info in SOLVENT_DATABASE.items():
        f_val = compute_multi_solvent_fiedler(smiles, solvent_key=key)
        res.append({
            "solvent_key": key,
            "solvent_name": info["name"],
            "type": info["type"],
            "dielectric": info["dielectric"],
            "fiedler_lambda2": f_val
        })
        print(f"{info['name']:<25} | {info['type']:<16} | {info['dielectric']:<10.1f} | {f_val:.6f}")
    print("-" * 75)
    return pd.DataFrame(res)

if __name__ == "__main__":
    compare_solute_across_solvents("チロシン (Tyr)", "NC(Cc1ccc(O)cc1)C(=O)O")
    compare_solute_across_solvents("テレフタル酸", "O=C(O)c1ccc(C(=O)O)cc1")
