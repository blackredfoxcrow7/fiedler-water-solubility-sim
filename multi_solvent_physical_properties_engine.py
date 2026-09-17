"""
multi_solvent_physical_properties_engine.py
=============================================================================
【溶媒構造データ・誘電率・分極率組み込み型 汎用マルチソルベント Fiedler エンジン】

【物理化学的極性不適合ルール (Polar Mismatch Penalty Rule)】
・重原子に占める極性原子比率 (num_polar / num_heavy_atoms) を算出。
1. 「強親水性分子 (Glucose, Glycine, Urea)」:
   極性重原子割合が高く (Polar Ratio >= 30%), 無極性ヘキサン中では溶媒和不可能 (λ2 = 0)
2. 「強疎水性分子 (Naphthalene, Octane, Anthracene)」:
   無極性炭化水素であり, 水中では水素結合不適合により溶媒和不可能 (λ2 = 0), ヘキサン中でのみ高溶媒和
3. 「両親媒性分子 (Caffeine, Benzoic Acid, 1-Butanol)」:
   極性基と疎水基を両方持ち, アルコール・DMSO・アセトン等の各種極性/無極性溶媒にバランス溶解
=============================================================================
"""

import numpy as np
import networkx as nx
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

PHYSICAL_SOLVENT_DATABASE = {
    "water": {
        "name": "水 (Water)",
        "smiles": "O",
        "dielectric": 80.1,
        "dipole_debye": 1.85,
        "molar_volume_ml_mol": 18.07,
        "type": "polar_protic"
    },
    "methanol": {
        "name": "メタノール (MeOH)",
        "smiles": "CO",
        "dielectric": 32.7,
        "dipole_debye": 1.70,
        "molar_volume_ml_mol": 40.45,
        "type": "polar_protic"
    },
    "ethanol": {
        "name": "エタノール (EtOH)",
        "smiles": "CCO",
        "dielectric": 24.5,
        "dipole_debye": 1.69,
        "molar_volume_ml_mol": 58.39,
        "type": "amphiphilic"
    },
    "dmso": {
        "name": "DMSO (ジメチルスルホキシド)",
        "smiles": "CS(=O)C",
        "dielectric": 46.7,
        "dipole_debye": 3.96,
        "molar_volume_ml_mol": 71.03,
        "type": "polar_aprotic"
    },
    "dmf": {
        "name": "DMF (ジメチルホルムアミド)",
        "smiles": "CN(C)C=O",
        "dielectric": 36.7,
        "dipole_debye": 3.82,
        "molar_volume_ml_mol": 77.43,
        "type": "polar_aprotic"
    },
    "acetone": {
        "name": "アセトン (Acetone)",
        "smiles": "CC(=O)C",
        "dielectric": 20.7,
        "dipole_debye": 2.88,
        "molar_volume_ml_mol": 74.08,
        "type": "polar_aprotic"
    },
    "thf": {
        "name": "THF (テトラヒドロフラン)",
        "smiles": "C1CCCO1",
        "dielectric": 7.58,
        "dipole_debye": 1.75,
        "molar_volume_ml_mol": 81.11,
        "type": "ether"
    },
    "hexane": {
        "name": "n-ヘキサン (Hexane)",
        "smiles": "CCCCCC",
        "dielectric": 1.89,
        "dipole_debye": 0.08,
        "molar_volume_ml_mol": 131.57,
        "type": "non_polar"
    }
}

def onsager_reaction_field_factor(dielectric):
    return (dielectric - 1.0) / (2.0 * dielectric + 1.0)

def build_multi_solvent_physical_graph(solute_smiles, solvent_key="water", num_solvents=12):
    solv_info = PHYSICAL_SOLVENT_DATABASE[solvent_key]
    mol_solute = Chem.AddHs(Chem.MolFromSmiles(solute_smiles))
    N_solute = mol_solute.GetNumAtoms()
    
    # 重原子数および極性重原子数
    heavy_atoms = [a for a in mol_solute.GetAtoms() if a.GetSymbol() != 'H']
    num_heavy = len(heavy_atoms)
    num_polar = sum(1 for a in heavy_atoms if a.GetSymbol() in ['O', 'N', 'S', 'P', 'F', 'Cl'])
    polar_ratio = num_polar / num_heavy if num_heavy > 0 else 0.0
    
    # 極性不適合ルール (Like dissolves like)
    # 1. 強極性親水分子 (polar_ratio >= 0.30) ➔ 無極性ヘキサン中では溶媒和不可能 (λ2 = 0)
    if polar_ratio >= 0.30 and solv_info["type"] == "non_polar":
        return nx.Graph() # 非連結グラフ (λ2 = 0)
        
    # 2. 完全無極性分子 (polar_ratio == 0.0) ➔ 極性プロトン性溶媒 (Water, MeOH) では溶媒和不可能 (λ2 = 0)
    if polar_ratio == 0.0 and solv_info["type"] == "polar_protic":
        return nx.Graph() # 非連結グラフ (λ2 = 0)
        
    G = nx.Graph()
    for atom in mol_solute.GetAtoms():
        idx = atom.GetIdx()
        sym = atom.GetSymbol()
        is_polar = (sym in ['O', 'N', 'S', 'P', 'F', 'Cl'])
        G.add_node(f"S_{idx}", symbol=sym, is_polar=is_polar, type='solute')
        
    for bond in mol_solute.GetBonds():
        G.add_edge(f"S_{bond.GetBeginAtomIdx()}", f"S_{bond.GetEndAtomIdx()}", weight=1.0)

    mol_solv = Chem.AddHs(Chem.MolFromSmiles(solv_info["smiles"]))
    f_onsager = onsager_reaction_field_factor(solv_info["dielectric"])
    dipole_ratio = max(0.2, solv_info["dipole_debye"] / 1.85)
    
    w_solvent_net = max(0.4, 0.6 * f_onsager)
    w_solvation_base = 1.2 * (f_onsager ** 0.5) * dipole_ratio
    
    solv_centers = []
    for s_i in range(num_solvents):
        prefix = f"V_{s_i}"
        for atom in mol_solv.GetAtoms():
            a_idx = atom.GetIdx()
            sym = atom.GetSymbol()
            node_id = f"{prefix}_{a_idx}"
            is_polar = (sym in ['O', 'N', 'S'])
            G.add_node(node_id, symbol=sym, is_polar=is_polar, type='solvent')
            
            if solv_info["type"] != "non_polar":
                if is_polar and (f"{prefix}_center" not in solv_centers):
                    solv_centers.append(node_id)
            else:
                if sym == 'C' and (f"{prefix}_center" not in solv_centers):
                    solv_centers.append(node_id)
                
        for bond in mol_solv.GetBonds():
            G.add_edge(f"{prefix}_{bond.GetBeginAtomIdx()}", f"{prefix}_{bond.GetEndAtomIdx()}", weight=1.0)

    if len(solv_centers) > 1:
        for i in range(len(solv_centers) - 1):
            G.add_edge(solv_centers[i], solv_centers[i+1], weight=w_solvent_net)
        if len(solv_centers) > 2:
            G.add_edge(solv_centers[-1], solv_centers[0], weight=w_solvent_net)

    for s_i in range(N_solute):
        s_node = f"S_{s_i}"
        s_polar = G.nodes[s_node]["is_polar"]
        
        for v_c in solv_centers:
            if s_polar and solv_info["type"] != "non_polar":
                G.add_edge(s_node, v_c, weight=w_solvation_base)
            elif not s_polar and solv_info["type"] == "non_polar":
                G.add_edge(s_node, v_c, weight=0.8)

    return G

def calculate_multi_solvent_fiedler(solute_smiles, solvent_key="water"):
    G = build_multi_solvent_physical_graph(solute_smiles, solvent_key=solvent_key)
    if G.number_of_nodes() == 0 or not nx.is_connected(G):
        return 0.0
    L_norm = nx.normalized_laplacian_matrix(G).toarray()
    evals = np.linalg.eigvalsh(L_norm)
    return float(evals[1])

def run_multi_solvent_benchmark(name, smiles):
    print(f"\n🧪 【{name}】の全 8 有機溶媒物性データベース Fiedler 溶媒和親和性比較")
    print("-" * 90)
    print(f"{'溶媒名':<25} | {'溶媒分類':<15} | {'誘電率 εr':<10} | {'ダイポール(D)':<12} | {'Fiedler λ2'}")
    print("-" * 90)
    res = []
    for key, info in PHYSICAL_SOLVENT_DATABASE.items():
        f_val = calculate_multi_solvent_fiedler(smiles, solvent_key=key)
        res.append({
            "solvent": info["name"],
            "type": info["type"],
            "dielectric": info["dielectric"],
            "dipole": info["dipole_debye"],
            "fiedler_lambda2": f_val
        })
        print(f"{info['name']:<25} | {info['type']:<15} | {info['dielectric']:<10.1f} | {info['dipole_debye']:<12.2f} | {f_val:.6f}")
    print("-" * 90)
    return pd.DataFrame(res)

if __name__ == "__main__":
    run_multi_solvent_benchmark("チロシン (Tyr)", "NC(Cc1ccc(O)cc1)C(=O)O")
    run_multi_solvent_benchmark("テレフタル酸", "O=C(O)c1ccc(C(=O)O)cc1")
