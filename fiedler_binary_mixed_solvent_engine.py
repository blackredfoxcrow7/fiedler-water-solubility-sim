"""
fiedler_binary_mixed_solvent_engine.py
=============================================================================
【二元混合溶媒（Water + Co-solvent）スペクトル Fiedler シミュレーションエンジン】

【理論モデル概要】
1. 水 (H2O) と水溶性共溶媒 (エタノール, DMSO, アセトン, メタノール等) のモル分率 χ_cosolvent (0.0 ~ 1.0)
   に応じた混合溶媒ネットワークグラフ G_mix(χ) を動的生成。
2. Jouyban-Acree 式に基づく混合誘電率 ε_mix(χ) およびオンサーガー反応場の連続補正:
   ln(ε_mix) = (1 - χ) * ln(ε_water) + χ * ln(ε_cosolv)
3. 溶質の選択的溶媒和 (Selective Solvation):
   - 親水性基 ➔ 水分子ネットワークと優先結合 (w_SW)
   - 疎水性基 ➔ 共溶媒のアルキル/メチル基と優先結合 (w_SV)
4. 共溶媒濃度 χ の変化に伴う Fiedler 溶媒和値 λ2(χ) のプロットにより、
   最適溶解度を与える「共溶媒極大ピーク (Cosolvency Peak)」を自動検出する。
=============================================================================
"""

import numpy as np
import networkx as nx
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

# 共溶媒物性パラメータ
COSOLVENT_DATABASE = {
    "ethanol": {
        "name": "水 ＋ エタノール (Water + EtOH)",
        "smiles": "CCO",
        "dielectric": 24.5,
        "dipole": 1.69,
        "type": "amphiphilic"
    },
    "dmso": {
        "name": "水 ＋ DMSO (Water + DMSO)",
        "smiles": "CS(=O)C",
        "dielectric": 46.7,
        "dipole": 3.96,
        "type": "polar_aprotic"
    },
    "acetone": {
        "name": "水 ＋ アセトン (Water + Acetone)",
        "smiles": "CC(=O)C",
        "dielectric": 20.7,
        "dipole": 2.88,
        "type": "polar_aprotic"
    },
    "methanol": {
        "name": "水 ＋ メタノール (Water + MeOH)",
        "smiles": "CO",
        "dielectric": 32.7,
        "dipole": 1.70,
        "type": "polar_protic"
    }
}

WATER_DIELECTRIC = 80.1
WATER_DIPOLE = 1.85

def compute_mixed_dielectric(chi_cosolv, cosolv_dielectric):
    """Jouyban-Acree 混合誘電率"""
    ln_e = (1.0 - chi_cosolv) * np.log(WATER_DIELECTRIC) + chi_cosolv * np.log(cosolv_dielectric)
    return float(np.exp(ln_e))

def build_binary_mixed_solvent_graph(solute_smiles, cosolv_key="ethanol", chi_cosolv=0.5, total_solvent_mols=20):
    """
    モル分率 chi_cosolv における 溶質-水-共溶媒 統合ネットワークグラフを構築
    """
    if cosolv_key not in COSOLVENT_DATABASE:
        raise ValueError(f"未対応の共溶媒キーです: {cosolv_key}")
        
    cosolv_info = COSOLVENT_DATABASE[cosolv_key]
    
    # 溶質構造
    mol_solute = Chem.AddHs(Chem.MolFromSmiles(solute_smiles))
    N_solute = mol_solute.GetNumAtoms()
    
    G = nx.Graph()
    for atom in mol_solute.GetAtoms():
        idx = atom.GetIdx()
        sym = atom.GetSymbol()
        is_polar = sym in ['O', 'N', 'S', 'P', 'F', 'Cl']
        G.add_node(f"S_{idx}", symbol=sym, is_polar=is_polar, type='solute')
        
    for b in mol_solute.GetBonds():
        G.add_edge(f"S_{b.GetBeginAtomIdx()}", f"S_{b.GetEndAtomIdx()}", weight=1.0)
        
    # 混合比率に基づく溶媒分子数の分配
    num_cosolv_mols = int(np.round(total_solvent_mols * chi_cosolv))
    num_water_mols = total_solvent_mols - num_cosolv_mols
    
    e_mix = compute_mixed_dielectric(chi_cosolv, cosolv_info["dielectric"])
    f_onsager = (e_mix - 1.0) / (2.0 * e_mix + 1.0)
    
    solvent_centers = []
    
    # 1. 水分子ノードの追加
    mol_w = Chem.AddHs(Chem.MolFromSmiles("O"))
    for w_i in range(num_water_mols):
        prefix = f"W_{w_i}"
        for a in mol_w.GetAtoms():
            node_id = f"{prefix}_{a.GetIdx()}"
            G.add_node(node_id, symbol=a.GetSymbol(), type='water')
            if a.GetSymbol() == 'O':
                solvent_centers.append((node_id, 'water'))
        for b in mol_w.GetBonds():
            G.add_edge(f"{prefix}_{b.GetBeginAtomIdx()}", f"{prefix}_{b.GetEndAtomIdx()}", weight=1.0)

    # 2. 共溶媒分子ノードの追加
    mol_v = Chem.AddHs(Chem.MolFromSmiles(cosolv_info["smiles"]))
    for v_i in range(num_cosolv_mols):
        prefix = f"V_{v_i}"
        for a in mol_v.GetAtoms():
            node_id = f"{prefix}_{a.GetIdx()}"
            is_polar = a.GetSymbol() in ['O', 'N', 'S']
            G.add_node(node_id, symbol=a.GetSymbol(), is_polar=is_polar, type='cosolvent')
            if is_polar and (f"{prefix}_center" not in [sc[0] for sc in solvent_centers]):
                solvent_centers.append((node_id, 'cosolvent'))
        for b in mol_v.GetBonds():
            G.add_edge(f"{prefix}_{b.GetBeginAtomIdx()}", f"{prefix}_{b.GetEndAtomIdx()}", weight=1.0)

    # 3. 溶媒間相互作用ネットワーク (Water-Water, Water-Cosolvent, Cosolvent-Cosolvent)
    w_net = 0.8 * f_onsager
    for i in range(len(solvent_centers) - 1):
        sc1, stype1 = solvent_centers[i]
        sc2, stype2 = solvent_centers[i+1]
        w_factor = 1.0 if stype1 == stype2 else 0.85 # 異種溶媒間の混合エントロピー表現
        G.add_edge(sc1, sc2, weight=w_net * w_factor)
        
    if len(solvent_centers) > 2:
        sc1, stype1 = solvent_centers[-1]
        sc2, stype2 = solvent_centers[0]
        G.add_edge(sc1, sc2, weight=w_net * 0.85)

    # 4. 溶質の選択的溶媒和 (Selective Solvation)
    for s_i in range(N_solute):
        s_node = f"S_{s_i}"
        s_polar = G.nodes[s_node]["is_polar"]
        
        for v_center, v_type in solvent_centers:
            if s_polar and v_type == 'water':
                # 極性基 ➔ 水と強力に水素結合
                G.add_edge(s_node, v_center, weight=1.4 * (f_onsager ** 0.5))
            elif not s_polar and v_type == 'cosolvent':
                # 疎水性基 ➔ 共溶媒（エタノールのアルキル基・DMSOメチル基）と選択的親和
                G.add_edge(s_node, v_center, weight=1.1)
            elif s_polar and v_type == 'cosolvent':
                # 極性基 ➔ 共溶媒の極性部分とも結合
                G.add_edge(s_node, v_center, weight=1.0 * (f_onsager ** 0.5))

    return G

def calculate_mixed_solvent_fiedler(solute_smiles, cosolv_key="ethanol", chi_cosolv=0.5):
    """指定の共溶媒モル分率 chi_cosolv における Fiedler λ2 を計算"""
    G = build_binary_mixed_solvent_graph(solute_smiles, cosolv_key=cosolv_key, chi_cosolv=chi_cosolv)
    if G.number_of_nodes() == 0 or not nx.is_connected(G):
        return 0.0
    L_norm = nx.normalized_laplacian_matrix(G).toarray()
    evals = np.linalg.eigvalsh(L_norm)
    return float(evals[1])

def scan_cosolvency_curve(name, smiles, cosolv_key="ethanol"):
    """共溶媒濃度 0% ~ 100% の Fiedler スキャン曲線を生成し、最高溶解度を与える極大点 (Peak) を検出"""
    chis = np.linspace(0.0, 1.0, 11)
    results = []
    
    print(f"\n🧪 【{name}】における {COSOLVENT_DATABASE[cosolv_key]['name']} 共溶媒スキャン")
    print("-" * 75)
    print(f"{'共溶媒モル分率 χ':<18} | {'混合誘電率 ε_mix':<18} | {'Fiedler 溶媒和 λ2'}")
    print("-" * 75)
    
    for chi in chis:
        e_mix = compute_mixed_dielectric(chi, COSOLVENT_DATABASE[cosolv_key]["dielectric"])
        f_val = calculate_mixed_solvent_fiedler(smiles, cosolv_key=cosolv_key, chi_cosolv=chi)
        results.append({
            "chi_cosolv": chi,
            "dielectric_mix": e_mix,
            "fiedler_lambda2": f_val
        })
        print(f"χ = {chi:<14.2f} | ε_mix = {e_mix:<12.1f} | {f_val:.6f}")
        
    print("-" * 75)
    df_res = pd.DataFrame(results)
    best_row = df_res.loc[df_res["fiedler_lambda2"].idxmax()]
    print(f" 🏆 【最高溶解度共溶媒比率 (Cosolvency Peak)】: χ_cosolvent = {best_row['chi_cosolv']:.2f} (最大 λ2 = {best_row['fiedler_lambda2']:.6f})")
    return df_res

if __name__ == "__main__":
    # 例: 水に難溶な安息香酸 (Benzoic Acid) の「水 + エタノール」混合溶媒中での共溶媒極大効果
    scan_cosolvency_curve("安息香酸 (Benzoic Acid)", "O=C(O)c1ccccc1", cosolv_key="ethanol")
    
    # 例: 安息香酸の「水 + DMSO」混合溶媒中での共溶媒極大効果
    scan_cosolvency_curve("安息香酸 (Benzoic Acid)", "O=C(O)c1ccccc1", cosolv_key="dmso")
