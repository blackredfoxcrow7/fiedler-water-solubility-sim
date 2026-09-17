"""
fiedler_two_stage_simulator.py
=============================================================================
二段階グラフ熱力学シミュレータ (Two-Stage Graph Physics Simulator)

【物理理論】
溶解度 ΔΔf = f_solvated (Stage 2: 水和相) - f_solute (Stage 1: 凝集・結晶相)
・Stage 1: 溶質分子同士の自己会合・結晶格子エネルギー (f_solute)
・Stage 2: 水分子を加えた統合水和ネットワーク (f_solvated)
・二段階実効溶解度 ΔΔf: 結晶解体ペナルティを差し引いた真の熱力学的溶解度尺度
=============================================================================
"""

import sys
import os
import networkx as nx
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import fiedler_sim_hbond_offset_clash_free as fsim

def run_two_stage_solubility_sim(smiles, num_solutes=2, water_per_atom=3, steps=100):
    """
    二段階グラフシミュレーションを実行し、
    Stage 1 (凝集相) と Stage 2 (水和相) の差分 ΔΔf を算出する。
    """
    # --- Stage 1: 溶質純粋相 (自己凝集・結晶格子相) ---
    G_solute, w_os_solute = fsim.create_integrated_system([smiles] * num_solutes, water_per_atom=0, is_loop=False)
    f_hist_solute, _, df_solute, _ = fsim.run_fiedler_opt_sim(G_solute, w_os_solute, steps=steps)
    f_final_solute = f_hist_solute[-1]
    
    # --- Stage 2: 水和溶解相 (水和統合ネットワーク相) ---
    G_solvated, w_os_solvated = fsim.create_integrated_system([smiles] * num_solutes, water_per_atom=water_per_atom, is_loop=True)
    f_hist_solvated, _, df_solvated, _ = fsim.run_fiedler_opt_sim(G_solvated, w_os_solvated, steps=steps)
    f_final_solvated = f_hist_solvated[-1]
    
    # 二段階真の実効溶解度 ΔΔf
    delta_delta_f = df_solvated - df_solute
    
    return {
        "smiles": smiles,
        "f_solute": f_final_solute,
        "f_solvated": f_final_solvated,
        "df_solute": df_solute,
        "df_solvated": df_solvated,
        "delta_delta_f": delta_delta_f
    }

if __name__ == "__main__":
    res = run_two_stage_solubility_sim("O=C(O)c1ccccc1C(=O)O", num_solutes=2)
    print(f"フタル酸 二段階溶解度 ΔΔf: {res['delta_delta_f']:.6f}")
