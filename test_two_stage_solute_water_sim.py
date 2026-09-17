import sys
import os

sys.path.append("/home/eldenring/waterMain")
import fiedler_sim_hbond_offset_clash_free as fsim
import networkx as nx
import numpy as np

# 二段階グラフシミュレーション (Two-Stage Graph Physics Simulator)
# Stage 1: 溶質のみの純粋相・凝集相ネットワークの最適化 (f_solute)
# Stage 2: 水分子を加えた水和・溶解相ネットワークの最適化 (f_solvated)
# 物理量: 二段階実効溶解度 ΔΔf = f_solvated - f_solute (水和安定化 - 結晶/自己凝集ペナルティ)

def run_two_stage_simulation(smiles, num_solutes=2, water_per_atom=3, steps=100):
    # --- Stage 1: 溶質のみの純粋相ネットワーク (Solute Aggregate Graph) ---
    G_solute, w_os_solute = fsim.create_integrated_system([smiles] * num_solutes, water_per_atom=0, is_loop=False)
    f_hist_solute, _, df_solute, _ = fsim.run_fiedler_opt_sim(G_solute, w_os_solute, steps=steps)
    f_final_solute = f_hist_solute[-1]
    
    # --- Stage 2: 水分子を加えた水和相ネットワーク (Solvated System Graph) ---
    G_solvated, w_os_solvated = fsim.create_integrated_system([smiles] * num_solutes, water_per_atom=water_per_atom, is_loop=True)
    f_hist_solvated, _, df_solvated, _ = fsim.run_fiedler_opt_sim(G_solvated, w_os_solvated, steps=steps)
    f_final_solvated = f_hist_solvated[-1]
    
    # 真の二段階溶解度 ΔΔf (Thermodynamic Solvation Delta)
    delta_delta_f = df_solvated - df_solute
    
    return f_final_solute, f_final_solvated, df_solute, df_solvated, delta_delta_f

print("=========================================================================")
print(" 🚀 2段階グラフシミュレータ (Stage 1: 溶質自己凝集 ➔ Stage 2: 水和溶解)")
print("=========================================================================")

pairs = [
    ("芳香族位置異性体", "フタル酸 (1,2-置換 / 可溶 16.2g/L)", "O=C(O)c1ccccc1C(=O)O", 
                        "テレフタル酸 (1,4-置換 / 超難溶 0.015g/L)", "O=C(O)c1ccc(C(=O)O)cc1"),
    ("Cis/Trans異性体",  "マレイン酸 (cis / 超高可溶 788g/L)", "O=C(O)/C=C\\C(=O)O", 
                        "フマル酸 (trans / 難溶 6.3g/L)", "O=C(O)/C=C/C(=O)O")
]

for cat, name1, smiles1, name2, smiles2 in pairs:
    print(f"\n📌 カテゴリ: {cat}")
    print("-" * 75)
    
    s1_sol, s1_solv, df1_sol, df1_solv, ddf1 = run_two_stage_simulation(smiles1, num_solutes=2)
    s2_sol, s2_solv, df2_sol, df2_solv, ddf2 = run_two_stage_simulation(smiles2, num_solutes=2)
    
    print(f"  分子1: {name1:<35}")
    print(f"        Stage 1 (凝集相 Δf_solute) = {df1_sol:.6f} | Stage 2 (水和相 Δf_solvated) = {df1_solv:.6f}")
    print(f"        ==> 二段階真の溶解度 ΔΔf = {ddf1:.6f}")
    
    print(f"  分子2: {name2:<35}")
    print(f"        Stage 1 (凝集相 Δf_solute) = {df2_sol:.6f} | Stage 2 (水和相 Δf_solvated) = {df2_solv:.6f}")
    print(f"        ==> 二段階真の溶解度 ΔΔf = {ddf2:.6f}")
    
    diff_ddf = ddf1 - ddf2
    print(f"  ==> 2段階評価の相対差 (ΔΔf1 - ΔΔf2) = {diff_ddf:+.6f}")

print("\n=========================================================================")
