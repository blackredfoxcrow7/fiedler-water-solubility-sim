import sys
import os

sys.path.append("/home/eldenring/waterMain")
import fiedler_sim_hbond_offset_clash_free as fsim

# 代表的なアミノ酸 (20種アミノ酸の主要分類グループから抽出)
amino_acids = [
    # 1. 脂肪族・非極性 (Aliphatic / Non-polar)
    ("グリシン (Gly)", "NCC(=O)O", "非極性・最少", "高可溶 (24.9 g/100g)"),
    ("アラニン (Ala)", "CC(N)C(=O)O", "非極性", "可溶 (16.7 g/100g)"),
    ("バリン (Val)", "CC(C)C(N)C(=O)O", "非極性・疎水性", "やや難溶 (8.85 g/100g)"),
    ("ロイシン (Leu)", "CC(C)CC(N)C(=O)O", "非極性・長鎖疎水", "難溶 (2.4 g/100g)"),
    ("プロリン (Pro)", "C1CC(NC1)C(=O)O", "環状アミノ酸", "極めて高可溶 (162 g/100g)"),
    
    # 2. 芳香族 (Aromatic)
    ("フェニルアラニン (Phe)", "NC(Cc1ccccc1)C(=O)O", "芳香族・強疎水性", "難溶 (2.96 g/100g)"),
    ("チロシン (Tyr)", "NC(Cc1ccc(O)cc1)C(=O)O", "芳香族・フェノールOH", "極めて難溶 (0.045 g/100g)"),
    ("トリプトファン (Trp)", "NC(Cc1c[nH]c2ccccc12)C(=O)O", "芳香族・インドール", "難溶 (1.14 g/100g)"),
    
    # 3. 極性・電荷なし (Polar Uncharged)
    ("セリン (Ser)", "OCC(N)C(=O)O", "極性・水酸基", "高可溶 (50.3 g/100g)"),
    ("トレオニン (Thr)", "CC(O)C(N)C(=O)O", "極性・水酸基", "高可溶 (20.5 g/100g)"),
    
    # 4. 親水性・酸性/塩基性 (Charged Acidic / Basic)
    ("グルタミン酸 (Glu)", "NC(CCC(=O)O)C(=O)O", "酸性・ジカルボン酸", "可溶 (0.86 g/100g)"),
    ("リシン (Lys)", "NCCCC(N)C(=O)O", "塩基性・アミノ基", "極めて高可溶 (>100 g/100g)")
]

print("=========================================================================")
print(" 🧬 代表的アミノ酸の水溶性・水和シミュレーション検証")
print("=========================================================================")
print(f"{'アミノ酸名':<22} | {'SMILES':<20} | {'実効変化量 Δf_final':<18} | {'実際の水溶性データ'}")
print("-" * 85)

for name, smiles, group, exp_sol in amino_acids:
    G_init, w_os = fsim.create_integrated_system([smiles], water_per_atom=3, is_loop=False)
    f_history, frames_data, delta_f_final, _ = fsim.run_fiedler_opt_sim(G_init, w_os, steps=100)
    print(f"{name:<22} | {smiles:<20} | {delta_f_final:<18.6f} | {exp_sol}")

print("=========================================================================")
