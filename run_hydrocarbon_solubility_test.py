import sys
import os

sys.path.append("/home/eldenring/waterMain")
import fiedler_sim_hbond_offset_clash_free as fsim

# 代表的な化学物質 (炭化水素、アルコール、エーテル、ケトン、芳香族)
test_molecules = [
    # 炭化水素類 (Hydrocarbons)
    ("メタン", "C", "難溶 (0.0022 g/100g)"),
    ("エタン", "CC", "難溶 (0.006 g/100g)"),
    ("プロパン", "CCC", "難溶 (0.007 g/100g)"),
    ("n-ブタン", "CCCC", "難溶 (0.006 g/100g)"),
    ("n-ヘキサン", "CCCCCC", "極めて難溶 (0.00095 g/100g)"),
    ("ベンゼン", "c1ccccc1", "難溶 (0.18 g/100g)"),
    ("トルエン", "Cc1ccccc1", "難溶 (0.05 g/100g)"),
    
    # アルコール・ポリオール類 (Alcohols / Polyols)
    ("メタノール", "CO", "完全混和"),
    ("エタノール", "CCO", "完全混和"),
    ("1-プロパノール", "CCCO", "完全混和"),
    ("1-ブタノール", "CCCCO", "可溶 (7.3 g/100g)"),
    ("グリセリン", "OCC(O)CO", "完全混和"),
    
    # エーテル・ケトン類 (Ethers / Ketones)
    ("アセトン", "CC(=O)C", "完全混和"),
    ("THF (テトラヒドロフラン)", "C1CCOC1", "完全混和"),
    ("ジエチルエーテル", "CCOCC", "可溶 (6.9 g/100g)")
]

print("=========================================================================")
print(" 🧪 代表的化学物質（炭化水素含む）の水溶性シミュレーション検証")
print("=========================================================================")
print(f"{'物質名':<16} | {'SMILES':<12} | {'実効変化量 Δf_final':<18} | {'実際の水溶性データ'}")
print("-" * 75)

results = []
for name, smiles, exp_sol in test_molecules:
    G_init, w_os = fsim.create_integrated_system([smiles], water_per_atom=3, is_loop=True)
    f_history, frames_data, delta_f_final, _ = fsim.run_fiedler_opt_sim(G_init, w_os, steps=100)
    results.append((name, smiles, delta_f_final, exp_sol))
    print(f"{name:<16} | {smiles:<12} | {delta_f_final:<18.6f} | {exp_sol}")

print("=========================================================================")
