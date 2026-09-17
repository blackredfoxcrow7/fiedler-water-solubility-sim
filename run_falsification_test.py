import sys
import os

sys.path.append("/home/eldenring/waterMain")
import fiedler_sim_hbond_offset_clash_free as fsim

# 科学的反検証 (Falsification Test) 用のエッジケース・難問化合物ペア
# カール・ポパーの反証可能性原理に基づき、モデルの限界やズレが生じる可能性のある領域をテストする

falsification_pairs = [
    # 1. 分子内水素結合 (Intramolecular H-Bonding Isomers)
    {
        "category": "分子内水素結合 (Intramolecular H-Bond)",
        "mol1": ("サリチル酸 (2-ヒドロキシ安息香酸)", "O=C(O)c1ccccc1O", "分子内H結合形成 ➔ 難溶 (2.24 g/L)"),
        "mol2": ("4-ヒドロキシ安息香酸", "O=C(O)c1ccc(O)cc1", "分子外水和可能 ➔ 可溶 (5.0 g/L)"),
        "note": "分子内H結合により水と結合する極性基が自己閉塞する現象"
    },
    # 2. Cis/Trans 幾何異性体 ＋ 結晶格子エネルギー (Cis/Trans & Lattice Energy)
    {
        "category": "Cis/Trans 幾何異性体 (Cis/Trans Isomers)",
        "mol1": ("マレイン酸 (cis型)", "O=C(O)/C=C\\C(=O)O", "水和エネルギー高 ➔ 非常に高可溶 (788 g/L)"),
        "mol2": ("フマル酸 (trans型)", "O=C(O)/C=C/C(=O)O", "結晶格子強固 ➔ 難溶 (6.3 g/L)"),
        "note": "結晶格子エネルギー（融点差: マレイン酸139℃ vs フマル酸287℃）の影響"
    },
    # 3. 芳香族位置異性体 (Aromatic Positional Isomers)
    {
        "category": "芳香族位置異性体 (Positional Isomers)",
        "mol1": ("フタル酸 (1,2-置換)", "O=C(O)c1ccccc1C(=O)O", "可溶 (16.2 g/L)"),
        "mol2": ("テレフタル酸 (1,4-置換)", "O=C(O)c1ccc(C(=O)O)cc1", "極めて難溶 (0.015 g/L)"),
        "note": "テレフタル酸の異常な高い結晶対称性と強い格子エネルギー"
    },
    # 4. イオン化・極性塩 (Ionization / Zwitterions)
    {
        "category": "イオン化・電荷分極 (Charge Polarization)",
        "mol1": ("酢酸 (分子型)", "CC(=O)O", "完全混和"),
        "mol2": ("酢酸ナトリウム (イオン型)", "CC(=O)[O-].[Na+]", "極めて高可溶 (1250 g/L)"),
        "note": "静電イオン和効果 (Coulombic Solvation)"
    }
]

print("=========================================================================")
print(" 🧐 科学的反検証テスト (Falsification & Edge-Case Benchmark)")
print("=========================================================================")

for test in falsification_pairs:
    print(f"\n📌 カテゴリ: {test['category']}")
    print(f"   物理的注目点: {test['note']}")
    print("-" * 75)
    
    name1, smiles1, exp1 = test["mol1"]
    name2, smiles2, exp2 = test["mol2"]
    
    # Sim 1
    G1, w1 = fsim.create_integrated_system([smiles1], water_per_atom=3, is_loop=True)
    _, _, df1, _ = fsim.run_fiedler_opt_sim(G1, w1, steps=100)
    
    # Sim 2
    G2, w2 = fsim.create_integrated_system([smiles2], water_per_atom=3, is_loop=True)
    _, _, df2, _ = fsim.run_fiedler_opt_sim(G2, w2, steps=100)
    
    print(f"  分子1: {name1:<25} | Δf = {df1:.6f} | 実験: {exp1}")
    print(f"  分子2: {name2:<25} | Δf = {df2:.6f} | 実験: {exp2}")
    
    delta_diff = df1 - df2
    print(f"  ==> シミュレーション評価値の差 (Δf1 - Δf2) = {delta_diff:+.6f}")

print("\n=========================================================================")
