import sys
import os

sys.path.append("/home/eldenring/waterMain")
import fiedler_sim_hbond_offset_clash_free as fsim

# 2分子（ダイマー / 多分子）同時投入による反証データ（結晶化・分子間スタッキング）の克服テスト

def run_dimer_falsification_test():
    print("=========================================================================")
    print(" 🧪 2分子（ダイマー）システムによる反証データの克服・検証シミュレーション")
    print("=========================================================================")
    
    # テスト 1: 芳香族位置異性体 (フタル酸 1,2-置換 vs テレフタル酸 1,4-置換) 2分子投入
    smiles_phthalic = "O=C(O)c1ccccc1C(=O)O"
    smiles_terephthalic = "O=C(O)c1ccc(C(=O)O)cc1"
    
    print("\n📌 実験 1: 芳香族位置異性体（2分子システムでの分子間スタッキング・自己会合）")
    print("   物理機構: テレフタル酸は2分子が対角線上に向き合い強固な結晶格子を作る")
    print("-" * 75)
    
    # 単一分子 vs 2分子
    G_ph_1, w_ph_1 = fsim.create_integrated_system([smiles_phthalic], water_per_atom=3, is_loop=True)
    _, _, df_ph_1, _ = fsim.run_fiedler_opt_sim(G_ph_1, w_ph_1, steps=100)
    
    G_ph_2, w_ph_2 = fsim.create_integrated_system([smiles_phthalic, smiles_phthalic], water_per_atom=3, is_loop=True)
    _, _, df_ph_2, _ = fsim.run_fiedler_opt_sim(G_ph_2, w_ph_2, steps=100)
    
    G_tp_1, w_tp_1 = fsim.create_integrated_system([smiles_terephthalic], water_per_atom=3, is_loop=True)
    _, _, df_tp_1, _ = fsim.run_fiedler_opt_sim(G_tp_1, w_tp_1, steps=100)
    
    G_tp_2, w_tp_2 = fsim.create_integrated_system([smiles_terephthalic, smiles_terephthalic], water_per_atom=3, is_loop=True)
    _, _, df_tp_2, _ = fsim.run_fiedler_opt_sim(G_tp_2, w_tp_2, steps=100)
    
    print(f"  フタル酸   (1分子) Δf = {df_ph_1:.6f} | (2分子ダイマー) Δf = {df_ph_2:.6f}")
    print(f"  テレフタル酸 (1分子) Δf = {df_tp_1:.6f} | (2分子ダイマー) Δf = {df_tp_2:.6f}")
    
    # 2分子化による水和抑制の強さ
    drop_ph = df_ph_2 - df_ph_1
    drop_tp = df_tp_2 - df_tp_1
    print(f"  ==> フタル酸の2分子化による変化   : {drop_ph:+.6f}")
    print(f"  ==> テレフタル酸の2分子化による変化 : {drop_tp:+.6f}")
    
    print("\n=========================================================================")

if __name__ == "__main__":
    run_dimer_falsification_test()
