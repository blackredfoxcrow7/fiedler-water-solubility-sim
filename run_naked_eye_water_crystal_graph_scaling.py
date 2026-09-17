"""
run_naked_eye_water_crystal_graph_scaling.py
=============================================================================
【肉眼でかすかに見える水結晶 (直径 0.1 mm = 100 μm) のグラフ理論スケーリング計算】

【物理・数学計算パラメータ】
1. 人間の肉眼視認限界: 直径 d = 0.1 mm (100 μm)
2. 水分子 1 個の直径: ~ 0.275 nm (2.75 Å)
3. 氷の密度: 0.917 g/cm3
4. アボガドロ定数: 6.022 x 10^23
=============================================================================
"""

import math

def calculate_naked_eye_water_crystal_graph():
    print("=========================================================================")
    print(" 🔬 肉眼限界水結晶 (直径 0.1 mm = 100 μm) のグラフ理論定量スケーリング")
    print("=========================================================================")
    
    d_naked_eye_mm = 0.1 # 0.1 mm (100 μm)
    r_naked_eye_cm = (d_naked_eye_mm / 2.0) / 10.0 # 0.005 cm
    
    rho_ice = 0.917 # g/cm3
    mw_water = 18.015 # g/mol
    avogadro = 6.02214076e23
    
    d_h2o_nm = 0.275 # nm
    
    vol_cm3 = (4.0 / 3.0) * math.pi * (r_naked_eye_cm ** 3)
    mass_g = vol_cm3 * rho_ice
    moles = mass_g / mw_water
    N_atoms = moles * avogadro
    
    E_bonds = 2.0 * N_atoms
    
    # 直径 100,000 nm / 0.275 nm = 363,636 ステップ
    D_steps = (d_naked_eye_mm * 1e6) / d_h2o_nm
    L_steps = D_steps / 3.0
    
    r_naked_eye_nm = (d_naked_eye_mm / 2.0) * 1e6
    surface_area_nm2 = 4.0 * math.pi * (r_naked_eye_nm ** 2)
    h2o_area_nm2 = math.pi * ((d_h2o_nm / 2.0) ** 2)
    N_surface = surface_area_nm2 / h2o_area_nm2
    
    print("-" * 80)
    print(f" 物理的寸法 (肉眼視認限界)     : 直径 d = {d_naked_eye_mm} mm ({d_naked_eye_mm*1000:.0f} μm)")
    print(f" 結晶の質量                    : m = {mass_g*1e6:.2f} μg (マイクログラム)")
    print("-" * 80)
    print(f" ① ノード数 N (構成水分子数)   : N ≈ {N_atoms:.3e} 個 (約 {N_atoms/1e16:.2f} 京個)")
    print(f" ② エッジ数 |E| (水素結合数)   : |E| ≈ {E_bonds:.3e} 本 (約 {E_bonds/1e16:.2f} 京本)")
    print(f" ③ グラフ直径 D_graph          : D ≈ {D_steps:.2e} ステップ (約 {D_steps/1e4:.1f} 万ステップ)")
    print(f" ④ 平均パス長 L_graph          : L ≈ {L_steps:.2e} ステップ (約 {L_steps/1e4:.1f} 万ステップ)")
    print(f" ⑤ 界面ノード数 |∂S| (チーガー) : |∂S| ≈ {N_surface:.3e} 個 (約 {N_surface/1e8:.1f} 億個)")
    print("-" * 80)

    out_md_path = "/home/eldenring/waterMain/NAKED_EYE_WATER_CRYSTAL_GRAPH_SCALING_REPORT.md"
    with open(out_md_path, "w", encoding="utf-8") as f:
        f.write("# 🔬 肉眼限界水結晶 (直径 0.1 mm = 100 μm) のグラフ理論定量スケーリング報告書\n\n")
        f.write(f"- **物理寸法**: 直径 `d = 0.1 mm = 100 μm`\n")
        f.write(f"- **質量**: `m = {mass_g*1e6:.2f} μg`\n\n")
        f.write(f"1. **ノード数 N (水分子数)**: `N = {N_atoms:.3e}` 個（**約 1.61 京個**）\n")
        f.write(f"2. **エッジ数 |E| (水素結合数)**: `|E| = {E_bonds:.3e}` 本（**約 3.21 京本**）\n")
        f.write(f"3. **グラフ直径 D (ネットワークスパン)**: `D = {D_steps:.2e}` ステップ（**約 36.4 万ステップ**）\n")
        f.write(f"4. **平均パス長 L (平均ステップ長)**: `L = {L_steps:.2e}` ステップ（**約 12.1 万ステップ**）\n")
        f.write(f"5. **界面ノード数 |∂S| (チーガーカット)**: `|∂S| = {N_surface:.3e}` 個（**約 5,289 億個**）\n")
    print(f"\n🎉 肉眼スケールスケーリングレポート作成完了: {out_md_path}")

if __name__ == "__main__":
    calculate_naked_eye_water_crystal_graph()
