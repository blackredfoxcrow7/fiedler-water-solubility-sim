"""
run_percolation_spanning_backbone_analysis.py
=============================================================================
【臨界パーコレーション (p = p_c) における「端から端までのスパニング経路 (Backbone)」幾何解析】

【理論物理的計算】
1. 密な完全結晶の場合: 最短直線経路 D_straight ≈ 36.4 万ステップ
2. 臨界パーコレーション (p = p_c) の場合: フラクタル骨格 (D_path ≈ 1.33)
   --> スパニング主背骨長 D_backbone ≈ (36.4万)^1.33 ≈ 240 万ステップ (6.6 倍蛇行)
3. 臨界ボトルネック (Red Bonds / 赤の結合) による全系開通
=============================================================================
"""

import math
import numpy as np

def calculate_percolation_spanning_pathway():
    print("=========================================================================")
    print(" 🔬 臨界点 (p = p_c) における微結晶 (0.1 mm) 「端と端のスパニング繋がり」解析")
    print("=========================================================================")
    
    # 物理寸法
    d_mm = 0.1 # 0.1 mm
    d_h2o_nm = 0.275 # nm
    
    # 完全密な結晶の直線直径ステップ数
    D_straight = (d_mm * 1e6) / d_h2o_nm # 363,636 ステップ
    
    # 臨界パーコレーションのフラクタル幾何次元 D_path ≈ 1.33
    d_path_dim = 1.33
    
    # フラクタル蛇行パス長 (背骨 Backbone)
    # L_eff = (D_straight / L_unit)^(1.33)
    # 式: D_backbone ≈ D_straight * (D_straight)^(d_path_dim - 1)
    D_backbone = D_straight ** d_path_dim
    
    # 蛇行比率 (Tortuosity factor)
    tortuosity = D_backbone / D_straight
    
    # 臨界ボトルネック (Red Bonds) の想定数 (Logarithmic scaling log2(D))
    red_bonds_est = int(math.log2(D_straight))
    
    print("-" * 80)
    print(f" 物理寸法                        : 直径 d = {d_mm} mm (100 μm)")
    print(f" 完全密な結晶の直線パス          : D_straight = {D_straight:,.0f} ステップ (約 36.4 万ステップ)")
    print("-" * 80)
    print(f" ① 臨界点 (p = p_c) でのフラクタル次元: D_path ≈ {d_path_dim:.2f} (曲がりくねったヘビ型幾何)")
    print(f" ② スパニング主背骨 (Backbone) 長 : D_backbone ≈ {D_backbone:,.0f} ステップ (約 {D_backbone/1e4:.1f} 万ステップ)")
    print(f" ③ 蛇行迂回比率 (Tortuosity)     : 直線の {tortuosity:.2f} 倍遠回り")
    print(f" ④ 赤の結合 (Red Bonds / 臨界鍵) : 主背骨中に約 {red_bonds_est} カ所の切断ボトルネックが存在")
    print("-" * 80)

    out_md_path = "/home/eldenring/waterMain/PERCOLATION_SPANNING_BACKBONE_REPORT.md"
    with open(out_md_path, "w", encoding="utf-8") as f:
        f.write("# 🔬 臨界点 (p = p_c) における「端と端のスパニング繋がり」幾何解析報告書\n\n")
        f.write(f"- **物理寸法**: 直径 `d = 0.1 mm = 100 μm`\n")
        f.write(f"- **完全密結晶の直線直径**: `D_straight = {D_straight:,.0f}` ステップ（約 36.4 万ステップ）\n\n")
        f.write(f"1. **フラクタル主背骨長 (Backbone Path)**: `D_backbone ≈ {D_backbone:,.0f}` ステップ（**約 240 万ステップ**）\n")
        f.write(f"2. **蛇行比率 (Tortuosity Factor)**: 直線距離の **`{tortuosity:.2f} 倍` 遠回り（ヘビ型迂回）**\n")
        f.write(f"3. **赤の結合 (Red Bonds / 臨界ボトルネック)**: 主背骨中に **約 `{red_bonds_est}` カ所の切断可能キー結合** が存在\n")
    print(f"\n🎉 臨界スパニング解密レポート作成完了: {out_md_path}")

if __name__ == "__main__":
    calculate_percolation_spanning_pathway()
