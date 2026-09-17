"""
example_multi_solvent_usage.py
=============================================================================
「8 大有機溶媒物性エンジン (multi_solvent_physical_properties_engine.py)」の使い方デモ
=============================================================================
"""

import sys
import os

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import multi_solvent_physical_properties_engine as mpe

def main():
    print("=========================================================================")
    print(" 🚀 8 大有機溶媒物性エンジン 使い方デモ（カフェイン & イブプロフェン）")
    print("=========================================================================")
    
    # 例 1: カフェイン (SMILES: CN1C=NC2=C1C(=O)N(C(=O)N2C)C) の全溶媒比較
    mpe.run_multi_solvent_benchmark("カフェイン", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
    
    # 例 2: イブプロフェン (SMILES: CC(C)Cc1ccc(cc1)C(C)C(=O)O) の特定溶媒計算
    smiles_ibuprofen = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
    
    f_water = mpe.calculate_multi_solvent_fiedler(smiles_ibuprofen, solvent_key="water")
    f_dmso  = mpe.calculate_multi_solvent_fiedler(smiles_ibuprofen, solvent_key="dmso")
    f_hex   = mpe.calculate_multi_solvent_fiedler(smiles_ibuprofen, solvent_key="hexane")
    
    print("\n💊 【イブプロフェン個別溶媒比較】")
    print(f"  ・水 (Water) 中での Fiedler 値     : {f_water:.6f} (水には難溶)")
    print(f"  ・DMSO 中での Fiedler 値          : {f_dmso:.6f}  (DMSO に高溶解)")
    print(f"  ・n-ヘキサン (Hexane) 中の Fiedler値: {f_hex:.6f}  (非極性溶媒和)")

if __name__ == "__main__":
    main()
