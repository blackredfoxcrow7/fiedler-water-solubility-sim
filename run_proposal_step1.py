import sys
import os
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import fiedler_sim_hbond_offset_clash_free as fsim
import visualize_toggle_water as vis

# 提案書 Step 1: 多価アルコールおよび環状 vs 直鎖分子の検証

proposals_molecules = {
    # 1. 多価アルコール (ポリオール)
    "Ethylene Glycol (OCCO)": ("OCCO", "Miscible (完全混和)", "Polyol C2"),
    "1,4-Butanediol (OCCCCO)": ("OCCCCO", "Miscible (完全混和)", "Polyol C4"),
    "Glycerin (OCC(O)CO)": ("OCC(O)CO", "Miscible (完全混和)", "Polyol C3"),
    
    # 2. 環状 vs 直鎖 比較ペア
    "THF (Tetrahydrofuran)": ("C1CCCO1", "Miscible (完全混和)", "Cyclic Ether C4"),
    "Diethyl Ether": ("CCOCC", "6.9 g / 100mL", "Linear Ether C4"),
    "Cyclohexanol": ("C1CCC(O)CC1", "3.6 g / 100mL", "Cyclic Alcohol C6"),
    "1-Hexanol": ("CCCCCCO", "0.59 g / 100mL", "Linear Alcohol C6"),
}

print("====================================================================================================")
print(" 🧪 Proposal Step 1: Polyols & Cyclic vs Linear Topological Solvation Simulation")
print("====================================================================================================")
print(f"{'Molecule Name':<28} | {'f_init (Offset)':<15} | {'f_peak':<10} | {'f_final':<10} | {'Δf_final':<12} | {'Physical Sol.':<20}")
print("-" * 105)

results = {}

for name, (smiles, sol, cat) in proposals_molecules.items():
    G_init, w_os = fsim.create_integrated_system([smiles], water_per_atom=3, is_loop=True)
    f_hist, frames_data, delta_f_final, _ = fsim.run_fiedler_opt_sim(G_init, w_os, steps=100)
    
    init = f_hist[0]
    peak = np.max(f_hist)
    final = f_hist[-1]
    
    results[name] = {
        'init': init, 'peak': peak, 'final': final, 'delta': delta_f_final, 'sol': sol, 'cat': cat
    }
    
    print(f"{name:<28} | {init:<15.6f} | {peak:<10.6f} | {final:<10.6f} | {delta_f_final:<12.6f} | {sol:<20}")

print("====================================================================================================")
