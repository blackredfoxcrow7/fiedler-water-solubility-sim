import sys
import os
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import fiedler_optimization_sim_old as fsim_old

isomers_c3 = {
    "1-Propanol": "CCCO",
    "2-Propanol": "CC(C)O",
    "Methoxyethane": "COCC"
}

isomers_c4 = {
    "1-Butanol": "CCCCO",
    "Isobutanol": "CC(C)CO",
    "2-Butanol": "CCC(O)C",
    "tert-Butanol": "CC(C)(C)O",
    "Diethyl Ether": "CCOCC"
}

print("=========================================================================")
print(" 🧪 Testing Offset Subtraction on Fiedler Values")
print("=========================================================================")

print("\n--- C3H8O Isomers ---")
print(f"{'Isomer':<15} | {'f_init (Offset)':<15} | {'f_peak':<10} | {'f_final':<10} | {'Δf_final (final - init)':<22} | {'Δf_drop (peak - final)':<22}")
print("-" * 100)

for name, smiles in isomers_c3.items():
    G_init, w_os = fsim_old.create_integrated_system([smiles], water_per_atom=3, is_loop=True)
    f_hist, _ = fsim_old.run_fiedler_opt_sim(G_init, w_os, steps=100)
    
    init = f_hist[0]
    peak = np.max(f_hist)
    final = f_hist[-1]
    
    delta_final = final - init
    delta_drop = peak - final
    
    print(f"{name:<15} | {init:<15.6f} | {peak:<10.6f} | {final:<10.6f} | {delta_final:<22.6f} | {delta_drop:<22.6f}")

print("\n--- C4H10O Isomers ---")
print(f"{'Isomer':<15} | {'f_init (Offset)':<15} | {'f_peak':<10} | {'f_final':<10} | {'Δf_final (final - init)':<22} | {'Δf_drop (peak - final)':<22}")
print("-" * 100)

for name, smiles in isomers_c4.items():
    G_init, w_os = fsim_old.create_integrated_system([smiles], water_per_atom=3, is_loop=True)
    f_hist, _ = fsim_old.run_fiedler_opt_sim(G_init, w_os, steps=100)
    
    init = f_hist[0]
    peak = np.max(f_hist)
    final = f_hist[-1]
    
    delta_final = final - init
    delta_drop = peak - final
    
    print(f"{name:<15} | {init:<15.6f} | {peak:<10.6f} | {final:<10.6f} | {delta_final:<22.6f} | {delta_drop:<22.6f}")

print("=========================================================================")
