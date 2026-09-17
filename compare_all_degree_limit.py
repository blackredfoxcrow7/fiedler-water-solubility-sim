import sys
import os
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

import fiedler_optimization_sim_degree_limit as fsim_deg
import fiedler_optimization_sim_refined as fsim_refined

# Define all C3 and C4 isomers
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
print(" 🧪 Comparison with Solute Degree Limit implemented in Old Simulator")
print("=========================================================================")

# 1. C3H8O Isomers
print("\n--- C3H8O Isomers ---")
print(f"{'Isomer':<20} | {'Refined Peak':<12} | {'Refined Final':<13} | {'Deg-Limit Peak':<14} | {'Deg-Limit Final':<15}")
print("-" * 85)
for name, smiles in isomers_c3.items():
    # Refined
    G_init_r, w_os_r = fsim_refined.create_integrated_system([smiles], water_per_atom=3, is_loop=True)
    f_hist_r, _ = fsim_refined.run_fiedler_opt_sim(G_init_r, w_os_r, steps=100)
    # Degree-Limited Old
    G_init_d, w_os_d = fsim_deg.create_integrated_system([smiles], water_per_atom=3, is_loop=True)
    f_hist_d, _ = fsim_deg.run_fiedler_opt_sim(G_init_d, w_os_d, steps=100)
    
    print(f"{name:<20} | {np.max(f_hist_r):<12.6f} | {f_hist_r[-1]:<13.6f} | {np.max(f_hist_d):<14.6f} | {f_hist_d[-1]:<15.6f}")

# 2. C4H10O Isomers
print("\n--- C4H10O Isomers ---")
print(f"{'Isomer':<20} | {'Refined Peak':<12} | {'Refined Final':<13} | {'Deg-Limit Peak':<14} | {'Deg-Limit Final':<15}")
print("-" * 85)
for name, smiles in isomers_c4.items():
    # Refined
    G_init_r, w_os_r = fsim_refined.create_integrated_system([smiles], water_per_atom=3, is_loop=True)
    f_hist_r, _ = fsim_refined.run_fiedler_opt_sim(G_init_r, w_os_r, steps=100)
    # Degree-Limited Old
    G_init_d, w_os_d = fsim_deg.create_integrated_system([smiles], water_per_atom=3, is_loop=True)
    f_hist_d, _ = fsim_deg.run_fiedler_opt_sim(G_init_d, w_os_d, steps=100)
    
    print(f"{name:<20} | {np.max(f_hist_r):<12.6f} | {f_hist_r[-1]:<13.6f} | {np.max(f_hist_d):<14.6f} | {f_hist_d[-1]:<15.6f}")

print("=========================================================================")
