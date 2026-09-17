import sys
import os
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

import fiedler_optimization_sim_old as fsim_old
import fiedler_optimization_sim_trajectory_corrected as fsim_traj

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
print(" 🧪 Verification of Purely Trajectory-Based Fiedler Correction")
print("=========================================================================")

print("\n--- C3H8O Isomers ---")
print(f"{'Isomer':<20} | {'Original Peak':<13} | {'Original Final':<14} | {'Traj-Corrected Score':<20}")
print("-" * 75)
for name, smiles in isomers_c3.items():
    # Original
    G_init_o, w_os_o = fsim_old.create_integrated_system([smiles], water_per_atom=3, is_loop=True)
    f_hist_o, _ = fsim_old.run_fiedler_opt_sim(G_init_o, w_os_o, steps=100)
    # Trajectory-Corrected
    G_init_t, w_os_t = fsim_traj.create_integrated_system([smiles], water_per_atom=3, is_loop=True)
    f_hist_t, _, score = fsim_traj.run_fiedler_opt_sim(G_init_t, w_os_t, steps=100)
    
    print(f"{name:<20} | {np.max(f_hist_o):<13.6f} | {f_hist_o[-1]:<14.6f} | {score:<20.6f}")

print("\n--- C4H10O Isomers ---")
print(f"{'Isomer':<20} | {'Original Peak':<13} | {'Original Final':<14} | {'Traj-Corrected Score':<20}")
print("-" * 75)
for name, smiles in isomers_c4.items():
    # Original
    G_init_o, w_os_o = fsim_old.create_integrated_system([smiles], water_per_atom=3, is_loop=True)
    f_hist_o, _ = fsim_old.run_fiedler_opt_sim(G_init_o, w_os_o, steps=100)
    # Trajectory-Corrected
    G_init_t, w_os_t = fsim_traj.create_integrated_system([smiles], water_per_atom=3, is_loop=True)
    f_hist_t, _, score = fsim_traj.run_fiedler_opt_sim(G_init_t, w_os_t, steps=100)
    
    print(f"{name:<20} | {np.max(f_hist_o):<13.6f} | {f_hist_o[-1]:<14.6f} | {score:<20.6f}")

print("=========================================================================")
