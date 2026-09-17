import sys
import os
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

import fiedler_optimization_sim_old as fsim_old
import fiedler_optimization_sim_refined as fsim_refined

# C4H10O isomers (Butanols and Diethyl ether)
isomers = {
    "1-Butanol (n-Butanol)": "CCCCO",
    "Isobutanol": "CC(C)CO",
    "2-Butanol (sec-Butanol)": "CCC(O)C",
    "tert-Butanol": "CC(C)(C)O",
    "Diethyl Ether": "CCOCC"
}

print("=========================================================")
print(" 🧪 Comparison of C4H10O Isomers (Butanols + Ether)")
print("=========================================================")

results_refined = {}
results_old = {}

# Run Refined Simulator
print("\n--- Running Refined Simulator ---")
for name, smiles in isomers.items():
    print(f"Simulating {name}...")
    G_init, w_os = fsim_refined.create_integrated_system([smiles], water_per_atom=3, is_loop=True)
    f_history, _ = fsim_refined.run_fiedler_opt_sim(G_init, w_os, steps=100, allow_long_distance_water=False)
    results_refined[name] = (f_history[0], np.max(f_history), f_history[-1])

# Run Original Simulator
print("\n--- Running Original Simulator ---")
for name, smiles in isomers.items():
    print(f"Simulating {name}...")
    G_init, w_os = fsim_old.create_integrated_system([smiles], water_per_atom=3, is_loop=True)
    f_history, _ = fsim_old.run_fiedler_opt_sim(G_init, w_os, steps=100, allow_long_distance_water=False)
    results_old[name] = (f_history[0], np.max(f_history), f_history[-1])

print("\n=========================================================")
print(" SUMMARY OF RESULTS (C4H10O Isomers)")
print("=========================================================")
print(f"{'Isomer':<36} | {'Refined Peak':<12} | {'Refined Final':<13} | {'Original Peak':<13} | {'Original Final':<14}")
print("-" * 100)
for name in isomers.keys():
    init_r, peak_r, final_r = results_refined[name]
    init_o, peak_o, final_o = results_old[name]
    print(f"{name:<36} | {peak_r:<12.6f} | {final_r:<13.6f} | {peak_o:<13.6f} | {final_o:<14.6f}")
print("=========================================================")
