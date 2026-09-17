import sys
import os
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

import fiedler_optimization_sim_old as fsim

# Propanol isomers and constitutional isomers (C3H8O)
isomers = {
    "1-Propanol (Propan-1-ol)": "CCCO",
    "2-Propanol (Isopropanol)": "CC(C)O",
    "Methoxyethane (Methyl Ethyl Ether)": "COCC"
}

print("=========================================================")
# Print the header of the comparison run (Original simulator)
print(" 🧪 Comparison of C3H8O Isomers using Original Fiedler Simulator")
print("=========================================================")

for name, smiles in isomers.items():
    print(f"\nRunning simulation for {name} (SMILES: {smiles})...")
    # Initialize system with water_per_atom=3
    G_init, w_os = fsim.create_integrated_system([smiles], water_per_atom=3, is_loop=True)
    f_history, frames_data = fsim.run_fiedler_opt_sim(G_init, w_os, steps=100, allow_long_distance_water=False)
    
    init_f = f_history[0]
    peak_f = np.max(f_history)
    final_f = f_history[-1]
    
    print(f"  Initial Fiedler Value : {init_f:.6f}")
    print(f"  Peak Fiedler Value    : {peak_f:.6f}")
    print(f"  Final Fiedler Value   : {final_f:.6f}")

print("\n=========================================================")
