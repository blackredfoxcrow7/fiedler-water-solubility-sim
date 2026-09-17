import sys
import os
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import fiedler_sim_hbond_offset as fsim

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

print("====================================================================================================")
print(" 🧪 Hydrogen-Conscious Valency & Offset Subtraction Simulation Results (fiedler_sim_hbond_offset)")
print("====================================================================================================")

print("\n--- C3H8O Isomers ---")
print(f"{'Isomer':<18} | {'f_init (Offset)':<15} | {'f_peak':<10} | {'f_final':<10} | {'Δf_final (final - init)':<22} | {'Physical Sol.':<15}")
print("-" * 100)

for name, smiles in isomers_c3.items():
    G_init, w_os = fsim.create_integrated_system([smiles], water_per_atom=3, is_loop=True)
    f_hist, _, delta_f_final, _ = fsim.run_fiedler_opt_sim(G_init, w_os, steps=100)
    
    init = f_hist[0]
    peak = np.max(f_hist)
    final = f_hist[-1]
    
    sol = "Miscible" if "Propanol" in name else "8.2 g/100mL"
    
    print(f"{name:<18} | {init:<15.6f} | {peak:<10.6f} | {final:<10.6f} | {delta_f_final:<22.6f} | {sol:<15}")

print("\n--- C4H10O Isomers ---")
print(f"{'Isomer':<18} | {'f_init (Offset)':<15} | {'f_peak':<10} | {'f_final':<10} | {'Δf_final (final - init)':<22} | {'Physical Sol.':<15}")
print("-" * 100)

sol_map = {
    "tert-Butanol": "Miscible",
    "2-Butanol": "29.0 g/100mL",
    "Isobutanol": "8.5 g/100mL",
    "1-Butanol": "7.3 g/100mL",
    "Diethyl Ether": "6.9 g/100mL"
}

for name, smiles in isomers_c4.items():
    G_init, w_os = fsim.create_integrated_system([smiles], water_per_atom=3, is_loop=True)
    f_hist, _, delta_f_final, _ = fsim.run_fiedler_opt_sim(G_init, w_os, steps=100)
    
    init = f_hist[0]
    peak = np.max(f_hist)
    final = f_hist[-1]
    
    print(f"{name:<18} | {init:<15.6f} | {peak:<10.6f} | {final:<10.6f} | {delta_f_final:<22.6f} | {sol_map[name]:<15}")

print("====================================================================================================")
