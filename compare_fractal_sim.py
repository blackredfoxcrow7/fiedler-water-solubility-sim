import sys
import os
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import fiedler_sim_hbond_offset_clash_free as fsim_standard
import fiedler_fractal_water_sim as fsim_fractal

isomers = {
    "tert-Butanol": "CC(C)(C)O",
    "Diethyl Ether": "CCOCC",
    "Glycerin": "OCC(O)CO"
}

print("==========================================================================================")
print(" 🧪 Standard vs Active Fractal Water Network Simulation Comparison")
print("==========================================================================================")
print(f"{'Molecule':<18} | {'Standard Δf_final':<20} | {'Active-Fractal Δf_final':<25} | {'Enhancement':<15}")
print("-" * 85)

for name, smiles in isomers.items():
    # Standard
    G_init1, w_os1 = fsim_standard.create_integrated_system([smiles], water_per_atom=3, is_loop=True)
    _, _, delta_std, _ = fsim_standard.run_fiedler_opt_sim(G_init1, w_os1, steps=100)
    
    # Active Fractal
    G_init2, w_os2 = fsim_fractal.fsim_base.create_integrated_system([smiles], water_per_atom=3, is_loop=True)
    _, _, delta_frac = fsim_fractal.run_fiedler_fractal_opt_sim(G_init2, w_os2, steps=100, alpha_powerlow=2.0)
    
    diff = delta_frac - delta_std
    print(f"{name:<18} | {delta_std:<20.6f} | {delta_frac:<25.6f} | {diff:<+15.6f}")

print("==========================================================================================")
