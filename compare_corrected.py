import sys
import os
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

import fiedler_optimization_sim_degree_limit as fsim

isomers_c3 = {
    "1-Propanol (CCCO)": "CCCO",
    "2-Propanol (CC(C)O)": "CC(C)O",
    "Methoxyethane (COCC)": "COCC"
}

isomers_c4 = {
    "1-Butanol (CCCCO)": "CCCCO",
    "Isobutanol (CC(C)CO)": "CC(C)CO",
    "2-Butanol (CCC(O)C)": "CCC(O)C",
    "tert-Butanol (CC(C)(C)O)": "CC(C)(C)O",
    "Diethyl Ether (CCOCC)": "CCOCC"
}

def run_and_calc_dsi(smiles):
    G_init, w_os = fsim.create_integrated_system([smiles], water_per_atom=3, is_loop=True)
    f_history, _ = fsim.run_fiedler_opt_sim(G_init, w_os, steps=100)
    
    f_history = np.array(f_history)
    peak = np.max(f_history)
    final = f_history[-1]
    
    # Calculate different correction formulas:
    # 1. DSI_RelativeDrop: final * (1 - final/peak)
    dsi_rel = final * (1.0 - final / peak) if peak > 0 else 0
    
    # 2. DSI_AbsDrop: final * (peak - final)
    dsi_abs = final * (peak - final)
    
    # 3. DSI_Fluctuation: final * std(f_history[10:])
    dsi_fluc = final * np.std(f_history[10:]) if len(f_history) > 10 else 0
    
    # 4. DSI_DynamicRatio: final / (peak / final) = final^2 / peak
    dsi_ratio = (final ** 2) / peak if peak > 0 else 0
    
    # 5. DSI_LogDrop: final * log(1 + peak - final)
    dsi_log = final * np.log(1.0 + peak - final)

    # 6. DSI_EntropyPenalized: final * (1 - final/peak)^2
    dsi_ent = final * ((1.0 - final / peak) ** 2) if peak > 0 else 0

    return peak, final, dsi_rel, dsi_abs, dsi_fluc, dsi_ratio, dsi_log, dsi_ent

print("==========================================================================================")
print(" Testing Trajectory-based Correction Formulas (Old Simulator w/ Solute Degree Limit)")
print("==========================================================================================")

for name, smiles in {**isomers_c3, **isomers_c4}.items():
    peak, final, dsi_rel, dsi_abs, dsi_fluc, dsi_ratio, dsi_log, dsi_ent = run_and_calc_dsi(smiles)
    print(f"\n{name}:")
    print(f"  Peak Fiedler: {peak:.6f} | Final Fiedler: {final:.6f}")
    print(f"  DSI_RelDrop   (final * (1 - final/peak))  : {dsi_rel:.6f}")
    print(f"  DSI_AbsDrop   (final * (peak - final))    : {dsi_abs:.6f}")
    print(f"  DSI_LogDrop   (final * log(1+peak-final)) : {dsi_log:.6f}")
    print(f"  DSI_Entropy   (final * (1 - final/peak)^2): {dsi_ent:.6f}")
    print(f"  DSI_Ratio     (final^2 / peak)            : {dsi_ratio:.6f}")
    print(f"  DSI_Fluc      (final * std(post-peak))    : {dsi_fluc:.6f}")
