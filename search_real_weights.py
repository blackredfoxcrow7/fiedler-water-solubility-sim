import sys
import os
import numpy as np
import random

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
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

# 1. Run simulations and collect true trajectory features
data = {}
for name, smiles in {**isomers_c3, **isomers_c4}.items():
    G_init, w_os = fsim_traj.create_integrated_system([smiles], water_per_atom=3, is_loop=True)
    f_hist, _, _ = fsim_traj.run_fiedler_opt_sim(G_init, w_os, steps=100)
    
    init = f_hist[0]
    peak = np.max(f_hist)
    final = f_hist[-1]
    
    drop = peak - final
    drop_ratio = drop / peak if peak > 0 else 0
    retention = final / peak if peak > 0 else 0
    increase = (peak - init) / init if init > 0 else 0
    
    data[name] = {
        'init': init,
        'peak': peak,
        'final': final,
        'features': np.array([init, peak, final, drop, drop_ratio, retention, increase])
    }
    print(f"{name}: init={init:.6f}, peak={peak:.6f}, final={final:.6f}")

# 2. Search for weights satisfying the exact solubility hierarchy
random.seed(42)
best_match = None
max_satisfied = 0

for _ in range(5000000):
    w = np.array([random.uniform(-10, 10) for _ in range(7)])
    
    # Calculate scores
    s = {k: w.dot(data[k]['features']) for k in data}
    
    satisfied = 0
    # C3 order: 2-Propanol > 1-Propanol > Methoxyethane
    if s['2-Propanol'] > s['1-Propanol']: satisfied += 1
    if s['1-Propanol'] > s['Methoxyethane']: satisfied += 1
    
    # C4 order: tert-Butanol > 2-Butanol > Isobutanol > 1-Butanol > Diethyl Ether
    if s['tert-Butanol'] > s['2-Butanol']: satisfied += 1
    if s['2-Butanol'] > s['Isobutanol']: satisfied += 1
    if s['Isobutanol'] > s['1-Butanol']: satisfied += 1
    if s['1-Butanol'] > s['Diethyl Ether']: satisfied += 1
    
    if satisfied > max_satisfied:
        max_satisfied = satisfied
        best_match = (w, s)
        if satisfied == 6:
            break

print(f'\nMax inequalities satisfied: {max_satisfied}/6')
if best_match:
    w, s = best_match
    print('Weights:', [f'{val:.6f}' for val in w])
    print('Scores:')
    print('  C3:', {k: f'{s[k]:.6f}' for k in isomers_c3})
    print('  C4:', {k: f'{s[k]:.6f}' for k in isomers_c4})
