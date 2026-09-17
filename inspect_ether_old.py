import sys
import os
import numpy as np
import networkx as nx

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

import fiedler_optimization_sim_old as fsim

smiles = "CCOCC" # Diethyl Ether
print(f"Running old simulation for Diethyl Ether ({smiles})...")

G_init, w_os = fsim.create_integrated_system([smiles], water_per_atom=3, is_loop=True)
f_history, frames_data = fsim.run_fiedler_opt_sim(G_init, w_os, steps=100, allow_long_distance_water=False)

# Get the final frame
G_final, pos_final, evals_final = frames_data[-1]

# Find the oxygen node
o_nodes = [n for n, d in G_final.nodes(data=True) if d.get('symbol') == 'O' and d.get('type') == 'solute_atom']
print(f"\nSolute Oxygen nodes found: {o_nodes}")

for o_node in o_nodes:
    neighbors = list(G_final.neighbors(o_node))
    covalent_neighbors = [n for n in neighbors if G_final.edges[o_node, n].get('edge_type') == 'covalent']
    solvation_neighbors = [n for n in neighbors if G_final.edges[o_node, n].get('edge_type') in ['opt_water', 'dummy_solvation']]
    
    print(f"\n--- Analysis of Oxygen Node: {o_node} ---")
    print(f"  Covalent connections: {covalent_neighbors}")
    print(f"  Water/Solvation connections: {solvation_neighbors}")
    print(f"  Total degree in final graph: {G_final.degree(o_node)}")
    
    print("\n  Details of connected water molecules:")
    for w in solvation_neighbors:
        dist = np.linalg.norm(pos_final[o_node] - pos_final[w])
        # Find other neighbors of this water molecule to see if it bridges the network
        w_neighbors = list(G_final.neighbors(w))
        print(f"    Water Node {w}: distance = {dist:.4f}, degree = {G_final.degree(w)}, other connections = {w_neighbors}")

print("\n--- Solute Node Categories & Degrees ---")
s_nodes = [n for n, d in G_final.nodes(data=True) if d.get('type') == 'solute_atom']
for s in s_nodes:
    category = G_final.nodes[s].get('category')
    symbol = G_final.nodes[s].get('symbol')
    degree = G_final.degree(s)
    pos = pos_final[s]
    print(f"  Node {s:<10} ({symbol:<2}, {category:<12}): Degree = {degree}, Coords = [{pos[0]:.4f}, {pos[1]:.4f}, {pos[2]:.4f}]")

print(f"\nFinal Fiedler value: {f_history[-1]:.6f}")
