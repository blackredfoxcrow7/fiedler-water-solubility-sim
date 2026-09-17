#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import numpy as np
import networkx as nx
import plotly.graph_objects as go

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

import fiedler_optimization_sim_old as fsim

# Premium Color Map for Solute & Water Atoms
COLOR_MAP = {
    'S_O': '#FF5E5E',   # 溶質酸素 (鮮やかな赤)
    'S_C': '#4A4A4A',   # 溶質炭素 (ダークグレー)
    'S_H': '#E0E0E0',   # 溶質水素 (ホワイト)
    'W_O': '#4A90E2',   # 水の酸素 (鮮やかな青)
    'W_H': '#F5F5F5',   # 水の水素 (ホワイト)
}

def generate_interactive_html(f_history, frames_data, output_path):
    print("Generating 3D interactive HTML for tert-Butanol (Old Simulator)...")
    frames = []

    for t, (G_t, pos_t, _) in enumerate(frames_data):
        # 1. Sort nodes into categories for plotting
        node_x, node_y, node_z = [], [], []
        node_colors = []
        node_sizes = []
        node_text = []

        for n in G_t.nodes():
            pos = pos_t[n]
            sym = G_t.nodes[n].get('symbol', 'C')
            ntype = G_t.nodes[n].get('type', 'solute_atom')

            # Identify node type for coloring/sizing
            if ntype == 'solute_atom':
                key = f"S_{sym}"
                size = 12 if sym == 'O' else (10 if sym == 'C' else 4)
            else:
                key = f"W_{sym}"
                size = 8 if sym == 'O' else 3

            color = COLOR_MAP.get(key, '#FF00FF')
            
            node_x.append(pos[0])
            node_y.append(pos[1])
            node_z.append(pos[2])
            node_colors.append(color)
            node_sizes.append(size)

            category = G_t.nodes[n].get('category', 'water')
            node_text.append(
                f"<b>Node:</b> {n}<br>"
                f"<b>Element:</b> {sym}<br>"
                f"<b>Type:</b> {ntype}<br>"
                f"<b>Category:</b> {category}"
            )

        # 2. Extract edge lines based on type
        cov_x, cov_y, cov_z = [], [], []         # Covalent bonds (Gray)
        wnet_x, wnet_y, wnet_z = [], [], []       # Water network (Dark blue)
        solv_x, solv_y, solv_z = [], [], []       # Solvation (Light blue)
        opt_x, opt_y, opt_z = [], [], []         # Optimized (Lime green)

        for u, v, d in G_t.edges(data=True):
            p1, p2 = pos_t[u], pos_t[v]
            etype = d.get('edge_type')

            if etype in ['covalent', 'covalent_w']:
                cov_x.extend([p1[0], p2[0], None])
                cov_y.extend([p1[1], p2[1], None])
                cov_z.extend([p1[2], p2[2], None])
            elif etype == 'water_net':
                wnet_x.extend([p1[0], p2[0], None])
                wnet_y.extend([p1[1], p2[1], None])
                wnet_z.extend([p1[2], p2[2], None])
            elif etype == 'dummy_solvation':
                solv_x.extend([p1[0], p2[0], None])
                solv_y.extend([p1[1], p2[1], None])
                solv_z.extend([p1[2], p2[2], None])
            elif etype == 'opt_water':
                opt_x.extend([p1[0], p2[0], None])
                opt_y.extend([p1[1], p2[1], None])
                opt_z.extend([p1[2], p2[2], None])

        plot_data = [
            # Covalent Bonds (Gray)
            go.Scatter3d(x=cov_x, y=cov_y, z=cov_z, mode='lines', line=dict(color='#888888', width=3), name='Covalent'),
            # Water Network (Blue)
            go.Scatter3d(x=wnet_x, y=wnet_y, z=wnet_z, mode='lines', line=dict(color='#1F77B4', width=2), opacity=0.5, name='Water Net'),
            # Solvation (Light Blue)
            go.Scatter3d(x=solv_x, y=solv_y, z=solv_z, mode='lines', line=dict(color='#AEC7E8', width=2, dash='dash'), opacity=0.7, name='Solvation'),
            # Optimized bonds (Lime Green)
            go.Scatter3d(x=opt_x, y=opt_y, z=opt_z, mode='lines', line=dict(color='#7ED321', width=3, dash='dash'), opacity=0.8, name='AI Optimized'),
            # Atoms
            go.Scatter3d(x=node_x, y=node_y, z=node_z, mode='markers', marker=dict(size=node_sizes, color=node_colors, line=dict(color='#1E1E1E', width=1)), text=node_text, hoverinfo='text', name='Atoms')
        ]

        frames.append(go.Frame(data=plot_data, name=f'step{t}'))

    # Config slider and buttons
    sliders = [dict(
        active=0, yanchor="top", xanchor="left",
        currentvalue=dict(font=dict(size=14, color='white'), prefix="Step: ", visible=True, xanchor="right"),
        transition=dict(duration=100, easing="cubic-in-out"),
        pad=dict(b=10, t=50), len=0.9, x=0.1, y=0,
        steps=[dict(
            args=[[f.name], dict(frame=dict(duration=100, redraw=True), mode="immediate", transition=dict(duration=0))],
            label=str(k), method="animate"
        ) for k, f in enumerate(frames)]
    )]

    updatemenus = [dict(
        type="buttons", showactive=False, y=0, x=0, xanchor="right", yanchor="top", pad=dict(t=50, r=10),
        buttons=[
            dict(label="▶ Play", method="animate", args=[None, dict(frame=dict(duration=100, redraw=True), fromcurrent=True, transition=dict(duration=0))]),
            dict(label="⏸ Pause", method="animate", args=[[None], dict(frame=dict(duration=0, redraw=False), mode="immediate", transition=dict(duration=0))])
        ]
    )]

    fig = go.Figure(
        data=frames[0].data if frames else [],
        layout=go.Layout(
            title=dict(text="tert-Butanol Trajectory (Original Old Simulator - Full System)", font=dict(color='white', size=18)),
            scene=dict(xaxis=dict(visible=False), yaxis=dict(visible=False), zaxis=dict(visible=False), bgcolor='black'),
            paper_bgcolor='black', font=dict(color='white'), template='plotly_dark',
            updatemenus=updatemenus, sliders=sliders, margin=dict(l=0, r=0, b=0, t=60)
        ),
        frames=frames
    )

    fig.write_html(output_path)
    print(f"[SUCCESS] Interactive HTML generated at: {output_path}")

if __name__ == "__main__":
    smiles = "CC(C)(C)O"
    G_init, w_os = fsim.create_integrated_system([smiles], water_per_atom=3, is_loop=True)
    f_history, frames_data = fsim.run_fiedler_opt_sim(G_init, w_os, steps=100, allow_long_distance_water=False)
    
    output_file = "/home/eldenring/waterMain/tert_butanol_old_simulation.html"
    generate_interactive_html(f_history, frames_data, output_file)
