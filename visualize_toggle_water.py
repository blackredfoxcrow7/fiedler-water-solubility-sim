import sys
import os
import networkx as nx
import numpy as np
import plotly.graph_objects as go

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import fiedler_sim_hbond_offset_clash_free as fsim

def generate_toggle_3d_visualization(smiles, mol_name, steps=100):
    G_init, w_os = fsim.create_integrated_system([smiles], water_per_atom=3, is_loop=True)
    f_hist, frames_data, delta_f_final, _ = fsim.run_fiedler_opt_sim(G_init, w_os, steps=steps)
    
    # 最終フレームのグラフと座標を取得
    G_final, pos_final, evals_final = frames_data[-1]
    
    # ノードの分類
    solute_nodes = [n for n, d in G_final.nodes(data=True) if d.get('type') == 'solute_atom']
    water_nodes = [n for n, d in G_final.nodes(data=True) if d.get('type') == 'water_atom']
    
    fig = go.Figure()

    # --- Trace 0: 溶質原子ノード ---
    s_x, s_y, s_z, s_color, s_text = [], [], [], [], []
    for n in solute_nodes:
        x, y, z = pos_final[n]
        s_x.append(x); s_y.append(y); s_z.append(z)
        s_color.append(G_final.nodes[n].get('color', 'gray'))
        s_text.append(f"{n} ({G_final.nodes[n].get('symbol')})")
        
    fig.add_trace(go.Scatter3d(
        x=s_x, y=s_y, z=s_z, mode='markers+text',
        marker=dict(size=10, color=s_color, opacity=0.95, line=dict(color='black', width=1)),
        text=[G_final.nodes[n].get('symbol') for n in solute_nodes],
        textposition="top center", name="Solute Atoms", visible=True
    ))

    # --- Trace 1: 溶質の共有結合エッジ ---
    cov_x, cov_y, cov_z = [], [], []
    for u, v, d in G_final.edges(data=True):
        if d.get('edge_type') == 'covalent':
            x0, y0, z0 = pos_final[u]
            x1, y1, z1 = pos_final[v]
            cov_x.extend([x0, x1, None]); cov_y.extend([y0, y1, None]); cov_z.extend([z0, z1, None])

    fig.add_trace(go.Scatter3d(
        x=cov_x, y=cov_y, z=cov_z, mode='lines',
        line=dict(color='black', width=5), name='Solute Covalent Bonds', visible=True
    ))

    # --- Trace 2: 水分子原子ノード ---
    w_x, w_y, w_z, w_color, w_size = [], [], [], [], []
    for n in water_nodes:
        x, y, z = pos_final[n]
        w_x.append(x); w_y.append(y); w_z.append(z)
        w_color.append(G_final.nodes[n].get('color', 'blue'))
        w_size.append(7 if G_final.nodes[n].get('symbol') == 'O' else 3)

    fig.add_trace(go.Scatter3d(
        x=w_x, y=w_y, z=w_z, mode='markers',
        marker=dict(size=w_size, color=w_color, opacity=0.8),
        name="Water Atoms", visible=True
    ))

    # --- Trace 3: 水分子共有結合 & 水-水 ネットワークエッジ ---
    ww_x, ww_y, ww_z = [], [], []
    for u, v, d in G_final.edges(data=True):
        etype = d.get('edge_type')
        if (etype in ['water_net', 'opt_water', 'covalent_w']) and (u in water_nodes and v in water_nodes):
            x0, y0, z0 = pos_final[u]
            x1, y1, z1 = pos_final[v]
            ww_x.extend([x0, x1, None]); ww_y.extend([y0, y1, None]); ww_z.extend([z0, z1, None])

    fig.add_trace(go.Scatter3d(
        x=ww_x, y=ww_y, z=ww_z, mode='lines',
        line=dict(color='deepskyblue', width=4), name='Water-Water Network', visible=True
    ))

    # --- Trace 4: 溶質-水 水素結合エッジ ---
    sw_x, sw_y, sw_z = [], [], []
    for u, v, d in G_final.edges(data=True):
        etype = d.get('edge_type')
        if etype in ['dummy_solvation', 'opt_water'] and ((u in solute_nodes and v in water_nodes) or (v in solute_nodes and u in water_nodes)):
            x0, y0, z0 = pos_final[u]
            x1, y1, z1 = pos_final[v]
            sw_x.extend([x0, x1, None]); sw_y.extend([y0, y1, None]); sw_z.extend([z0, z1, None])

    fig.add_trace(go.Scatter3d(
        x=sw_x, y=sw_y, z=sw_z, mode='lines',
        line=dict(color='lime', width=3, dash='dash'), name='Solute-Water H-Bonds', visible=True
    ))

    # --- 3モード切替ボタン (Show All / Solute Only / Water Only) の配置 ---
    updatemenus = [
        dict(
            type="buttons",
            direction="left",
            buttons=[
                dict(
                    label="🌐 すべて表示 (Show All)",
                    method="update",
                    args=[{"visible": [True, True, True, True, True]},
                          {"title": f"Final State (All): {mol_name} | Δf_final: {delta_f_final:.4f}"}]
                ),
                dict(
                    label="🧪 溶質のみ (Solute Only)",
                    method="update",
                    args=[{"visible": [True, True, False, False, False]},
                          {"title": f"Final State (Solute Only): {mol_name} | Δf_final: {delta_f_final:.4f}"}]
                ),
                dict(
                    label="💧 水分子のみ (Water Only)",
                    method="update",
                    args=[{"visible": [False, False, True, True, False]},
                          {"title": f"Final State (Water Network Only): {mol_name} | Δf_final: {delta_f_final:.4f}"}]
                )
            ],
            pad={"r": 10, "t": 10},
            showactive=True,
            x=0.02,
            xanchor="left",
            y=1.12,
            yanchor="top",
            font=dict(color="black", size=12)
        )
    ]

    fig.update_layout(
        title=f"Final State (All): {mol_name} ({smiles}) | Δf_final: {delta_f_final:.4f}",
        updatemenus=updatemenus,
        scene=dict(
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            zaxis=dict(visible=False),
            bgcolor='rgb(15, 23, 42)' # シックなダーク背景
        ),
        margin=dict(l=0, r=0, b=0, t=80)
    )

    out_html = f"/home/eldenring/waterMain/toggle_3d_{mol_name.lower()}.html"
    fig.write_html(out_html)
    print(f"  ➜ 3モード切替付き 3D HTML 出力完了: {out_html}")
    return out_html

# 再実行
generate_toggle_3d_visualization("CC(C)(C)O", "tert_butanol")
generate_toggle_3d_visualization("CCOCC", "diethyl_ether")
generate_toggle_3d_visualization("CCCCO", "1_butanol")
