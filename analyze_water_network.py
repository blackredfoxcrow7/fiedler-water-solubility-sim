import sys
import os
import networkx as nx
import numpy as np
import plotly.graph_objects as go

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import fiedler_sim_hbond_offset as fsim

def analyze_and_visualize_water_network(smiles, mol_name, steps=100):
    G_init, w_os = fsim.create_integrated_system([smiles], water_per_atom=3, is_loop=True)
    f_hist, frames_data, delta_f_final, _ = fsim.run_fiedler_opt_sim(G_init, w_os, steps=steps)
    
    # 最終フレームのグラフと座標を取得
    G_final, pos_final, evals_final = frames_data[-1]
    
    # ノードの分類
    solute_nodes = [n for n, d in G_final.nodes(data=True) if d.get('type') == 'solute_atom']
    water_nodes = [n for n, d in G_final.nodes(data=True) if d.get('type') == 'water_atom']
    water_o_nodes = [n for n in water_nodes if G_final.nodes[n].get('symbol') == 'O']
    
    # 水分子のみからなるサブグラフ (Water-only Subgraph)
    G_water_only = G_final.subgraph(water_nodes).copy()
    
    # 水ネットワークの定量的解析
    num_water_components = nx.number_connected_components(G_water_only)
    water_degrees = [G_water_only.degree(n) for n in water_o_nodes]
    avg_water_degree = np.mean(water_degrees) if water_degrees else 0
    
    # 溶質と結合している水酸素の数 (Solvation Water Count)
    bound_waters = set()
    for u, v, d in G_final.edges(data=True):
        if (u in solute_nodes and v in water_nodes) or (v in solute_nodes and u in water_nodes):
            w_node = u if u in water_nodes else v
            bound_waters.add(w_node)
            
    num_bound_waters = len(bound_waters)
    num_free_waters = len(water_o_nodes) - num_bound_waters
    
    print(f"\n==================================================")
    print(f" 💧 Water Network Analysis: {mol_name} ({smiles})")
    print(f"==================================================")
    print(f"  ・総水分子数 (O)            : {len(water_o_nodes)}")
    print(f"  ・溶質に結合中の水分子数     : {num_bound_waters} ({num_bound_waters/len(water_o_nodes)*100:.1f}%)")
    print(f"  ・フリー/排除された水分子数  : {num_free_waters} ({num_free_waters/len(water_o_nodes)*100:.1f}%)")
    print(f"  ・水クラスタの独立コンポーネント数: {num_water_components}")
    print(f"  ・水酸素間の平均結合次数    : {avg_water_degree:.2f}")
    print(f"  ・実効水和変化量 (Δf_final) : {delta_f_final:.6f}")

    # --- 3D インタラクティブ Plotly HTML の作成 (水ネットワーク強調表示) ---
    fig = go.Figure()

    # 1. 溶質原子ノードの描画
    s_x, s_y, s_z, s_color, s_text = [], [], [], [], []
    for n in solute_nodes:
        x, y, z = pos_final[n]
        s_x.append(x); s_y.append(y); s_z.append(z)
        s_color.append(G_final.nodes[n].get('color', 'gray'))
        s_text.append(f"{n} ({G_final.nodes[n].get('symbol')})")
        
    fig.add_trace(go.Scatter3d(
        x=s_x, y=s_y, z=s_z, mode='markers+text',
        marker=dict(size=10, color=s_color, opacity=0.9, line=dict(color='black', width=1)),
        text=[G_final.nodes[n].get('symbol') for n in solute_nodes],
        textposition="top center", name="Solute Atoms"
    ))

    # 2. 水ノード（酸素/水素）の描画
    w_x, w_y, w_z, w_color, w_size = [], [], [], [], []
    for n in water_nodes:
        x, y, z = pos_final[n]
        w_x.append(x); w_y.append(y); w_z.append(z)
        w_color.append(G_final.nodes[n].get('color', 'blue'))
        w_size.append(7 if G_final.nodes[n].get('symbol') == 'O' else 3)

    fig.add_trace(go.Scatter3d(
        x=w_x, y=w_y, z=w_z, mode='markers',
        marker=dict(size=w_size, color=w_color, opacity=0.8),
        name="Water Network Atoms"
    ))

    # 3. エッジの描画 (エッジの種類ごとに色分け)
    # A. 水-水 ネットワーク結合 (シアン太線)
    ww_x, ww_y, ww_z = [], [], []
    # B. 溶質-水 水和結合 (緑破線)
    sw_x, sw_y, sw_z = [], [], []
    # C. 共有結合 (黒太線)
    cov_x, cov_y, cov_z = [], [], []

    for u, v, d in G_final.edges(data=True):
        x0, y0, z0 = pos_final[u]
        x1, y1, z1 = pos_final[v]
        etype = d.get('edge_type')
        
        if etype in ['water_net', 'opt_water'] and u in water_nodes and v in water_nodes:
            ww_x.extend([x0, x1, None]); ww_y.extend([y0, y1, None]); ww_z.extend([z0, z1, None])
        elif etype in ['dummy_solvation', 'opt_water'] and ((u in solute_nodes and v in water_nodes) or (v in solute_nodes and u in water_nodes)):
            sw_x.extend([x0, x1, None]); sw_y.extend([y0, y1, None]); sw_z.extend([z0, z1, None])
        elif etype in ['covalent', 'covalent_w']:
            cov_x.extend([x0, x1, None]); cov_y.extend([y0, y1, None]); cov_z.extend([z0, z1, None])

    # 共有結合
    fig.add_trace(go.Scatter3d(x=cov_x, y=cov_y, z=cov_z, mode='lines', line=dict(color='black', width=4), name='Covalent Bonds'))
    # 水-水 結合
    fig.add_trace(go.Scatter3d(x=ww_x, y=ww_y, z=ww_z, mode='lines', line=dict(color='deepskyblue', width=5), name='Water-Water Network'))
    # 溶质-水 結合
    fig.add_trace(go.Scatter3d(x=sw_x, y=sw_y, z=sw_z, mode='lines', line=dict(color='lime', width=3, dash='dash'), name='Solute-Water H-Bonds'))

    fig.update_layout(
        title=f"Final Water Network Topology: {mol_name} ({smiles})<br>Δf_final: {delta_f_final:.4f} | Bound Waters: {num_bound_waters}/{len(water_o_nodes)}",
        scene=dict(
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            zaxis=dict(visible=False),
            bgcolor='rgb(15, 23, 42)' # 宇宙/真空を表現するシックなダークブルー背景
        ),
        margin=dict(l=0, r=0, b=0, t=50)
    )

    out_html = f"/home/eldenring/waterMain/water_network_{mol_name.lower()}.html"
    fig.write_html(out_html)
    print(f"  ➜ 3D 可視化 HTML 保存完了: {out_html}")
    return out_html

# 代表的な3つの分子で解析と可視化 HTML 生成を実行
analyze_and_visualize_water_network("CC(C)(C)O", "tert_butanol")
analyze_and_visualize_water_network("CCOCC", "diethyl_ether")
analyze_and_visualize_water_network("CCCCO", "1_butanol")
