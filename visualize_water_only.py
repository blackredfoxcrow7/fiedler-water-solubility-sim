#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
visualize_water_only.py

このスクリプトは、fiedler_optimization_sim_old.py のシミュレーションをそのまま実行し、
「水分子（溶媒）部分」をメインに抽出して3Dインタラクティブグラフィックス（Plotly）で黒背景に可視化します。
（幾何的修正や力場クリーンアップを適用する前の、オリジナルに忠実な挙動です）

タンパク質（溶質アトム）は、水分子の配置が分かりやすいよう、背景にうっすらと影（低不透明度）として描画されます。
その際、トリプトファン（Trp）やチロシン（Tyr）の原子は薄いピンクや紫で色分けされ、水との相対配置が視覚的に分かります。
"""

import os
import sys
import numpy as np
import networkx as nx
import random
import plotly.graph_objects as go
import matplotlib
from rdkit import Chem

# ディスプレイのない環境（GUIなし）でも動作するようにMatplotlibのバックエンドを調整
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# インポートパスにスクリプトのあるディレクトリを追加
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

try:
    import fiedler_optimization_sim_old as fsim
except ImportError:
    print("Error: fiedler_optimization_sim_old.py が見つかりません。同じディレクトリに配置してください。")
    sys.exit(1)

def extract_and_visualize_water_only(f_history, frames_data, smiles_list):
    """
    シミュレーション結果から水分子部分のみを抽出し、
    タンパク質の影（トリプトファン・チロシン着色）を背景に配置した3Dアニメーションを作成します。
    """
    print("\n水分子部分の抽出と3D可視化モデルの構築を開始します...")
    
    # RDKitを用いた特定アミノ酸側鎖の原子インデックス抽出（トポロジーマッチング）
    trp_atom_indices = set()
    tyr_atom_indices = set()
    
    for mol_idx, smiles in enumerate(smiles_list):
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        
        # トリプトファンのインドール環 (SMARTS)
        trp_pattern = Chem.MolFromSmarts("c1c[nH]c2ccccc12")
        trp_matches = mol.GetSubstructMatches(trp_pattern)
        for match in trp_matches:
            for idx in match:
                trp_atom_indices.add(f"S_{mol_idx}_{idx}")
                
        # チロシンのフェノール環 (SMARTS)
        tyr_pattern = Chem.MolFromSmarts("c1ccc(O)cc1")
        tyr_matches = mol.GetSubstructMatches(tyr_pattern)
        for match in tyr_matches:
            for idx in match:
                tyr_atom_indices.add(f"S_{mol_idx}_{idx}")
                
    frames = []
    
    for t, (G_t, pos_t, _) in enumerate(frames_data):
        # --- ノードの分類抽出 ---
        water_nodes = [n for n, d in G_t.nodes(data=True) if d.get('type') == 'water_atom']
        solute_nodes = [n for n, d in G_t.nodes(data=True) if d.get('type') == 'solute_atom']

        # 1. 水分子ノードのプロットデータ
        o_nodes = [n for n in water_nodes if G_t.nodes[n].get('symbol') == 'O']
        h_nodes = [n for n in water_nodes if G_t.nodes[n].get('symbol') == 'H']

        o_x = [pos_t[n][0] for n in o_nodes]
        o_y = [pos_t[n][1] for n in o_nodes]
        o_z = [pos_t[n][2] for n in o_nodes]
        
        h_x = [pos_t[n][0] for n in h_nodes]
        h_y = [pos_t[n][1] for n in h_nodes]
        h_z = [pos_t[n][2] for n in h_nodes]

        # 2. 水関連エッジの抽出
        cov_w_x, cov_w_y, cov_w_z = [], [], []  # 水分子内の共有結合 (O-H)
        net_w_x, net_w_y, net_w_z = [], [], []  # 初期水ネットワーク結合
        opt_w_x, opt_w_y, opt_w_z = [], [], []  # 最適化された水-水水素結合
        
        for u, v, d in G_t.edges(data=True):
            if u in water_nodes and v in water_nodes:
                p1, p2 = pos_t[u], pos_t[v]
                etype = d.get('edge_type')
                
                if etype == 'covalent_w':  # 水分子内共有結合 (細い薄グレー線)
                    cov_w_x.extend([p1[0], p2[0], None])
                    cov_w_y.extend([p1[1], p2[1], None])
                    cov_w_z.extend([p1[2], p2[2], None])
                elif etype == 'water_net':  # 水の初期連なり (シアン)
                    net_w_x.extend([p1[0], p2[0], None])
                    net_w_y.extend([p1[1], p2[1], None])
                    net_w_z.extend([p1[2], p2[2], None])
                elif etype == 'opt_water':  # 最適化された水素結合 (黄緑の点線)
                    opt_w_x.extend([p1[0], p2[0], None])
                    opt_w_y.extend([p1[1], p2[1], None])
                    opt_w_z.extend([p1[2], p2[2], None])

        # 3. 背景参照用：タンパク質アトム（トリプトファンとチロシンを着色）
        solute_x = [pos_t[n][0] for n in solute_nodes]
        solute_y = [pos_t[n][1] for n in solute_nodes]
        solute_z = [pos_t[n][2] for n in solute_nodes]
        
        solute_colors = []
        for n in solute_nodes:
            symbol = G_t.nodes[n].get('symbol', 'C')
            if n in trp_atom_indices and symbol == 'C':
                solute_colors.append('#FF1744')  # トリプトファン側鎖（赤・ピンク）
            elif n in tyr_atom_indices and symbol == 'C':
                solute_colors.append('#D500F9')  # チロシン側鎖（紫）
            else:
                solute_colors.append('#888888')  # 通常の原子（グレー）

        solute_edge_x, solute_edge_y, solute_edge_z = [], [], []
        for u, v, d in G_t.edges(data=True):
            if u in solute_nodes and v in solute_nodes:
                p1, p2 = pos_t[u], pos_t[v]
                solute_edge_x.extend([p1[0], p2[0], None])
                solute_edge_y.extend([p1[1], p2[1], None])
                solute_edge_z.extend([p1[2], p2[2], None])

        # プロット用データの構築
        plot_data = [
            # 背景：タンパク質の結合影 (やや見えやすくするため不透明度を調整)
            go.Scatter3d(
                x=solute_edge_x, y=solute_edge_y, z=solute_edge_z,
                mode='lines',
                line=dict(color='#666666', width=1.5),
                opacity=0.15,
                name='Protein Backbone (Shadow)',
                hoverinfo='skip'
            ),
            # 背景：タンパク質原子影（特定アミノ酸のカラー化）
            go.Scatter3d(
                x=solute_x, y=solute_y, z=solute_z,
                mode='markers',
                marker=dict(size=4.5, color=solute_colors, opacity=0.22),
                name='Protein Atoms (Shadow)',
                hoverinfo='skip'
            ),
            # 水分子内の共有結合 (O-H)
            go.Scatter3d(
                x=cov_w_x, y=cov_w_y, z=cov_w_z,
                mode='lines',
                line=dict(color='#E5E5E5', width=2),
                opacity=0.6,
                name='Water O-H Bond',
                hoverinfo='skip'
            ),
            # 水分子間のネットワーク結合 (water_net)
            go.Scatter3d(
                x=net_w_x, y=net_w_y, z=net_w_z,
                mode='lines',
                line=dict(color='#00E5FF', width=2),
                opacity=0.5,
                name='Water Net Loop',
                hoverinfo='skip'
            ),
            # AI最適化された水-水水素結合 (opt_water)
            go.Scatter3d(
                x=opt_w_x, y=opt_w_y, z=opt_w_z,
                mode='lines',
                line=dict(color='#7ED321', width=3, dash='dash'),
                opacity=0.8,
                name='AI Optimized H-Bond (W-W)',
                hoverinfo='skip'
            ),
            # 水素原子ノード (H)
            go.Scatter3d(
                x=h_x, y=h_y, z=h_z,
                mode='markers',
                marker=dict(size=4, color='#F0F0F0', opacity=0.9),
                name='Hydrogen (Water)',
                hoverinfo='name'
            ),
            # 酸素原子ノード (O)
            go.Scatter3d(
                x=o_x, y=o_y, z=o_z,
                mode='markers',
                marker=dict(
                    size=8,
                    color='#00A8FF',
                    opacity=0.95,
                    line=dict(color='#111111', width=1)
                ),
                name='Oxygen (Water)',
                hoverinfo='name'
            )
        ]

        frames.append(go.Frame(data=plot_data, name=f'step{t}'))

    # スライダーの挙動設定
    sliders = [dict(
        active=0,
        yanchor="top",
        xanchor="left",
        currentvalue=dict(font=dict(size=14, color='white'), prefix="Simulation Step: ", visible=True, xanchor="right"),
        transition=dict(duration=100, easing="cubic-in-out"),
        pad=dict(b=10, t=50),
        len=0.9,
        x=0.1,
        y=0,
        steps=[dict(
            args=[[f.name], dict(frame=dict(duration=100, redraw=True), mode="immediate", transition=dict(duration=0))],
            label=str(k),
            method="animate"
        ) for k, f in enumerate(frames)]
    )]

    # 再生と一時停止のボタン
    updatemenus = [dict(
        type="buttons",
        showactive=False,
        y=0,
        x=0,
        xanchor="right",
        yanchor="top",
        pad=dict(t=50, r=10),
        buttons=[
            dict(
                label="▶ Play",
                method="animate",
                args=[None, dict(frame=dict(duration=100, redraw=True), fromcurrent=True, transition=dict(duration=0))]
            ),
            dict(
                label="⏸ Pause",
                method="animate",
                args=[[None], dict(frame=dict(duration=0, redraw=False), mode="immediate", transition=dict(duration=0))]
            )
        ]
    )]

    # フィギュア全体のレイアウト設定
    fig = go.Figure(
        data=frames[0].data if frames else [],
        layout=go.Layout(
            title=dict(
                text="Water-Only 3D Folding Simulation (Original Physics)",
                font=dict(size=18, color='white')
            ),
            scene=dict(
                xaxis=dict(visible=False),
                yaxis=dict(visible=False),
                zaxis=dict(visible=False),
                bgcolor='black'
            ),
            paper_bgcolor='black',
            font=dict(color='white'),
            template='plotly_dark',
            updatemenus=updatemenus,
            sliders=sliders,
            margin=dict(l=0, r=0, b=0, t=60)
        ),
        frames=frames
    )

    # HTMLファイルとしても保存し、ブラウザで自在に開けるようにする
    html_filename = "water_only_simulation.html"
    fig.write_html(html_filename)
    print(f"\n[SUCCESS] インタラクティブな3Dモデル（HTML形式）を保存しました:\n          {os.path.abspath(html_filename)}")
    
    # GUI表示（headless環境等でエラーが発生した場合はスキップ）
    try:
        fig.show()
    except Exception as e:
        print("[NOTE] 環境により図の自動起動はスキップされました。作成された上記のHTMLファイルをブラウザで直接開いてください。")

if __name__ == "__main__":
    # ユーザーが指定したトリプトファン含有のペプチド配列 (Gly-Tyr-Asp-Pro-Glu-Thr-Gly-Thr-Trp-Gly)
    smiles_list = ["NCC(=O)N[C@@H](Cc1ccc(cc1)O)C(=O)N[C@@H](CC(=O)O)C(=O)N2CCC[C@H]2C(=O)N[C@@H](CCC(=O)O)C(=O)N[C@@H]([C@@H](C)O)C(=O)NCC(=O)N[C@@H]([C@@H](C)O)C(=O)N[C@@H](Cc3c[nH]c4ccccc34)C(=O)NCC(=O)O"]
    
    print("=" * 60)
    print(" 🧪 Fiedler-driven Topology Optimization Simulation (Water Only)")
    print("=" * 60)
    print(f"入力SMILES (ペプチド): {smiles_list[0][:80]}...")
    
    # シミュレーション設定
    ALLOW_LONG_DISTANCE = False
    
    # 1. シミュレーションの実行 (fiedler_optimization_sim_old.py のオリジナル関数をそのまま呼び出す)
    print("\n[1/3] システム構築とシミュレーションを開始します...")
    # 計算時間短縮のためデフォルトの配位水数を 2 に設定（4 に増やすと水和殻が密になりますが計算時間が伸びます）
    G_init, w_os = fsim.create_integrated_system(smiles_list, water_per_atom=2, is_loop=True)
    f_history, frames_data = fsim.run_fiedler_opt_sim(G_init, w_os, steps=100, allow_long_distance_water=ALLOW_LONG_DISTANCE)
    print("[SUCCESS] シミュレーションが完了しました。")
    
    # 2. グラフ保存 (Matplotlib)
    print("\n[2/3] Fiedler最大化による安定化過程グラフを作成中...")
    best_step = np.argmax(f_history)
    _, _, best_evals = frames_data[best_step]

    plt.figure(figsize=(10, 4))
    plt.subplot(1, 2, 1)
    plt.plot(f_history, color='blue', linewidth=2)
    plt.axvline(x=10, color='gray', linestyle=':', label='Optimization Start')
    plt.axvline(x=best_step, color='red', linestyle='--', label=f'Peak Stability (Step {best_step})')
    plt.title("Evolution by Fiedler Maximization")
    plt.xlabel("Step")
    plt.ylabel("Fiedler Value")
    plt.legend()
    plt.grid(True)
    
    plt.subplot(1, 2, 2)
    plt.hist(best_evals, bins=80, color='purple', alpha=0.7)
    plt.title(f"Eigenvalue Spikes (Best State at Step {best_step})")
    plt.xlabel("Eigenvalue")
    plt.ylabel("Frequency")
    plt.grid(True)
    
    plt.tight_layout()
    
    plot_filename = "fiedler_evolution_water.png"
    plt.savefig(plot_filename)
    print(f"[SUCCESS] 統計グラフ画像を保存しました: {os.path.abspath(plot_filename)}")
    
    # 3. 水分子部分のみを抽出して3D Plotly可視化
    extract_and_visualize_water_only(f_history, frames_data, smiles_list)
