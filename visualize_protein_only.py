#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
visualize_protein_only.py

このスクリプトは、fiedler_optimization_sim_old.py のシミュレーションをそのまま実行し、
水分子（溶媒）を取り除いた「入力タンパク質（溶質アトム）部分のみ」を抽出して
黒背景に映えるプレミアムカラー（最初の改良コードのデザイン）で3D可視化します。
（幾何的修正や力場クリーンアップを入れる前の、オリジナルに忠実な挙動です）
"""

import os
import sys
import numpy as np
import networkx as nx
import plotly.graph_objects as go
import matplotlib

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

# 原子別のプレミアムカラーマップ
PREMIUM_COLORS = {
    'H': '#F0F0F0',   # 水素 (薄いグレー)
    'C': '#4A4A4A',   # 炭素 (ダークグレー)
    'O': '#FF5E5E',   # 酸素 (鮮やかな赤)
    'N': '#4A90E2',   # 窒素 (鮮やかな青)
    'P': '#F5A623',   # リン (オレンジ)
    'S': '#F8E71C',   # 硫黄 (イエロー)
    'F': '#7ED321',   # フッ素 (ライトグリーン)
    'Cl': '#417505',  # 塩素 (グリーン)
    'Br': '#2E5B02',  # 臭素 (ダークグリーン)
    'Na': '#9013FE',  # ナトリウム (パープル)
    'K': '#BD10E0'    # カリウム (バイオレット)
}

def extract_and_visualize_protein_only(f_history, frames_data):
    """
    シミュレーション結果から溶質（タンパク質）部分のみを抽出し、
    黒背景に映えるプレミアムカラー（グレーの結合線、太い結合）で3Dアニメーション化します。
    """
    print("\nタンパク質部分の抽出と3D可視化モデルの構築を開始します...")
    frames = []
    
    for t, (G_t, pos_t, _) in enumerate(frames_data):
        # 1. 溶質（タンパク質）アトムノードのみを抽出
        solute_nodes = [n for n, d in G_t.nodes(data=True) if d.get('type') == 'solute_atom']
        
        if not solute_nodes:
            continue
            
        node_x = [pos_t[n][0] for n in solute_nodes]
        node_y = [pos_t[n][1] for n in solute_nodes]
        node_z = [pos_t[n][2] for n in solute_nodes]
        
        node_colors = []
        node_sizes = []
        node_text = []
        
        for n in solute_nodes:
            symbol = G_t.nodes[n].get('symbol', 'C')
            color = PREMIUM_COLORS.get(symbol, G_t.nodes[n].get('color', '#FF00FF'))
            size = G_t.nodes[n].get('size', 10)
            
            # 見やすさ調整（水素アトムをやや認識しやすく）
            adjusted_size = size * 1.5 if symbol == 'H' else size
            
            node_colors.append(color)
            node_sizes.append(adjusted_size)
            
            # ホバー時の詳細情報
            category = G_t.nodes[n].get('category', 'unknown')
            charge = G_t.nodes[n].get('charge', 0.0)
            mol_id = G_t.nodes[n].get('molecule_id', 'MOL')
            node_text.append(
                f"<b>Node:</b> {n}<br>"
                f"<b>Element:</b> {symbol}<br>"
                f"<b>Category:</b> {category}<br>"
                f"<b>Charge:</b> {charge:.3f}<br>"
                f"<b>Molecule:</b> {mol_id}"
            )

        # 2. 溶質ノード同士を繋ぐエッジのみを抽出
        edge_x, edge_y, edge_z = [], [], []
        agg_x, agg_y, agg_z = [], [], []
        opt_x, opt_y, opt_z = [], [], []
        
        for u, v, d in G_t.edges(data=True):
            if u in solute_nodes and v in solute_nodes:
                p1, p2 = pos_t[u], pos_t[v]
                etype = d.get('edge_type')
                
                if etype == 'agg':  # 疎水性凝集結合 (オレンジの点線)
                    agg_x.extend([p1[0], p2[0], None])
                    agg_y.extend([p1[1], p2[1], None])
                    agg_z.extend([p1[2], p2[2], None])
                elif etype == 'opt_water':  # 最適化された相互作用結合 (黄緑の点線)
                    opt_x.extend([p1[0], p2[0], None])
                    opt_y.extend([p1[1], p2[1], None])
                    opt_z.extend([p1[2], p2[2], None])
                else:  # コバレント結合など (グレーの実線)
                    edge_x.extend([p1[0], p2[0], None])
                    edge_y.extend([p1[1], p2[1], None])
                    edge_z.extend([p1[2], p2[2], None])

        # プロット用データの構築
        plot_data = [
            # 共有結合 (黒背景に映える太いグレーの線)
            go.Scatter3d(
                x=edge_x, y=edge_y, z=edge_z,
                mode='lines',
                line=dict(color='#888888', width=3),
                opacity=0.7,
                name='Covalent Bond',
                hoverinfo='skip'
            ),
            # 疎水性凝集 (オレンジの点線)
            go.Scatter3d(
                x=agg_x, y=agg_y, z=agg_z,
                mode='lines',
                line=dict(color='#FFA500', width=3, dash='dash'),
                opacity=0.8,
                name='Hydrophobic Aggregation',
                hoverinfo='skip'
            ),
            # AI最適化結合 (黄緑の点線)
            go.Scatter3d(
                x=opt_x, y=opt_y, z=opt_z,
                mode='lines',
                line=dict(color='#7ED321', width=3, dash='dash'),
                opacity=0.8,
                name='AI Optimized Bond',
                hoverinfo='skip'
            ),
            # 原子
            go.Scatter3d(
                x=node_x, y=node_y, z=node_z,
                mode='markers',
                marker=dict(
                    size=node_sizes,
                    color=node_colors,
                    opacity=0.95,
                    line=dict(color='#1E1E1E', width=1.5)
                ),
                text=node_text,
                hoverinfo='text',
                name='Atoms'
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

    # フィギュア全体のレイアウト設定（黒背景・文字白色）
    fig = go.Figure(
        data=frames[0].data if frames else [],
        layout=go.Layout(
            title=dict(
                text="Protein-Only 3D Folding Simulation (Original Physics)",
                font=dict(color='white', size=18)
            ),
            scene=dict(
                xaxis=dict(visible=False),
                yaxis=dict(visible=False),
                zaxis=dict(visible=False),
                bgcolor='black'  # 3D空間の背景色を黒に設定
            ),
            paper_bgcolor='black',  # 周辺全体の背景色を黒に設定
            font=dict(color='white'),  # 文字色を白に設定
            template='plotly_dark',  # ダークテーマの適用
            updatemenus=updatemenus,
            sliders=sliders,
            margin=dict(l=0, r=0, b=0, t=60)
        ),
        frames=frames
    )

    # HTMLファイルとして保存
    html_filename = "protein_only_simulation.html"
    fig.write_html(html_filename)
    print(f"[SUCCESS] インタラクティブな3Dモデル（HTML形式）を保存しました:\n          {os.path.abspath(html_filename)}")
    
    # 表示試行
    try:
        fig.show()
    except Exception as e:
        print("[NOTE] 環境により図の自動起動はスキップされました。作成された上記のHTMLファイルをブラウザで直接開いてください。")

if __name__ == "__main__":
    # ユーザーが指定したトリプトファン含有 of Peptide sequence (Gly-Tyr-Asp-Pro-Glu-Thr-Gly-Thr-Trp-Gly)H-Gly-Tyr-Asp-Pro-Glu-Thr-Gly-Thr-Trp-Gly-OH
    smiles_list = ["NCC(=O)N[C@@H](Cc1ccc(cc1)O)C(=O)N[C@@H](CC(=O)O)C(=O)N2CCC[C@H]2C(=O)N[C@@H](CCC(=O)O)C(=O)N[C@@H]([C@@H](C)O)C(=O)NCC(=O)N[C@@H]([C@@H](C)O)C(=O)N[C@@H](Cc3c[nH]c4ccccc34)C(=O)NCC(=O)O"]
    #  Chignolin
    print("=" * 60)
    print(" 🧪 Fiedler-driven Topology Optimization Simulation (Protein Only)")
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
    
    # 2. Fiedler値のグラフ保存 (Matplotlib)
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
    
    plot_filename = "fiedler_evolution.png"
    plt.savefig(plot_filename)
    print(f"[SUCCESS] 統計グラフ画像を保存しました: {os.path.abspath(plot_filename)}")
    
    # 3. タンパク質部分のみを抽出して3D Plotly可視化
    extract_and_visualize_protein_only(f_history, frames_data)
