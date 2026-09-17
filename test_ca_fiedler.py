import numpy as np
import networkx as nx

# 1D セル・オートマトン (Elementary CA) の時空間グラフを構築し Fiedler 値を計算
def build_ca_spacetime_graph(rule_number, steps=30, width=61):
    # ルールを8ビットバイナリに変換
    rule_bin = [int(x) for x in f"{rule_number:08b}"]
    
    # 初期状態 (中央のみ1)
    state = np.zeros(width, dtype=int)
    state[width // 2] = 1
    
    history = [state.copy()]
    
    for _ in range(steps - 1):
        new_state = np.zeros(width, dtype=int)
        for i in range(width):
            # 周期境界条件
            left = state[(i - 1) % width]
            center = state[i]
            right = state[(i + 1) % width]
            neighborhood = (left << 2) | (center << 1) | right
            # 7 - neighborhood
            new_state[i] = rule_bin[7 - neighborhood]
        state = new_state
        history.append(state.copy())
        
    history = np.array(history)
    
    # ON状態(1)のセルをノードとし、時空間で隣接(近傍および次ステップ)するノード間にエッジを張る
    G = nx.Graph()
    active_cells = []
    for t in range(steps):
        for x in range(width):
            if history[t, x] == 1:
                node_id = f"{t}_{x}"
                G.add_node(node_id, t=t, x=x)
                active_cells.append((t, x))
                
    # エッジの張込み (空間隣接 ＆ 時間隣接)
    for t, x in active_cells:
        curr_id = f"{t}_{x}"
        # 同一時間ステップの隣
        if (t, (x + 1) % width) in active_cells:
            G.add_edge(curr_id, f"{t}_{(x + 1) % width}")
        # 次の時間ステップの同位置・両隣
        if t + 1 < steps:
            for dx in [-1, 0, 1]:
                next_x = (x + dx) % width
                if (t + 1, next_x) in active_cells:
                    G.add_edge(curr_id, f"{t+1}_{next_x}")
                    
    return G, history

# 代表的なCAルールでのFiedler値の計算テスト
rules_to_test = {
    "Rule 250 (均一/周期 - Class 1/2)": 250,
    "Rule 90 (シェルピンスキー・フラクタル - Class 3)": 90,
    "Rule 30 (カオス/ランダム - Class 3)": 30,
    "Rule 110 (複雑/計算普遍性 - Class 4)": 110
}

print("=========================================================================")
print(" 🧪 Stephen Wolfram Cellular Automata vs Graph Fiedler Value Test")
print("=========================================================================")

for name, rule in rules_to_test.items():
    G, _ = build_ca_spacetime_graph(rule, steps=30, width=61)
    if nx.is_connected(G) and len(G) > 2:
        L = nx.laplacian_matrix(G).toarray()
        evals = np.linalg.eigvalsh(L)
        fiedler = evals[1]
        density = nx.density(G)
        print(f"{name:<45} | ノード数={len(G):<4} | エッジ密度={density:.4f} | Fiedler値(λ2)={fiedler:.6f}")
    else:
        # 非連結の場合最大コンポーネントで計算
        components = list(nx.connected_components(G))
        largest = max(components, key=len)
        subG = G.subgraph(largest)
        L = nx.laplacian_matrix(subG).toarray()
        evals = np.linalg.eigvalsh(L)
        fiedler = evals[1]
        print(f"{name:<45} | ノード数={len(subG):<4} | (最大連結成分) | Fiedler値(λ2)={fiedler:.6f}")

print("=========================================================================")
