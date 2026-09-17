import numpy as np
import networkx as nx

# ロジスティックマップ x_{n+1} = r * x_n * (1 - x_n) から リカレンス・ネットワークを構築
def build_logistic_recurrence_network(r_param, steps=300, eps=0.05, transient=100):
    x = 0.5
    # 過渡応答を捨てる
    for _ in range(transient):
        x = r_param * x * (1 - x)
        
    trajectory = []
    for _ in range(steps):
        x = r_param * x * (1 - x)
        trajectory.append(x)
        
    trajectory = np.array(trajectory)
    
    # 距離 ε 以内の点を結合してリカレンス・ネットワーク (Recurrence Network) を作成
    G = nx.Graph()
    for i in range(steps):
        G.add_node(i, val=trajectory[i])
        
    for i in range(steps):
        for j in range(i + 1, steps):
            if abs(trajectory[i] - trajectory[j]) < eps:
                G.add_edge(i, j)
                
    return G, trajectory

r_values = {
    "r = 3.2 (周期1 / リミットサイクル)": 3.2,
    "r = 3.5 (周期4 / 周期倍化分岐)": 3.5,
    "r = 3.83 (カオスの中の周期3の窓)": 3.83,
    "r = 4.00 (完全カオス / カオス状態)": 4.00
}

print("=========================================================================")
print(" 🌀 Deterministic Chaos (Logistic Map) vs Recurrence Network Fiedler Value")
print("=========================================================================")

for name, r in r_values.items():
    G, _ = build_logistic_recurrence_network(r, steps=300, eps=0.05)
    
    # 最大連結成分で Fiedler 値を算出
    components = list(nx.connected_components(G))
    largest = max(components, key=len)
    subG = G.subgraph(largest)
    
    L = nx.laplacian_matrix(subG).toarray()
    evals = np.linalg.eigvalsh(L)
    fiedler = evals[1]
    avg_deg = np.mean([d for _, d in subG.degree()])
    
    print(f"{name:<35} | 結節点数={len(subG):<3} | 平均次数={avg_deg:<5.2f} | Fiedler値(λ2)={fiedler:.6f}")

print("=========================================================================")
