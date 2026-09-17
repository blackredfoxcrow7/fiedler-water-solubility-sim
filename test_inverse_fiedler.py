import numpy as np
import networkx as nx

# 逆問題 (Inverse Problem): 目的の巨視的機能 (Target Fiedler Value λ2) から、
# 微視的ルール (どの原子ペアのエッジ A_ij を追加・削除すべきか) をFiedler感度定理で逆算・設計するテスト

def inverse_fiedler_design(target_fiedler=0.25, num_nodes=10, max_iter=20):
    # 初期状態: 10ノードの単純な一重鎖 (伸長ペプチド鎖)
    G = nx.path_graph(num_nodes)
    
    print(f"🎯 目標の巨視的機能 (Target λ2) : {target_fiedler:.6f}")
    
    for iteration in range(max_iter):
        L = nx.laplacian_matrix(G).toarray().astype(float)
        evals, evecs = np.linalg.eigh(L)
        
        curr_fiedler = evals[1]
        fiedler_vec = evecs[:, 1]
        
        print(f"Iter {iteration:2d} | 現在の Fiedler値 (λ2) = {curr_fiedler:.6f} | ノード数={len(G)} | エッジ数={G.number_of_edges()}")
        
        if abs(curr_fiedler - target_fiedler) < 0.005:
            print("✨ 逆計算成功！目標の巨視的機能に到達しました！")
            break
            
        # Fiedler感度定理 (Fiedler Sensitivity Theorem):
        # エッジ (i, j) を追加したときの λ2 の変化量 Δλ2 ≈ (v_i - v_j)^2
        # エッジ (i, j) を削除したときの λ2 の変化量 Δλ2 ≈ -(v_i - v_j)^2
        
        best_edge = None
        best_delta = 0
        action = None
        
        if curr_fiedler < target_fiedler:
            # λ2 を高めたい ➔ (v_i - v_j)^2 が最も大きい非存在エッジを追加するルールを逆算
            max_diff = -1
            for i in range(num_nodes):
                for j in range(i + 1, num_nodes):
                    if not G.has_edge(i, j):
                        diff = (fiedler_vec[i] - fiedler_vec[j]) ** 2
                        if diff > max_diff:
                            max_diff = diff
                            best_edge = (i, j)
            action = "add"
        else:
            # λ2 を下げたい ➔ (v_i - v_j)^2 が最も大きい既存エッジを削除するルールを逆算
            max_diff = -1
            for i, j in G.edges():
                diff = (fiedler_vec[i] - fiedler_vec[j]) ** 2
                if diff > max_diff:
                    max_diff = diff
                    best_edge = (i, j)
            action = "remove"
            
        if best_edge is None:
            break
            
        if action == "add":
            G.add_edge(*best_edge)
        else:
            G.remove_edge(*best_edge)
            
    return G

print("=========================================================================")
print(" 🔬 Inverse Fiedler Design: Macro Function -> Micro Rule Inverse Engine")
print("=========================================================================")
G_designed = inverse_fiedler_design(target_fiedler=0.25, num_nodes=10)
print("=========================================================================")
