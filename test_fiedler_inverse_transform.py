import numpy as np
import networkx as nx
import sys
import os

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import fiedler_dimer_interaction_optimizer as fdio

def test_direct_pinv(name, smiles):
    mol, coords = fdio.get_3d_conformer(smiles)
    N = mol.GetNumAtoms()
    G = nx.Graph()
    for i in range(N): G.add_node(i)
    for b in mol.GetBonds(): G.add_edge(b.GetBeginAtomIdx(), b.GetEndAtomIdx(), weight=1.0)
    L_norm = nx.normalized_laplacian_matrix(G).toarray()
    evals, evecs = np.linalg.eigh(L_norm)
    v = evecs[:, 1]
    
    X_c = coords - np.mean(coords, axis=0)
    
    # 順変換: d_phase = X_c^T * v   ((3, N) * (N,) -> (3,))
    d_phase = np.dot(X_c.T, v)
    
    # 逆変換: v_recon = (X_c^T)^+ * d_phase   (pinv of (3, N) is (N, 3))
    # pinv(X_c.T) is the exact Moore-Penrose pseudo-inverse!
    v_recon = np.dot(np.linalg.pinv(X_c.T), d_phase)
    
    # 再順変換: d_reconstructed = X_c^T * v_recon
    d_recon = np.dot(X_c.T, v_recon)
    
    err = np.linalg.norm(d_phase - d_recon)
    corr = np.corrcoef(v, v_recon)[0, 1]
    
    print(f"{name:<25} | d_phase Error: {err:.2e} | v Correlation (r): {corr:.4f}")
    return err < 1e-12

for name, smiles in [
    ("チロシン (Tyr)", "NC(Cc1ccc(O)cc1)C(=O)O"),
    ("テレフタル酸 (平面分子)", "O=C(O)c1ccc(C(=O)O)cc1"),
    ("1-ブタノール", "CCCCO"),
    ("フマル酸 (平面分子)", "O=C(O)/C=C/C(=O)O")
]:
    test_direct_pinv(name, smiles)

