import numpy as np
from scipy.optimize import minimize
import random

c3_data = {
    '1-Propanol': {'init': 0.182022, 'peak': 0.182022, 'final': 0.049296},
    '2-Propanol': {'init': 0.170477, 'peak': 0.176527, 'final': 0.057929},
    'Methoxyethane': {'init': 0.178808, 'peak': 0.181276, 'final': 0.074263}
}

c4_data = {
    '1-Butanol': {'init': 0.145123, 'peak': 0.145123, 'final': 0.033329},
    'Isobutanol': {'init': 0.149859, 'peak': 0.149859, 'final': 0.036259},
    '2-Butanol': {'init': 0.145899, 'peak': 0.145899, 'final': 0.042233},
    'tert-Butanol': {'init': 0.145802, 'peak': 0.145802, 'final': 0.047211},
    'Diethyl Ether': {'init': 0.141299, 'peak': 0.149233, 'final': 0.061575}
}

# Features: [init, peak, final, peak/init, final/init, peak/final, final/peak]
def get_features(d):
    init, peak, final = d['init'], d['peak'], d['final']
    return np.array([
        init,
        peak,
        final,
        peak / init,
        final / init,
        peak / final,
        final / peak
    ])

features = {k: get_features(d) for k, d in {**c3_data, **c4_data}.items()}

# We want:
# C3: 2-Propanol > 1-Propanol > Methoxyethane
# C4: tert-Butanol > 2-Butanol > Isobutanol > 1-Butanol > Diethyl Ether

constraints = [
    ('2-Propanol', '1-Propanol'),
    ('1-Propanol', 'Methoxyethane'),
    ('tert-Butanol', '2-Butanol'),
    ('2-Butanol', 'Isobutanol'),
    ('Isobutanol', '1-Butanol'),
    ('1-Butanol', 'Diethyl Ether')
]

# Define loss function to maximize margins
def loss_func(w):
    loss = 0.0
    for high, low in constraints:
        diff = np.dot(w, features[high]) - np.dot(w, features[low])
        # We want diff > 0.1. Penalize if diff < 0.1.
        if diff < 0.1:
            loss += (0.1 - diff)**2 + 10.0 * (diff < 0)  # Heavy penalty for wrong order
    # L2 regularization on weights
    loss += 0.01 * np.sum(w**2)
    return loss

# Run optimization from multiple random starting points
best_w = None
best_loss = float('inf')

random.seed(42)
for _ in range(500):
    w0 = np.array([random.uniform(-5, 5) for _ in range(7)])
    res = minimize(loss_func, w0, method='BFGS')
    if res.fun < best_loss:
        best_loss = res.fun
        best_w = res.x

w = best_w
s = {k: np.dot(w, features[k]) for k in features}

# Count satisfied
satisfied = 0
for high, low in constraints:
    if s[high] > s[low]:
        satisfied += 1

print(f"Optimal weights found! Loss: {best_loss:.6f}, Satisfied: {satisfied}/6")
print("Weights:", [f"{val:.4f}" for val in w])
print("Scores:")
print("  C3:")
for k in c3_data:
    print(f"    {k:<15}: {s[k]:.4f}")
print("  C4:")
for k in c4_data:
    print(f"    {k:<15}: {s[k]:.4f}")
