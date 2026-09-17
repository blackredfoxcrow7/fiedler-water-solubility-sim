import sys
import os

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import visualize_toggle_water as vis

# 提案書 Step 1 の注目分子について 3D 切替可視化 HTML を一括生成
vis.generate_toggle_3d_visualization("OCC(O)CO", "glycerin")
vis.generate_toggle_3d_visualization("C1CCC(O)CC1", "cyclohexanol")
vis.generate_toggle_3d_visualization("CCCCCCO", "1_hexanol")
