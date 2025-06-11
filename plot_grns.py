import pandas as pd

# file1 = pd.read_csv(r"C:\Users\86188\Desktop\method\Gold-network\Gold-network\Non-Specific-ChIP-seq-network.csv", usecols=[0, 1], header=None)
# file2 = pd.read_csv(r"C:\Users\86188\Desktop\method\Gold-network\Gold-network\mHSC-ChIP-seq-network.csv", usecols=[0, 1], header=None)
# file3 = pd.read_csv(r"C:\Users\86188\Desktop\method\Gold-network\Gold-network\STRING-network.csv", usecols=[0, 1], header=None)

# file1 = pd.read_csv(r"C:\Users\86188\Desktop\method\Gold-network\Gold-network\human\Non-specific-ChIP-seq-network.csv", usecols=[0, 1], header=None)
# file2 = pd.read_csv(r"C:\Users\86188\Desktop\method\Gold-network\Gold-network\human\HepG2-ChIP-seq-network.csv", usecols=[0, 1], header=None)
# file3 = pd.read_csv(r"C:\Users\86188\Desktop\method\Gold-network\Gold-network\human\STRING-network.csv", usecols=[0, 1], header=None)

# file1_pairs = set(tuple(sorted([row[0], row[1]])) for row in file1.values)
# file2_pairs = set(tuple(sorted([row[0], row[1]])) for row in file2.values)
# file3_pairs = set(tuple(sorted([row[0], row[1]])) for row in file3.values)
#
# file = pd.read_excel(r"C:\Users\86188\Desktop\method\predict\hHep-500.xlsx", usecols=[0, 1], header=None)
#
# file_pairs = set(tuple(sorted([row[0], row[1]])) for row in file.values)
#
# match_file1 = len(file_pairs.intersection(file1_pairs))
# match_file2 = len(file_pairs.intersection(file2_pairs))
# match_file3 = len(file_pairs.intersection(file3_pairs))
#
# print(f"{match_file1}")
# print(f"{match_file2}")
# print(f"{match_file3}")
#
#
# matches = []
#
# for _, row in file.iterrows():
#     pair = tuple(sorted([row[0], row[1]]))
#     if pair in file1_pairs:
#         matches.append(1)
#     elif pair in file2_pairs:
#         matches.append(2)
#     elif pair in file3_pairs:
#         matches.append(3)
#     else:
#         matches.append(0)
#
# file['Matched'] = matches
#
# file.to_excel(r"C:\Users\86188\Desktop\method\predict\hHep-500_with_match.xlsx", index=False,header=False)
#
# print("done")

######################
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from collections import Counter
import numpy as np

import matplotlib.pyplot as plt
from matplotlib import font_manager
font = font_manager.FontProperties(fname='C:/Users/86188/Desktop/simsun.ttf')
plt.rc('font', family=['Songti SC','Times New Roman'])

file = pd.read_excel(r"C:\Users\86188\Desktop\method\predict\hHep-500_with_match.xlsx", header=None)

G = nx.DiGraph()

tf_counts = Counter(file[0])

tf_color = '#FF6347'
target_color = '#32CD32'

# 添加节点和边
for _, row in file.iterrows():
    tf = row[0]
    target = row[1]
    match = row[2]

    if not G.has_node(tf):
        G.add_node(tf, type='TF', size=tf_counts[tf] * 100, color=tf_color)


    if not G.has_node(target):

        if target in tf_counts:
            G.add_node(target, type='TF', size=tf_counts[target] * 100, color=tf_color)
        else:
            G.add_node(target, type='Target', size=200, color=target_color)


    if match == 1:
        edge_color = '#FFA07A'
    elif match == 2:
        edge_color = '#20B2AA'
    elif match == 3:
        edge_color = '#9370DB'
    else:
        edge_color = '#A9A9A9'
    G.add_edge(tf, target, color=edge_color)


tf_nodes = [node for node, data in G.nodes(data=True) if data['type'] == 'TF']
target_nodes = [node for node, data in G.nodes(data=True) if data['type'] == 'Target']

# 定义圆盘半径
radius_tf = 1.0
radius_target = 0.8


pos = {}
for i, tf_node in enumerate(tf_nodes):
    angle = 2 * np.pi * i / len(tf_nodes)
    pos[tf_node] = (radius_tf * np.cos(angle), radius_tf * np.sin(angle))


for i, target_node in enumerate(target_nodes):
    angle = 2 * np.pi * i / len(target_nodes)
    pos[target_node] = (radius_target * np.cos(angle), radius_target * np.sin(angle))

# 绘图
node_sizes = [G.nodes[node]['size'] for node in G.nodes]
node_colors = [G.nodes[node]['color'] for node in G.nodes]
edge_colors = [G[u][v]['color'] for u, v in G.edges]


font_size = 4
font_size_tf = 36

plt.figure(figsize=(18, 18))
nx.draw(G, pos, with_labels=False, node_size=node_sizes, node_color=node_colors, edge_color=edge_colors,
        font_size=font_size, font_weight='bold', width=2, alpha=0.7, arrowstyle='-|>', arrowsize=20)  # 增加箭头（方向）

# 获取TF节点的标签并设置更大的字体
for node, (x, y) in pos.items():
    if G.nodes[node]['type'] == 'TF':
        plt.text(x, y, node, fontsize=font_size_tf, ha='center', va='center')
    else:
        plt.text(x, y, node, fontsize=font_size, ha='center', va='center')

plt.title("Gene Regulatory Network with Dispersed TF Nodes in Circular Layout", fontsize=18)

plt.show()

#################################
