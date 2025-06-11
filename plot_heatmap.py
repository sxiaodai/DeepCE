import numpy as np
import os
from dataPreprocess_1 import *
import csv
import pandas as pd
import time

import matplotlib.pyplot as plt
from matplotlib import font_manager
font = font_manager.FontProperties(fname='C:/Users/86188/Desktop/simsun.ttf')
plt.rc('font', family=['Songti SC','Times New Roman'])

start = time.perf_counter()

# input
gene_pair_list_path='C:/Users/86188/Desktop/method/LSY/1000genes/GenePairList/mHSC-E-NonSpecific-GPL.csv'
gene_expression_path='C:/Users/86188/Desktop/method/scRNA-Seq/mHSC-E/mHSC-E/ExpressionDataOrdered.csv'
gold_network_path='C:/Users/86188/Desktop/method/Gold-network/Gold-network/Non-specific-ChIP-seq-network.csv'
# output
save_dir = 'C:/Users/86188/Desktop/method/LSY/1000genes/mHSC-E-Nonspecific-LogPearson_new_clr/WindowSize441-TF5-Target5-Lag64/FullMatrix_TF'
if not os.path.isdir(save_dir):
    os.makedirs(save_dir)

Window_size=441
TF_Gap=5
Target_Gap=5
Max_lag=64

# Load gene expression data
origin_expression_record,cells=get_normalized_expression_data(gene_expression_path)
origin_expression_record_df=pd.DataFrame(origin_expression_record)


# Load gold_pair_record
all_gene_list=[]
gold_pair_record={}
f_genePairList=open(gene_pair_list_path)

for single_pair in list(csv.reader(f_genePairList))[1:]:
    if single_pair[2]=='1':
        if single_pair[0] not in gold_pair_record:
            gold_pair_record[single_pair[0]]=[single_pair[1]]
        else:
            gold_pair_record[single_pair[0]].append(single_pair[1])
        # count all genes in gold edges
        if single_pair[0] not in all_gene_list:
            all_gene_list.append(single_pair[0])
        if single_pair[1] not in all_gene_list:
            all_gene_list.append(single_pair[1])
f_genePairList.close()
# print dataset statistics
print('All genes:'+str(len(all_gene_list)))
print('TFs:'+str(len(gold_pair_record.keys())))

# [Added Section]
############################################
# Gene expression smoothing via top-10 correlation averaging
# Create an empty DataFrame to store filtered expression data
filtered_expression_data = pd.DataFrame()

for gene in all_gene_list:
    if gene in origin_expression_record:
        gene_expression = origin_expression_record[gene].tolist()
        filtered_expression_data[gene] = gene_expression
    else:
        print(f"{gene} is not found")

filtered_expression_data.T.to_csv(os.path.join(save_dir, 'filtered_expression_data.csv'), index=True)

end = time.perf_counter()
print(f"time: {end - start:.2f}")

filtered_expression_data = pd.read_csv(os.path.join(save_dir, 'filtered_expression_data.csv'))

genes = filtered_expression_data.iloc[:, 0]
expression_data = filtered_expression_data.iloc[:, 1:]
expression_data.index = genes

print(expression_data)

correlation_matrix = expression_data.T.corr(method='pearson')
print(correlation_matrix)

top_genes = {}
for gene in correlation_matrix.index:
    sorted_genes = correlation_matrix[gene].sort_values(ascending=False)
    top_genes[gene] = sorted_genes.head(10).index.tolist()

def replace_with_average(expression_data, top_genes):
    new_expression_data = pd.DataFrame(index=expression_data.index, columns=expression_data.columns)
    for gene, genes_list in top_genes.items():
        selected_genes = expression_data.loc[genes_list]
        new_expression_data.loc[gene] = selected_genes.mean()
    return new_expression_data

new_expression_data = replace_with_average(expression_data, top_genes)

print(new_expression_data)

for gene in new_expression_data.index:
    if gene in origin_expression_record:
        origin_expression_record[gene] = new_expression_data.loc[gene].to_numpy()
    else:
        print(f"基因 {gene}is not found")
############################################


# Generate Pearson matrix
label_list=[]
pair_list=[]
total_matrix=[]
num_tf=-1 #为什么是-1啊
num_label1=0
num_label0=0

#regulate
#MCM2-CDC7
#MCM4-CDK1

## not regulate
#MCM2-ALG9
#MCM4-MYCN

tf_name = 'MCM2'
target_name = 'CDK1'
tf_data=origin_expression_record[tf_name]
target_data = origin_expression_record[target_name]
tf_list = list(tf_data)
target_list = list(target_data)
tf_data = tf_list
target_data = target_list

#shuffle
# combined_data = np.vstack([tf_list, target_list])  #
# shuffled_data = combined_data[:, np.random.permutation(combined_data.shape[1])]
# tf_data = shuffled_data[0, :]  # 第一行是 tf_data
# target_data = shuffled_data[1, :]

pair_matrix_1 = None
pair_matrix_2 = None

#control cell numbers by means of timepoints
timepoints=len(cells)

pair_matrix = None

for i in range(64):
    tf_window_left = i * TF_Gap
    if tf_window_left + Window_size + Target_Gap * (Max_lag - 1) > timepoints:
        break
    tf_window_data = tf_data[tf_window_left:tf_window_left + Window_size]
    single_tf_list = []
    tf_pd = pd.Series(tf_window_data)
    if (tf_window_data == np.zeros(Window_size)).all():
        exit('tf data error!')
        continue
    for lag in range(0, Max_lag):
        target_window_left = tf_window_left + lag * Target_Gap
        if target_window_left + Window_size > timepoints:
            exit('Planning error!')
        target_window_data = target_data[target_window_left:target_window_left + Window_size]
        target_pd = pd.Series(target_window_data)
        pcc = tf_pd.corr(target_pd, method="pearson")
        if np.isnan(pcc) or np.isinf(pcc):
            single_tf_list.append(0.0)
        else:
            single_tf_list.append(pcc)
    if pair_matrix is None:
        pair_matrix = np.array(single_tf_list).copy()
    else:
        pair_matrix = np.vstack((pair_matrix, np.array(single_tf_list)))

pair_matrix_new = np.zeros((64, 64))
for u in range(64):
    row_mean = np.mean(pair_matrix[u, :])
    row_std = np.std(pair_matrix[u, :])
    for v in range(64):
        col_mean = np.mean(pair_matrix[:, v])
        col_std = np.std(pair_matrix[:, v])
        # auv_u = (pair_matrix[u, v] - row_mean) / row_std
        # auv_v = (pair_matrix[u, v] - col_mean) / col_std
        auv_u = max(0, (pair_matrix[u, v] - row_mean) / row_std)
        auv_v = 0 if col_std == 0 else max(0, (pair_matrix[u, v] - col_mean) / col_std)

        auv_new = np.sqrt(auv_u ** 2 + auv_v ** 2)
        pair_matrix_new[u, v] = auv_new
# total_matrix.append(pair_matrix.copy())
total_matrix.append(pair_matrix_new.copy())

# paint heatmap
import numpy as np
import matplotlib.pyplot as plt

levels = np.linspace(0, 6.4, 20)
# levels = np.linspace(-0.15, 0.3, 20)
# levels = np.linspace(-0.4, 0.7, 20)
# levels = np.linspace(0, 4.5, 20)
# levels = np.linspace(0, 4.5, 20)
plt.figure(figsize=(6.4, 4.8))

contour_plot = plt.contourf(pair_matrix_new, cmap='bwr', levels=levels, vmin=0, vmax=6.4)

cbar = plt.colorbar(contour_plot)

cbar.set_ticks([0 ,0.8,1.6,2.4,3.2,4,4.8,5.6,6.4])
# cbar.set_ticks([0,0.5,1,1.5,2,2.5,3,3.5,4,4.5])
# cbar.set_ticks([-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7])
# cbar.set_ticks([-0.15,-0.1,-0.05,0,0.05,0.1,0.15,0.2,0.25,0.3])

cbar.ax.tick_params(labelsize=20)
plt.title(f'{tf_name} - {target_name}',fontsize=20)
plt.xlabel('target window',fontproperties=font,fontsize=20)
plt.ylabel('TF window',fontproperties=font,fontsize=20)

plt.gca().invert_yaxis()

plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
plt.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)

plt.show()