import numpy as np
import os
from dataPreprocess import *
import csv
import pandas as pd
import time

# calculate running time
start = time.perf_counter()

# input
gene_pair_list_path='/home/qqwu/paper/LSY/1000genes/GenePairList/mHSC-E-NonSpecific-GPL.csv'
gene_expression_path='/home/qqwu/paper/scRNA-Seq/mHSC-E/mHSC-E/ExpressionDataOrdered.csv'
gold_network_path='/home/qqwu/paper/Gold-network/Gold-network/Non-specific-ChIP-seq-network.csv'
# output
save_dir = '/home/qqwu/paper/LSY/1000genes/mHSC-E-Nonspecific-LogPearson/WindowSize441-TF5-Target5-Lag64/FullMatrix_TF'
if not os.path.isdir(save_dir):
    os.makedirs(save_dir)

Window_size=441
TF_Gap=5
Target_Gap=5
Max_lag=64

# Load gene expression data
origin_expression_record,cells=get_normalized_expression_data(gene_expression_path)

# Load gold_pair_record
all_gene_list=[]
gold_pair_record={}
f_genePairList=open(gene_pair_list_path)


# Read each gene pair from the gene pair list
for single_pair in list(csv.reader(f_genePairList))[1:]:
    # If this is a positive sample (label == '1')
    if single_pair[2]=='1':
        # Update the gold_pair_record dictionary: TF -> list of target genes
        if single_pair[0] not in gold_pair_record:
            gold_pair_record[single_pair[0]]=[single_pair[1]]
        else:
            gold_pair_record[single_pair[0]].append(single_pair[1])
        # Record all genes involved in gold (positive) edges
        if single_pair[0] not in all_gene_list:
            all_gene_list.append(single_pair[0])
        if single_pair[1] not in all_gene_list:
            all_gene_list.append(single_pair[1])
f_genePairList.close()
# Print dataset statistics
print('All genes:'+str(len(all_gene_list)))
print('TFs:'+str(len(gold_pair_record.keys())))

# [Added Section]
############################################
# Gene expression smoothing via top-10 correlation averaging
# Create an empty DataFrame to store filtered expression data
filtered_expression_data = pd.DataFrame()

# Extract expression data for genes in all_gene_list from origin_expression_record
for gene in all_gene_list:
    if gene in origin_expression_record:
        gene_expression = origin_expression_record[gene].tolist()
        filtered_expression_data[gene] = gene_expression
    else:
        print(f"{gene} is not found")

# Save the filtered expression data (transpose: genes as rows) to a CSV file
filtered_expression_data.T.to_csv(os.path.join(save_dir, 'filtered_expression_data.csv'), index=True)

# Record end time and print the elapsed time
end = time.perf_counter()
print(f"time: {end - start:.2f}")

# Read back the saved expression data
filtered_expression_data = pd.read_csv(os.path.join(save_dir, 'filtered_expression_data.csv'))

# Set gene names as index and remove the first column
genes = filtered_expression_data.iloc[:, 0]
expression_data = filtered_expression_data.iloc[:, 1:]
expression_data.index = genes

# Compute Pearson correlation between all gene pairs
correlation_matrix = expression_data.T.corr(method='pearson')
print(correlation_matrix)

# For each gene, find the 10 genes with the highest correlation
top_genes = {}
for gene in correlation_matrix.index:
    sorted_genes = correlation_matrix[gene].sort_values(ascending=False)
    top_genes[gene] = sorted_genes.head(10).index.tolist()

# Replace each gene's expression values with the average of its top 10 most correlated genes
def replace_with_average(expression_data, top_genes):
    new_expression_data = pd.DataFrame(index=expression_data.index, columns=expression_data.columns)
    for gene, genes_list in top_genes.items():
        selected_genes = expression_data.loc[genes_list]
        new_expression_data.loc[gene] = selected_genes.mean()
    return new_expression_data
new_expression_data = replace_with_average(expression_data, top_genes)

# Output the updated expression data
print(new_expression_data)

# Replace the original expression values in origin_expression_record with the new averaged values
for gene in new_expression_data.index:
    if gene in origin_expression_record:
        origin_expression_record[gene] = new_expression_data.loc[gene].to_numpy()
    else:
        print(f"{gene}is not found")
############################################


# Generate Pearson matrix
label_list=[]
pair_list=[]
total_matrix=[]
num_tf=-1
num_label1=0
num_label0=0

#control cell numbers by means of timepoints
timepoints=len(cells)

pair_matrix = None

for i in gold_pair_record:
    num_tf+=1
    for j in range(len(all_gene_list)):
        print ('Generating matrix of gene pair '+str(num_tf)+' '+str(j))
        tf_name=i
        target_name= all_gene_list[j]
        if tf_name in gold_pair_record and target_name in gold_pair_record[tf_name]:
            label=1
            num_label1+=1
        else:
            label=0
            num_label0+=1
        label_list.append(label)
        pair_list.append(tf_name + ',' + target_name)
        tf_data=origin_expression_record[tf_name]
        target_data = origin_expression_record[target_name]

        # Calculate matrix for one gene pair based on sliding windows and lag correlations
        pair_matrix = None
        for tf_window_left in range(0, timepoints - Window_size + 1, TF_Gap):
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
                #pcc = tf_pd.corr(target_pd, method="spearman")
                #pcc = tf_pd.corr(target_pd, method="kendall")
                if np.isnan(pcc) or np.isinf(pcc):
                    single_tf_list.append(0.0)
                else:
                    single_tf_list.append(pcc)
            if pair_matrix is None:
                pair_matrix = np.array(single_tf_list).copy()
            else:
                pair_matrix = np.vstack((pair_matrix, np.array(single_tf_list)))
        pair_matrix_new = np.zeros((64,64))
        # [Added Function]
        ############################################
        # Normalize the pair_matrix using row-wise and column-wise statistics
        # to emphasize significant correlations and suppress noise.
        for u in range(64): 
            row_mean = np.mean(pair_matrix[u, :])
            row_std = np.std(pair_matrix[u, :])
            for v in range(64):
                col_mean = np.mean(pair_matrix[:, v])
                col_std = np.std(pair_matrix[:, v])
                #auv_u = (pair_matrix[u, v] - row_mean) / row_std
                #auv_v = (pair_matrix[u, v] - col_mean) / col_std
                auv_u = max(0, (pair_matrix[u, v] - row_mean) / row_std)
                auv_v = 0 if col_std == 0 else max(0, (pair_matrix[u, v] - col_mean) / col_std)

                auv_new = np.sqrt(auv_u ** 2 + auv_v ** 2)
                pair_matrix_new[u, v] = auv_new
        # total_matrix.append(pair_matrix.copy())
        total_matrix.append(pair_matrix_new.copy())
        ############################################

end = time.perf_counter()
print('Running time:'+str(end-start))


if (len(total_matrix)>0):
    total_matrix = np.array(total_matrix)[:, :, :, np.newaxis]
else:
    exit('Save error.')
np.save(save_dir +'/matrix.npy', total_matrix)
np.save(save_dir+'/label.npy', np.array(label_list))
np.save(save_dir+'/gene_pair.npy', np.array(pair_list))
print('PCC matrix generation finish.')
print('Positive edges:'+str(num_label1))
print('Negative edges:'+str(num_label0))
print('Density='+str(num_label1/(num_label1+num_label0)))

