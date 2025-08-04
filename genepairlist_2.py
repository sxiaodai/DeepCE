from dataPreprocess_1 import *
import os
#input
gene_expression_path = 'C:/Users/86188/Desktop/method/scRNA-Seq/mHSC-E/mHSC-E/ExpressionDataOrdered.csv' #现有
gene_order_path = 'C:/Users/86188/Desktop/method/scRNA-Seq/mHSC-E/mHSC-E/GeneOrdering.csv' #现有
gold_network_path='C:/Users/86188/Desktop/method/Gold-network/Gold-network/Non-Specific-ChIP-seq-network.csv' #现有
#output
filtered_path= 'C:/Users/86188/Desktop/method/LSY/1000genes/FilteredGoldNetwork/'
FGN_file_name= 'mHSC-E-NonSpecific-FGN.csv'
rank_path= 'C:/Users/86188/Desktop/method/LSY/1000genes/'
Rank_file_name='mHSC-E-1000genes-NonSpecific-rank.csv'
genePairList_path= 'C:/Users/86188/Desktop/method/LSY/1000genes/GenePairList/'
GPL_file_name= 'mHSC-E-NonSpecific-GPL.csv'
Rank_num=1000

# 1. Reads and filters original gene expression data.
origin_expression_record,cells=get_origin_expression_data(gene_expression_path)
Expression_gene_num=len(origin_expression_record)
Expression_cell_num=len(cells)

# 2. Identifies lowly expressed genes (expressed in ≤10% of cells).
low_express_gene_list=get_low_express_gene(origin_expression_record,len(cells))
print(str(len(low_express_gene_list))+' genes in low expression.')
for gene in low_express_gene_list:
    origin_expression_record.pop(gene)

# 3. Ranks genes based on ordering file and selects top-k.
if not os.path.isdir(rank_path):
    os.makedirs(rank_path)
rank_list,significant_gene_list = get_gene_ranking(gene_order_path, low_express_gene_list,Rank_num, rank_path + Rank_file_name, False)

# 4. Filters the gold standard network to retain only ranked genes.
gold_pair_record, gold_score_record ,unique_gene_list= get_filtered_gold(gold_network_path, rank_list,rank_path+Rank_file_name,True)
if not os.path.isdir(filtered_path):
    os.makedirs(filtered_path)
generate_filtered_gold( gold_pair_record, gold_score_record,filtered_path + FGN_file_name)

# 5. Generates gene pair list.
if not os.path.isdir(genePairList_path):
    os.makedirs(genePairList_path)
get_gene_pair_list(unique_gene_list, gold_pair_record, gold_score_record, genePairList_path + GPL_file_name)




