# This file uses code from the DGRNS project.
# Original repository: https://github.com/guofei-tju/DGRNS
# Some modifications and adaptations have been made for this work.

import pandas as pd
import numpy as np
import csv
import os
import random
import seaborn as sns
import matplotlib.pyplot as plt


# Load a list of transcription factors (TFs) from a CSV file.
# Assumes the first column contains TF names and skips the header row.
def get_tf_list(tf_path):
    f_tf = open(tf_path)
    tf_reader = list(csv.reader(f_tf))
    tf_list=[]
    for single in tf_reader[1:]:
        tf_list.append(single[0])
    print('Load '+str(len(tf_list))+' TFs successfully!')
    return tf_list


# Load original gene expression data from a CSV file.
# Each row corresponds to a gene and each column (after the first) to a cell.
# Returns a dictionary mapping gene names to expression values and a list of cell names.
# Exits if duplicate gene names are found.
def get_origin_expression_data(gene_expression_path):
    f_expression = open(gene_expression_path)
    expression_reader = list(csv.reader(f_expression))
    cells = expression_reader[0][1:]
    num_cells = len(cells)
    expression_record = {}
    num_genes = 0
    for single_expression_reader in expression_reader[1:]:
        if single_expression_reader[0] in expression_record:
            exit('Gene name '+single_expression_reader[0]+' repeat!')
        expression_record[single_expression_reader[0]] = list(map(float, single_expression_reader[1:]))
        num_genes += 1
    print(str(num_genes) + ' genes and ' + str(num_cells) + ' cells are included in origin expression data.')
    return expression_record,cells


# Load original gene expression data from a CSV file.
# Returns a dictionary mapping each gene to its expression values across cells, and a list of cell names.
# Checks for duplicate gene names and exits if any are found.
def get_normalized_expression_data(gene_expression_path):
    expression_record,cells=get_origin_expression_data(gene_expression_path)
    expression_matrix = np.zeros((len(expression_record), len(cells)))
    index_row=0
    for gene in expression_record:
        expression_record[gene]=np.log10(np.array(expression_record[gene])+10**-2)
        expression_matrix[index_row]=expression_record[gene]
        index_row+=1

    return expression_record, cells

############################################
# [Added Function]
# Load and normalize gene expression data from a CSV file.
# Applies log10 transformation to each gene's expression values: log10(x + 1e-2 + 1).
# Returns the normalized expression data as a dictionary and the list of cell names.
def get_normalized_expression_data_new(gene_expression_path):
    expression_record,cells=get_origin_expression_data(gene_expression_path)
    expression_matrix = np.zeros((len(expression_record), len(cells)))
    index_row=0
    for gene in expression_record:
        expression_record[gene]=np.log10(np.array(expression_record[gene])+1+10**-2)
        expression_matrix[index_row]=expression_record[gene]
        index_row+=1

    return expression_record, cells
############################################

# Rank genes based on their variance from a given CSV file.
# Filters out genes with p-value >= 0.01 or in the low expression gene list.
# Detects and handles genes with duplicated variance values.
# If `flag` is True, writes the top `gene_num` ranked genes to `output_path`.
# Returns the list of top-ranked genes and all significant genes after filtering.
def get_gene_ranking(gene_order_path,low_express_gene_list,gene_num,output_path,flag):
    f_order = open(gene_order_path)
    order_reader = list(csv.reader(f_order))
    if flag:
        f_rank = open(output_path, 'w', newline='\n')
        f_rank_writer = csv.writer(f_rank)
    variance_record = {}
    variance_list = []
    significant_gene_list=[]
    for single_order_reader in order_reader[1:]:
        if float(single_order_reader[1]) >= 0.01:
            continue
        if single_order_reader[0] in low_express_gene_list:
            continue
        variance = float(single_order_reader[2])
        if variance not in variance_record:
            variance_record[variance] = single_order_reader[0]
        else:
            print(str(variance_record[variance]) + ' and ' + single_order_reader[0] + ' variance repeat!')
            variance_record[variance]=[variance_record[variance]]
            variance_record[variance].append(single_order_reader[0])
        variance_list.append(variance)
        significant_gene_list.append(single_order_reader[0])
    print('After delete genes with p-value>=0.01 or low expression, '+str(len(variance_list))+' genes left.')
    variance_list.sort(reverse=True)

    variance_range = {
        "max": max(variance_list),
        "min": min(variance_list),
        "mean": sum(variance_list) / len(variance_list),
        "std": (sum([(x - sum(variance_list) / len(variance_list)) ** 2 for x in variance_list]) / len(
            variance_list)) ** 0.5,
    }

    print(
        f"Variance Range: Max = {variance_range['max']}, Min = {variance_range['min']}, Mean = {variance_range['mean']}, Std = {variance_range['std']}")

    duplicate_genes_count = sum([len(v) - 1 for v in variance_record.values() if isinstance(v, list)])
    print(f"Number of duplicate genes (variance repeated): {duplicate_genes_count}")

    gene_rank = []
    for single_variance_list in variance_list[0:gene_num]:
        if type(variance_record[single_variance_list]) is str:
            gene_rank.append(variance_record[single_variance_list])
        else:
            gene_rank.append(variance_record[single_variance_list][0])
            del variance_record[single_variance_list][0]
            if len(variance_record[single_variance_list])==1:
                variance_record[single_variance_list]=variance_record[single_variance_list][0]
        if flag:
            f_rank_writer.writerow([variance_record[single_variance_list]])
    f_order.close()
    if flag:
        f_rank.close()
    print(len(gene_rank))
    return gene_rank,significant_gene_list


# Filter the gold standard gene regulatory network based on a ranked gene list.
# Only retains gene pairs where both genes are present in the `rank_list`.
# Calculates and prints statistics on the number of TFs, edges, gene overlap, network density,
# average target gene diversity, and target gene regulation count.
# If `flag` is True, writes the list of unique genes involved in the filtered network to `output_path`.
# Returns:
# - gold_pair_record: dict mapping TFs to their target genes
# - gold_score_record: dict mapping gene pairs to their associated scores (or default score if not provided)
# - unique_gene_list: list of unique genes found in the filtered network
def get_filtered_gold(gold_network_path, rank_list, output_path, flag):
    f_gold = open(gold_network_path)
    gold_reader = list(csv.reader(f_gold))
    has_score = True
    if len(gold_reader[0]) < 3:
        has_score = False
    gold_pair_record = {}
    gold_score_record = {}
    unique_gene_list = []
    tf_target_diversity = {}
    gene_regulation_count = {}

    for single_gold_reader in gold_reader[1:]:
        if (single_gold_reader[0] not in rank_list) or (single_gold_reader[1] not in rank_list):
            continue
        gene_pair = [single_gold_reader[0], single_gold_reader[1]]
        str_gene_pair = single_gold_reader[0] + ',' + single_gold_reader[1]

        if single_gold_reader[0] not in unique_gene_list:
            unique_gene_list.append(single_gold_reader[0])
        if single_gold_reader[1] not in unique_gene_list:
            unique_gene_list.append(single_gold_reader[1])

        if single_gold_reader[0] not in tf_target_diversity:
            tf_target_diversity[single_gold_reader[0]] = set()
        tf_target_diversity[single_gold_reader[0]].add(single_gold_reader[1])

        if single_gold_reader[1] not in gene_regulation_count:
            gene_regulation_count[single_gold_reader[1]] = 0
        gene_regulation_count[single_gold_reader[1]] += 1

        if str_gene_pair in gold_score_record:
            continue
        if has_score:
            gold_score_record[str_gene_pair] = float(single_gold_reader[2])
        else:
            gold_score_record[str_gene_pair] = 999
        if gene_pair[0] not in gold_pair_record:
            gold_pair_record[gene_pair[0]] = [gene_pair[1]]
        else:
            gold_pair_record[gene_pair[0]].append(gene_pair[1])

    print(str(len(gold_pair_record)) + ' TFs and ' + str(
            len(gold_score_record)) + ' edges in gold_network consisted of genes in rank_list.')
    print(str(len(unique_gene_list)) + ' genes are common in rank_list and gold_network.')
    rank_density = len(gold_score_record) / (len(gold_pair_record) * (len(rank_list)))
    gold_density = len(gold_score_record) / (len(gold_pair_record) * (len(unique_gene_list)))
    print('Rank genes density = edges/(TFs*(len(rank_gene)-1))=' + str(rank_density))
    print('Gold genes density = edges/(TFs*len(unique_gene_list))=' + str(gold_density))

    tf_target_diversity_avg = sum(len(targets) for targets in tf_target_diversity.values()) / len(tf_target_diversity)
    print('Average Target Gene Diversity = ' + str(tf_target_diversity_avg))

    gene_regulation_count_values = list(gene_regulation_count.values())
    print('Gene Regulation Count Distribution:')
    print('Max:', max(gene_regulation_count_values))
    print('Min:', min(gene_regulation_count_values))
    print('Average:', sum(gene_regulation_count_values) / len(gene_regulation_count_values))

    if flag:
        f_unique = open(output_path, 'w', newline='\n')
        f_unique_writer = csv.writer(f_unique)
        out_unique = np.array(unique_gene_list).reshape(len(unique_gene_list), 1)
        f_unique_writer.writerows(out_unique)
        f_unique.close()

    return gold_pair_record, gold_score_record, unique_gene_list


# Write the filtered gold standard network to a CSV file.
# Outputs each TF-target pair along with its corresponding score.
# The file will have a header: ['TF', 'Target', 'Score'].
def generate_filtered_gold(gold_pair_record,gold_score_record,output_path):
    f_filtered = open(output_path, 'w', newline='\n')
    f_filtered_writer = csv.writer(f_filtered)
    f_filtered_writer.writerow(['TF', 'Target', 'Score'])

    for tf in gold_pair_record:
        once_output = []
        for target in gold_pair_record[tf]:
            single_output = [tf, target, gold_score_record[tf + ',' + target]]
            once_output.append(single_output)
        f_filtered_writer.writerows(once_output)
    f_filtered.close()


# Generate a list of gene pairs with balanced positive and negative samples.
# For each transcription factor (TF), matches each known positive TF-target interaction with one negative sample
# (a gene that is not its known target), ensuring a balanced dataset for supervised learning tasks.
# If there are not enough negatives for a TF, negatives are borrowed from other TFs if available.
# The resulting CSV file includes columns: ['TF', 'Target', 'Label', 'Score'],
# where Label is 1 for positive samples and 0 for negatives.
def get_gene_pair_list(unique_gene_list, gold_pair_record, gold_score_record, output_file):

    all_tf_negtive_record = {}
    for tf in gold_pair_record:
        all_tf_negtive_record[tf] = []
        for target in unique_gene_list:
            if target in gold_pair_record[tf]:
                continue
            all_tf_negtive_record[tf].append(target)

    rank_negtive_record = {}
    for tf in gold_pair_record:
        num_positive = len(gold_pair_record[tf])
        if num_positive > len(all_tf_negtive_record[tf]):
            rank_negtive_record[tf] = all_tf_negtive_record[tf]
            all_tf_negtive_record[tf] = []
        else:
            rank_negtive_record[tf] = all_tf_negtive_record[tf][:num_positive]
            all_tf_negtive_record[tf] = all_tf_negtive_record[tf][num_positive:]

    f_gpl = open(output_file, 'w', newline='\n')
    f_gpl_writer = csv.writer(f_gpl)
    f_gpl_writer.writerow(['TF', 'Target', 'Label', 'Score'])
    stop_flag=False
    for tf in gold_pair_record:
        once_output = []
        for target in gold_pair_record[tf]:
            single_output = [tf, target, '1', gold_score_record[tf + ',' + target]]
            once_output.append(single_output)
            if len(rank_negtive_record[tf]) == 0:
                find_negtive = False
                for borrow_tf in all_tf_negtive_record:
                    if len(all_tf_negtive_record[borrow_tf]) > 0:
                        find_negtive=True
                        single_output = [borrow_tf, all_tf_negtive_record[borrow_tf][0], 0, 0]
                        del all_tf_negtive_record[borrow_tf][0]
                        break
                if not find_negtive:
                    stop_flag = True
                    break
            else:

                single_output = [tf, rank_negtive_record[tf][0], 0, 0]
                del rank_negtive_record[tf][0]
            once_output.append(single_output)
        if stop_flag:
            f_gpl_writer.writerows(once_output[:-1])
            print('Negtive not enough!')
            break
        f_gpl_writer.writerows(once_output)
    f_gpl.close()


# Identify lowly expressed genes based on their expression across cells.
# A gene is considered lowly expressed if it has non-zero expression in at most 10% of the cells.
# Returns a list of lowly expressed gene names.
def get_low_express_gene(origin_expression_record,num_cells):

    gene_list=[]
    threshold=num_cells//10
    for gene in origin_expression_record:
        num=0
        for expression in origin_expression_record[gene]:
            if expression !=0:
                num+=1
                if num>threshold:
                    break
        if num<=threshold:
            gene_list.append(gene)
    return gene_list


