import numpy as np
from sklearn.model_selection import train_test_split
import random
import csv

# Function to randomly shuffle and split indices into train/validation/test sets
def split_index(all_index):
    random.shuffle(all_index)
    part = len(all_index) // 5
    train_index = all_index[:3 * part]
    val_index = all_index[3 * part:4 * part]
    test_index = all_index[4 * part:]
    return train_index,val_index,test_index

# Path to the directory containing input data files
data_path ='/home/qqwu/paper/LSY/1000genes/mHSC-E-Nonspecific-LogPearson_new_clr/WindowSize441-TF5-Target5-Lag64/FullMatrix_TF'

# Load the matrix (features), labels (0/1), and gene pair information
matrix_data = np.load(data_path + '/matrix.npy')
label_data = np.load(data_path + '/label.npy')
gene_pair = np.load(data_path + '/gene_pair.npy')
num_pairs = len(label_data)
pos_index=[index for index,value in enumerate(label_data) if value==1]
neg_index=[index for index,value in enumerate(label_data) if value==0]

# Total number of gene pairs
pos_count = len(pos_index)
neg_count = len(neg_index)

# Separate positive and negative sample indices based on the label
pos_train_index,pos_val_index,pos_test_index=split_index(pos_index)
neg_train_index,neg_val_index,neg_test_index=split_index(neg_index)

# Combine positive and negative indices for each dataset split
train_index=pos_train_index+neg_train_index
val_index=pos_val_index+neg_val_index
test_index=pos_test_index+neg_test_index

# Write indices to file
with open(data_path+'/train_index.txt','w',newline='') as f_train:
    csv_w=csv.writer(f_train,delimiter='\n')
    csv_w.writerow(train_index)
with open(data_path+'/val_index.txt','w',newline='') as f_val:
    csv_w=csv.writer(f_val,delimiter='\n')
    csv_w.writerow(val_index)
with open(data_path+'/test_index.txt','w',newline='') as f_test:
    csv_w=csv.writer(f_test,delimiter='\n')
    csv_w.writerow(test_index)

