# This file uses code from the DGRNS project.
# Original repository: https://github.com/guofei-tju/DGRNS
# Some modifications and adaptations have been made for this work.
from __future__ import print_function
import keras
from keras.preprocessing.image import ImageDataGenerator
from keras.models import Sequential
from keras.models import load_model
from keras.layers import Dense, Dropout, Activation, Flatten
from keras.layers import Conv2D, MaxPooling2D
from keras.optimizers import SGD
from keras.callbacks import EarlyStopping,ModelCheckpoint
import openpyxl
from keras import backend as K
import os,sys
import numpy as np
import matplotlib.pyplot as plt
from sklearn import metrics
import itertools
import csv

# Load gene expression data from CSV and filter to keep only specified genes
def get_filtered_expression_data(gene_expression_path,gene_name_list):
    f_expression = open(gene_expression_path)
    expression_reader = list(csv.reader(f_expression))

    expression_record = {}
    Num_genes = 0
    for single_expression_reader in expression_reader[1:]:
        if single_expression_reader[0] in gene_name_list:
            if single_expression_reader[0] not in expression_record:
                expression_record[single_expression_reader[0]] = np.array(list(map(float, single_expression_reader[1:])))
    return expression_record

# Load matrix, label, and gene pair data from specified directory and return as arrays
def load_data(data_path):
    matrix_array=np.load(data_path+'/matrix.npy')
    label_array=np.load(data_path+'/label.npy')
    pair_array=np.load(data_path+'/gene_pair.npy')
    print(matrix_array.shape)
    return(matrix_array,label_array,pair_array)

# Function to load test data and a trained model, then predict and save results
def predict(matrix_path,model_path,save_dir,fold_index):
    (x_test, y_test,pair_test) = load_data(matrix_path)
    print(x_test.shape, 'x_test samples')
    model = load_model(model_path)
    if not os.path.isdir(save_dir):
        os.makedirs(save_dir)
    print ('load model and predict')
    y_predict = model.predict(x_test)
    f_record.write(model_path)
    f_record.write('\n')
    output_predict = pair_test.copy()
    file = open(save_dir + '/fold'+str(fold_index)+'_y_predict.txt', 'w')
    output=''
    for index_y_predict in range(len(y_predict)):
        tf,target=output_predict[index_y_predict].split(',')
        output=output+tf+'\t'+target+'\t'+str(y_predict[index_y_predict][0])+'\n'
    file.write(output)
    file.close()
    return y_test,y_predict,output_predict

# Convert continuous data to binary classification labels based on a given threshold
def continuousToClassification(sourcedata,threshold):
    result=sourcedata.copy()
    for i in range(len(sourcedata)):
        if sourcedata[i]<threshold:
            result[i]=0
        else:
            result[i] = 1
    return result

# Estimate model performance by predicting and evaluating on test data
def estimation(matrix_path, model_path, save_dir, f_record,fold_index):
    # Get true labels, predicted scores, and gene pairs from prediction function
    y_test, y_predict,output_predict = predict(matrix_path, model_path, save_dir,fold_index)
    # Convert predictions and true labels to lists for further processing
    y_predict=y_predict[:,0].tolist()
    y_test=y_test.tolist()
    # Get indices that would sort predictions in descending order
    map_index=np.argsort(-np.array(y_predict))
    # Initialize binary prediction array with zeros
    binary_y_predict=np.zeros(len(y_predict))
    # Count the number of positive samples in the test set
    true_positive=0
    for i in y_test:
        if i == 1:
            true_positive+=1
    # Calculate the density (proportion of positive edges)
    density = true_positive / len(y_test)
    print(str(len(y_test)) + ' edges in test set.')
    print(str(true_positive)+' positive edges in test set.')
    print('Density=' + str(density))
    # Assign top-ranked predictions as positive based on true positive count
    for i in range(true_positive):
        binary_y_predict[map_index[i]]=1
    y_predict_classification=binary_y_predict

    # [Added Section]
    ############################################
    # Identify false positive pairs where true label is 0 but predicted as 1
    false_positive_pairs = []
    for i in range(len(y_test)):
        if y_test[i] == 0 and y_predict_classification[i] == 1:
            tf, target = output_predict[i].split(',')
            false_positive_pairs.append(f'{tf}\t{target}\t1')
    # Save these potential false positive pairs to a file
    with open(save_dir + '/fold' + str(fold_index) + '_possible_positive.txt', 'w') as f_fp:
        for line in false_positive_pairs:
            f_fp.write(line + '\n')

    # Collect all predicted positive pairs (regardless of true label)
    false_positive_pairs = []
    for i in range(len(y_test)):
        if y_predict_classification[i] == 1:
            tf, target = output_predict[i].split(',')
            false_positive_pairs.append(f'{tf}\t{target}\t1')
    # Save all predicted positives to another file
    with open(save_dir + '/fold' + str(fold_index) + 'predict_positive.txt', 'w') as f_fp:
        for line in false_positive_pairs:
            f_fp.write(line + '\n')
    ############################################

    # Print and write evaluation metrics header
    table_head = 'EPR\tAUPRC\taccuracy\tprecision\trecall\tf1-score\tAUC\tAUPR'
    print(table_head)
    f_record.write(table_head + '\n')
    # Convert true labels to numpy array for metric calculations
    array_y_test=np.array(y_test)
    # Calculate common classification metrics
    accuracy= metrics.accuracy_score(array_y_test, y_predict_classification)
    precision= metrics.precision_score(array_y_test, y_predict_classification,pos_label=1)
    recall = metrics.recall_score(array_y_test, y_predict_classification,pos_label=1)
    f1 = metrics.f1_score(array_y_test, y_predict_classification)

    # ROC curve and AUC
    fpr, tpr, _ = metrics.roc_curve(y_test, y_predict)
    auc = metrics.auc(fpr, tpr)

    # Precision-Recall curve and AUPR
    list_precision, list_recall, _ = metrics.precision_recall_curve(y_test, y_predict)
    aupr = metrics.auc(list_recall, list_precision)

    # Calculate enrichment metrics by normalizing precision and AUPR with density
    single_class_result = str(round(precision / density, 6)) + '\t' + str(round(aupr / density, 6)) + '\t' + str(
        round(accuracy, 6)) + '\t' + str(
        round(precision, 6)) + '\t' + str(round(recall, 6)) + '\t' + str(round(f1, 6)) + '\t' + str(
        round(auc, 6)) + '\t' + str(round(aupr, 6))
    print(single_class_result)
    f_record.write(single_class_result + '\n')

    # Plot ROC curve
    plt.figure()
    colors = itertools.cycle(['aqua', 'darkorange', 'cornflowerblue'])
    plt.plot(fpr, tpr, label='ROC curve(area = {})'.format(auc))
    plt.plot([0, 1], [0, 1], 'k--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC curve')
    plt.legend(loc="lower right")
    plt.savefig(save_dir + '/fold'+str(fold_index)+'-ROC.pdf', bbox_inches='tight')
    # plt.show()
    plt.close()

    # Plot Precision-Recall curve
    plt.figure()
    plt.plot(list_precision, list_recall, label='PR curve(area = {})'.format(aupr))
    plt.plot([0, 1], [1, 0], 'k--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title('PR curve')
    plt.legend(loc="lower right")
    plt.savefig(save_dir + '/fold'+str(fold_index)+ '-PR.pdf', bbox_inches='tight')
    # plt.show()
    plt.close()

# This script performs model evaluation using a saved deep learning model.
# It loads test data from `matrix_path`, evaluates model predictions on the data,
# and saves classification results, evaluation metrics, and plots (ROC, PR curves) to `save_dir`.
matrix_path = '/home/qqwu/paper/LSY/500genes/hHep-Nonspecific-LogPearson/WindowSize173-TF2-Target2-Lag64/FullMatrix_TF'
save_dir ='/home/qqwu/paper/LSY/500genes/hHep-Nonspecific-LogPearson/WindowSize173-TF2-Target2-Lag64/model/FullResult_TF'
if not os.path.isdir(save_dir):
    os.makedirs(save_dir)
f_record = open(save_dir + '/record.txt', 'w')
for index_fold in range(0,1):
    model_path = '/home/qqwu/paper/LSY/500genes/hHep-Nonspecific-LogPearson_new_clr/WindowSize173-TF2-Target2-Lag64/DGRNS-model-Independent/final_model.h5'
    estimation(matrix_path, model_path, save_dir, f_record, index_fold)
print('done')
f_record.close()
