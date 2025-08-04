import os
import time
import numpy as np
import matplotlib.pyplot as plt
import tensorflow as tf
from tensorflow.keras.models import Model
from tensorflow.keras.layers import Input, GRU, Reshape, Conv2D, MaxPooling2D, Dropout, Flatten, Dense, Attention
from tensorflow.keras.optimizers import SGD
from tensorflow.keras.optimizers.schedules import ExponentialDecay
from sklearn import metrics
from tensorflow.keras.regularizers import l2
from tensorflow.keras.layers import LSTM
from tensorflow.keras.layers import Bidirectional, GRU
from tensorflow.keras.layers import SimpleRNN, Reshape
import random
from tensorflow.keras.layers import UpSampling2D
from keras.layers import Concatenate

# Load integer indices from a CSV file
def load_index(file):
    import csv
    with open(file, 'r') as f:
        csv_r = list(csv.reader(f, delimiter='\n'))
    return np.array(csv_r).flatten().astype(int)

# Calculate AUROC metric using sklearn, wrapped for TensorFlow
def auroc(y_true, y_pred):
    return tf.py_function(metrics.roc_auc_score, (y_true, y_pred), tf.double)

# Helper function to compute AUPR using sklearn metrics
def aupr_cal(y_true, y_pred):
    precision, recall, thresholds_PR = metrics.precision_recall_curve(y_true, y_pred)
    return metrics.auc(recall, precision)

# Compute AUPR metric wrapped for TensorFlow
def aupr(y_true, y_pred):
    return tf.py_function(aupr_cal, (y_true, y_pred), tf.double)

batch_size = 1024
num_classes = 2
epochs = 50
#epochs = 25
#epochs = 30

# Define the model filename and directory paths for data and model saving
model_name = 'final_model.h5'
dir = '/home/qqwu/paper/LSY/500genes/mHSC-E-Nonspecific-LogPearson/WindowSize441-TF5-Target5-Lag64/'
data_path = dir + 'FullMatrix_TF'
model_save_dir = dir + 'DGRNS-model-Independent/'

# Create the model saving directory if it doesn't exist
if not os.path.isdir(model_save_dir):
    os.makedirs(model_save_dir)

# Load the input feature matrix and label data from numpy files
matrix_data = np.load(data_path + '/matrix.npy')
label_data = np.load(data_path + '/label.npy')
num_pairs = len(label_data)

# Get the total number of gene pairs (samples) from the label data length
index_data_path='/home/qqwu/paper/LSY/500genes/mHSC-E-Nonspecific-LogPearson/WindowSize441-TF5-Target5-Lag64/FullMatrix_TF'
train_index=load_index(index_data_path+'/train_index.txt')
val_index=load_index(index_data_path+'/val_index.txt')
test_index=load_index(index_data_path+'/test_index.txt')

# Extract data and labels using indices
x_train = matrix_data[train_index]
y_train = label_data[train_index]
x_val = matrix_data[val_index]
y_val = label_data[val_index]
x_test = matrix_data[test_index]
y_test = label_data[test_index]

start = time.perf_counter()

# [Modified Section]
############################################
# BiGRU + CNN Architecture
input = Input(shape=x_train.shape[1:], name='input')

# rnn structure
# with tf.name_scope('rnn_section'):
#     input = Reshape((x_train.shape[1], x_train.shape[2]))(input)
#     rnn_output = GRU(128, dropout=0.3, return_sequences=True)(input)
#     rnn_output = LSTM(128, dropout=0.3, return_sequences=True)(input)
#     rnn_output = SimpleRNN(128, dropout=0.3, return_sequences=True)(input)
#     rnn_output = Bidirectional(GRU(64, dropout=0.3, return_sequences=False))(input)
#     output = Dense(1, activation='sigmoid')(rnn_output)

# bigru
with tf.name_scope('bigru_section'):
    gru_input = Reshape((x_train.shape[1], x_train.shape[2]))(input)
    rnn_output = Bidirectional(GRU(64, dropout=0.3, return_sequences=True))(gru_input)

# cnn structure
# with tf.name_scope('cnn_section'):
#     cnn_input = input
#     cnn_input = Reshape((x_train.shape[1], 128, 1))(rnn_output)
#     layer1_output = Conv2D(64, (3, 3), padding='same', activation='relu')(cnn_input)
#     layer2_output = Dropout(0.3)(MaxPooling2D(pool_size=(2, 2))(Conv2D(32, (3, 3), activation='relu')(layer1_output)))
#     layer3_output = Conv2D(64, (3, 3), padding='same', activation='relu')(layer2_output)
#     layer4_output = Dropout(0.3)(MaxPooling2D(pool_size=(2, 2))(Conv2D(32, (3, 3), activation='relu')(layer3_output)))
#     cnn_output = Dropout(0.5)(Dense(512, activation='relu')(Flatten()(layer4_output)))
#     output = Dense(1, activation='sigmoid')(Dense(128, activation='relu')(cnn_output))

# cnn
with tf.name_scope('cnn_section'):
    cnn_input = Reshape((x_train.shape[1], 128, 1))(rnn_output)
    layer1_output = Conv2D(64, (3, 3), padding='same', activation='relu')(cnn_input)
    layer2_output = Dropout(0.3)(MaxPooling2D(pool_size=(2, 2))(Conv2D(32, (3, 3), activation='relu')(layer1_output)))
    cnn_output = Dropout(0.5)(Dense(512, activation='relu')(Flatten()(layer2_output)))
    output = Dense(1, activation='sigmoid')(Dense(128, activation='relu')(cnn_output))

# Set initial learning rate and define exponential decay schedule for learning rate
initial_learning_rate = 0.01
lr_schedule = ExponentialDecay(initial_learning_rate, decay_steps=100000, decay_rate=0.96, staircase=True)
############################################

# Create the model using previously defined input and output layers
model = Model(inputs=input, outputs=output)
sgd = SGD(learning_rate=lr_schedule, momentum=0.9, nesterov=True)

# Compile the model with binary crossentropy loss and accuracy metric
model.compile(optimizer=sgd, loss='binary_crossentropy', metrics=['accuracy'])

# Print model architecture summary
model.summary()

# Train the model on training data with validation data, specified batch size and epochs, shuffling each epoch
history = model.fit(x_train, y_train, validation_data=(x_val, y_val), batch_size=batch_size, epochs=epochs, shuffle=True)

# Save the trained model to the specified directory
model_path = os.path.join(model_save_dir, model_name)
model.save(model_path)
print('Saved trained model at %s ' % model_path)

# Predict labels on test set using the trained model and save test labels and predictions for later analysis
y_predict = model.predict(x_test)
np.save(model_save_dir + 'y_test.npy', y_test)
np.save(model_save_dir + 'y_predict.npy', y_predict)

# Record and print total running time of the process
end = time.perf_counter()
print('Running time:'+str(end-start))

# Plot training and validation accuracy and loss curves
plt.figure(figsize=(10, 6))
plt.subplot(1, 2, 1)
plt.plot(history.history['accuracy'])
plt.plot(history.history['val_accuracy'])
plt.title('model accuracy')
plt.ylabel('accuracy')
plt.xlabel('epoch')
plt.grid()
plt.legend(['train', 'val'], loc='upper left')
plt.subplot(1, 2, 2)
plt.plot(history.history['loss'])
plt.plot(history.history['val_loss'])
plt.title('model loss')
plt.ylabel('loss')
plt.xlabel('epoch')
plt.legend(['train', 'val'], loc='upper left')
plt.grid()
plt.savefig(model_save_dir + 'train_process.pdf')

# Plot ROC curve for test predictions
plt.figure(figsize=(10, 6))
fpr, tpr, thresholds = metrics.roc_curve(y_test, y_predict, pos_label=1)
auc = metrics.auc(fpr, tpr)
plt.plot(fpr, tpr, label='AUROC = %0.5f)' % auc)
plt.grid()
plt.plot([0, 1], [0, 1])
plt.title('ROC curve')
plt.xlabel('FP')
plt.ylabel('TP')
plt.ylim([0, 1])
plt.xlim([0, 1])
plt.legend(loc="lower right")
print('AUC in figure:', auc)
plt.savefig(model_save_dir + 'AUCofTest.pdf')

# Calculate and plot Precision-Recall (PR) curve for test predictions
precision, recall, thresholds_PR = metrics.precision_recall_curve(y_test, y_predict)
AUPR = metrics.auc(recall, precision)
plt.figure(1)
plt.plot(precision, recall, label='Area Under the Curve (AUPR = %0.5f)' % AUPR)
plt.title('PR curve')
plt.xlabel("Recall")
plt.ylabel("Precision")
plt.legend(loc="lower right")
plt.ylim([0, 1])
plt.xlim([0, 1])
plt.savefig(model_save_dir + 'AUPRofTest.pdf')

print("AUPR in figure:", AUPR)
