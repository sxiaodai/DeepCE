import matplotlib.pyplot as plt
from matplotlib import font_manager
font = font_manager.FontProperties(fname='C:/Users/86188/Desktop/simsun.ttf')
plt.rc('font', family=['Songti SC','Times New Roman'])


## Line plot of AUROC and AUPR changes with different mean aggregation numbers k.
# import matplotlib.pyplot as plt
#
# aggregation_numbers = [5, 10, 15, 20]
# auroc_values = [0.706983333, 0.714283333, 0.685166667, 0.640666667]
# aupr_values = [0.180816667, 0.188316667, 0.147, 0.129266667]
#
# color_auroc = '#66c2a5'
# color_aupr = '#fc8d62'
#
# plt.figure(figsize=(12, 6))
#
# plt.subplot(2, 1, 1)
# plt.plot(aggregation_numbers, auroc_values, label='AUROC', marker='o', color=color_auroc)
# # plt.xlabel('Aggeration Number (k)', color='black', fontsize=12)
# plt.xticks(aggregation_numbers, color='black', fontsize=18)
# plt.yticks(color='black', fontsize=18)
#
# for i in range(len(aggregation_numbers)):
#     plt.text(aggregation_numbers[i]-0.3, auroc_values[i]-0.003, f'{auroc_values[i]:.4f}', color='black', fontsize=18)
#
# plt.legend(fontsize=18)
# plt.grid(True)
# # plt.title('AUROC')
#
# plt.subplot(2, 1, 2)  # (行数, 列数, 当前子图的位置)
# plt.plot(aggregation_numbers, aupr_values, label='AUPR', marker='s', color=color_aupr)
# plt.xlabel('Aggeration Number (k)', color='black', fontsize=18)
# plt.xticks(aggregation_numbers, color='black', fontsize=18)
# plt.yticks(color='black', fontsize=18)
#
# for i in range(len(aggregation_numbers)):
#     plt.text(aggregation_numbers[i]-0.3, aupr_values[i]-0.003, f'{aupr_values[i]:.4f}', color='black',  fontsize=18)
#
# plt.legend(fontsize=18)
# plt.grid(True)
# # plt.title('AUPR')
#
# plt.tight_layout()
# plt.show()


## Line plot of AUROC and AUPR changes with different data processing methods
# import matplotlib.pyplot as plt
#
# aggregation_numbers = ['Baseline', 'Aggregation', 'Adjustment', 'Combined']
# auroc_values = [0.610116667, 0.714283333, 0.699283333, 0.817483333]
# aupr_values = [0.133833333, 0.188316667, 0.180566667, 0.362383333]
#
# color_auroc = '#66c2a5'
# color_aupr = '#fc8d62'
#
# plt.figure(figsize=(12, 6))
#
# plt.subplot(2, 1, 1)
# plt.plot(aggregation_numbers, auroc_values, label='AUROC', marker='o', color=color_auroc)
# # plt.xlabel('Data Processing Methods', color='black', fontsize=12)
# plt.xticks(aggregation_numbers,  color='black', fontsize=18)
# plt.yticks(color='black', fontsize=18)
#
# for i in range(len(aggregation_numbers)):
#     plt.text(i-0.08, auroc_values[i] - 0.009, f'{auroc_values[i]:.4f}', color='black', fontsize=18)
#
# plt.legend(fontsize=18)
# plt.grid(True)
# # plt.title('AUROC')
#
# plt.subplot(2, 1, 2)
# plt.plot(aggregation_numbers, aupr_values, label='AUPR', marker='s', color=color_aupr)
# plt.xlabel('Method', color='black', fontsize=18)
# plt.xticks(aggregation_numbers, color='black', fontsize=18)
# plt.yticks(color='black', fontsize=18)
#
# for i in range(len(aggregation_numbers)):
#     plt.text(i-0.08, aupr_values[i] - 0.009, f'{aupr_values[i]:.4f}', color='black', fontsize=18)
#
# plt.legend(fontsize=18)
# plt.grid(True)
# # plt.title('AUPR')
#
# plt.tight_layout()
# plt.show()


## Comparison of model performance across different methods
# import matplotlib.pyplot as plt
# import numpy as np
#
#
# data = {
#     'mHSC-E(500)': {
#         'DeepSEM': {'Auroc': [0.5145], 'Aupr': [0.0235]},
#         'GENELink': {'Auroc': [0.6591], 'Aupr': [0.0527]},
#         'DGRNS': {'Auroc': [0.7624], 'Aupr': [0.2644]},
#         'DeepCE': {'Auroc': [0.8559], 'Aupr': [0.5281]},
#     },
#
#     'mHSC-E(1000)':{
#         'DeepSEM': {'Auroc': [0.5155], 'Aupr': [0.0205]},
#         'GENELink': {'Auroc': [0.7604], 'Aupr': [0.1653]},
#         'DGRNS': {'Auroc': [0.7354], 'Aupr': [0.2257]},
#         'DeepCE': {'Auroc': [0.8574], 'Aupr': [0.3891]},
#     },
#
#     'mHSC-L(500)': {
#         'DeepSEM': {'Auroc': [0.4989], 'Aupr': [0.0290]},
#         'GENELink': {'Auroc': [0.6271], 'Aupr': [0.0847]},
#         'DGRNS': {'Auroc': [0.6804], 'Aupr': [0.1189]},
#         'DeepCE': {'Auroc': [0.814], 'Aupr': [0.3358]},
#     },
#
#     'mHSC-L(1000)': {
#         'DeepSEM': {'Auroc': [0.5669], 'Aupr': [0.0641]},
#         'GENELink': {'Auroc': [0.5802], 'Aupr': [0.1158]},
#         'DGRNS': {'Auroc': [0.6162], 'Aupr': [0.1121]},
#         'DeepCE': {'Auroc': [0.7734], 'Aupr': [0.2656]},
#     },
#
#     'mHSC-GM(500)': {
#         'DeepSEM': {'Auroc': [0.5194], 'Aupr': [0.0324]},
#         'GENELink': {'Auroc': [0.7512], 'Aupr': [0.2890]},
#         'DGRNS': {'Auroc': [0.6612], 'Aupr': [0.1380]},
#         'DeepCE': {'Auroc': [0.7751], 'Aupr': [0.3976]},
#     },
#
#     'mHSC-GM(1000)': {
#         'DeepSEM': {'Auroc': [0.5046], 'Aupr': [0.0293]},
#         'GENELink': {'Auroc': [0.7788], 'Aupr': [0.3295]},
#         'DGRNS': {'Auroc': [0.6323], 'Aupr': [0.1117]},
#         'DeepCE': {'Auroc': [0.8291], 'Aupr': [0.2581]},
#     },
#
#     'hESC(500)': {
#         'DeepSEM': {'Auroc': [0.5070], 'Aupr': [0.0164]},
#         'GENELink': {'Auroc': [0.7015], 'Aupr': [0.0473]},
#         'DGRNS': {'Auroc': [0.7517], 'Aupr': [0.1871]},
#         'DeepCE': {'Auroc': [0.8256], 'Aupr': [0.4008]},
#     },
#
#     'hESC(1000)': {
#         'DeepSEM': {'Auroc': [0.5058], 'Aupr': [0.0140]},
#         'GENELink': {'Auroc': [0.7024], 'Aupr': [0.0398]},
#         'DGRNS': {'Auroc': [0.7342], 'Aupr': [0.1294]},
#         'DeepCE': {'Auroc': [0.7839], 'Aupr': [0.1807]},
#     },
#
#     'hHep(500)': {
#         'DeepSEM': {'Auroc': [0.5020], 'Aupr': [0.0155]},
#         'GENELink': {'Auroc': [0.7138], 'Aupr': [0.0366]},
#         'DGRNS': {'Auroc': [0.6322], 'Aupr': [0.2075]},
#         'DeepCE': {'Auroc': [0.7774], 'Aupr': [0.2784]},
#     },
#
#     'hHep(1000)': {
#         'DeepSEM': {'Auroc': [0.5011], 'Aupr': [0.0132]},
#         'GENELink': {'Auroc': [0.711], 'Aupr': [0.0309]},
#         'DGRNS': {'Auroc': [0.671], 'Aupr': [0.0876]},
#         'DeepCE': {'Auroc': [0.7361], 'Aupr': [0.1437]},
#     },
# }
#
# # Set up data for the plot
# # methods = list(data['mHSC-GM(500)'].keys())
# # categories = ['mHSC-E (1000) AUROC', 'mHSC-E (1000) AUPR', 'mHSC-L (1000) AUROC', 'mHSC-L (1000) AUPR', 'mHSC-GM (1000) AUROC', 'mHSC-GM (1000) AUPR']
# # values = [
# #     [data['mHSC-E(1000)'][method]['Auroc'][0] for method in methods],
# #     [data['mHSC-E(1000)'][method]['Aupr'][0] for method in methods],
# #     [data['mHSC-L(1000)'][method]['Auroc'][0] for method in methods],
# #     [data['mHSC-L(1000)'][method]['Aupr'][0] for method in methods],
# #     [data['mHSC-GM(1000)'][method]['Auroc'][0] for method in methods],
# #     [data['mHSC-GM(1000)'][method]['Aupr'][0] for method in methods],
# # ]
#
# methods = list(data['mHSC-GM(500)'].keys())
# categories = ['hHep (500) AUROC', 'hHep (500) AUPR', 'hESC (500) AUROC', 'hESC (500) AUPR']
# values = [
#     [data['hHep(500)'][method]['Auroc'][0] for method in methods],
#     [data['hHep(500)'][method]['Aupr'][0] for method in methods],
#     [data['hESC(500)'][method]['Auroc'][0] for method in methods],
#     [data['hESC(500)'][method]['Aupr'][0] for method in methods],
# ]
#
# # Set up x-axis positions and bar width
# x = np.arange(len(categories))
# width = 0.15
#
# fig, ax = plt.subplots(figsize=(12, 5))
# # fig, ax = plt.subplots(figsize=(8, 5))
#
# # fig, ax = plt.subplots(figsize=(10, 6))
#
# for i, (method, color) in enumerate(zip(methods, plt.cm.Set2.colors)):
#     bars = ax.bar(x + i * width - width * (len(methods) - 1) / 2,
#                   [values[j][i] for j in range(len(values))],
#                   width, label=method, color=color)
#
#     # Add the value labels on top of each bar
#     for bar in bars:
#         yval = bar.get_height()
#         ax.text(bar.get_x() + bar.get_width() / 2, yval+0.005, f'{yval:.4f}', ha='center', va='bottom',rotation=90,fontsize=12)
#
# for i, (method, color) in enumerate(zip(methods, plt.cm.Set2.colors)):
#     ax.bar(x + i * width - width * (len(methods) - 1) / 2,
#            [values[j][i] for j in range(len(values))],
#            width, label=method, color=color)
#
# # Set labels and title
# # ax.set_xlabel('Dataset and Metric')
# # ax.set_ylabel('Value')
# # ax.set_title('Performance Comparison of Methods on Different Datasets')
#
# ax.set_xticks(x)
# ax.set_xticklabels(categories, rotation=0,fontsize=12)
# ax.legend(methods,fontsize=12)
# ax.axvline(x=1.5, color='grey', linewidth=1.5, linestyle='--')
# # ax.axvline(x=3.5, color='grey', linewidth=1.5, linestyle='--')
# ax.set_ylim(0, 1)
# ax.tick_params(axis='y', labelsize=12)
# ax.tick_params(axis='x', labelsize=12)
# # Show the plot
# plt.tight_layout()
# plt.show()
#


# import matplotlib.pyplot as plt
# import matplotlib.image as mpimg
#
# img_paths = [
#     r"C:\Users\86188\Desktop\plot3.png",
#     r"C:\Users\86188\Desktop\plot4.png",
# ]
#
# fig, axes = plt.subplots(2, 1, figsize=(12, 8))
#
# for ax, img_path in zip(axes, img_paths):
#     img = mpimg.imread(img_path)
#     ax.imshow(img)
#     ax.axis('off')  #
#
# plt.tight_layout()
# plt.show()
