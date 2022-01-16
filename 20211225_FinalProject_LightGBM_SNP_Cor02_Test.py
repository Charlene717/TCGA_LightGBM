# -*- coding: utf-8 -*-
"""
Created on Sat Oct  2 16:32:21 2021

@author: user
"""
## 設定路徑

Path = 'D:/Dropbox/1-成大醫工/博士班/1-課程/機器學習/Final Projects/Data_LUAD/'
Path2 = 'D:/Dropbox/1-成大醫工/博士班/1-課程/機器學習/Final Projects/Data_LUAD/Output20211226_SNP_Cor02_V3/'
RunCondition = 'LightGBMTar'

#%%
## 匯入檔案
# data=pd.read_csv('C:/workspace/00690.csv')
"""
import OS
FolderPath = 'D:/Ch_Bioinformatics Dropbox/Chang Charlene/1-成大醫工/博士班/1-課程/機器學習/HW1/hw1_dataset/' 
allFileList = OS.listdir(FolderPath)
"""
import math
import pandas as pd

Train_data_list = [
     Path + 'TCGA_Com_TP53_SNP_cor02_Train.csv'
#     Path + 'data1/train.csv',
#     Path + 'data2/train.csv',
#     Path + 'data3/train.csv', 
#     Path + 'data4/train.csv',
#     Path + 'data5/train.csv',
#     Path + 'data6/train.csv',
#     Path + 'data7/train.csv'

    ]


Test_data_list = [
     Path + 'TCGA_Com_TP53_SNP_cor02_Test.csv'    
#     Path + 'data1/test.csv',
#     Path + 'data2/test.csv',
#     Path + 'data3/test.csv',
#     Path + 'data4/test.csv',
#     Path + 'data5/test.csv',
#     Path + 'data6/test.csv',
#     Path + 'data7/test.csv'

    ]

TestAns_data_list = [
     Path + 'TCGA_Com_TP53_SNP_cor02_TestAns.csv'    
#     Path + 'data1/test.csv',
#     Path + 'data2/test.csv',
#     Path + 'data3/test.csv',
#     Path + 'data4/test.csv',
#     Path + 'data5/test.csv',
#     Path + 'data6/test.csv',
#     Path + 'data7/test.csv'

    ]



Train_dataset = []
for filepath_Train in Train_data_list:
    Train_dataset.append(pd.read_csv(filepath_Train, header=0))
del filepath_Train

Test_dataset = []
for filepath_Test in Test_data_list:
    Test_dataset.append(pd.read_csv(filepath_Test, header=0))
del filepath_Test  
  
    
TestAns_dataset = []
for filepath_TestAns in TestAns_data_list:
    TestAns_dataset.append(pd.read_csv(filepath_TestAns, header=0))
del filepath_TestAns  


#%%
## Load Package
# Grid Search
from sklearn.model_selection import GridSearchCV
from sklearn import metrics
import lightgbm as lgb
import numpy

## 參數設定
config = {#'trials': 10,
      #'train_ratio': 0.7,
      #'val_ratio': 0.10,
      'cv': 10,
      'test_ratio': 0.1,
      'n_estimators': 1000,
      'bagging_fraction': 0.8,
      'feature_fraction':0.8#,
      #'reg_lambda':0.0001
}


#%%
## Hyperparameter for Target
# https://www.itread01.com/content/1543214717.html
boosting_Target = ['rf'] #['gbdt','dart','goss','rf']
learning_rate_Target = list(range(1, 5))
learning_rate_Target_new_list = [i * 0.01 for i in learning_rate_Target]
max_depth_Target = list(range(3, 7, 2))
#num_leaves_Target = list(range(math.ceil(len(cdata_feature_array_copy)/10), math.ceil(len(cdata_feature_array_copy)/8)))
num_leaves_Target= list(range(10, 50, 10))
min_data_in_leaf_Target = list(range(20, 70, 10))


reg_lambda = list(range(0, 5))
reg_lambda_new_list = [i * 0.0001 for i in reg_lambda]

#reg_lambda_new_list = list(range(0, 5, 1))

# create a parameter grid: map the parameter names to the values that should be searched
param_grid_Target = dict(learning_rate = learning_rate_Target_new_list , 
                         max_depth = max_depth_Target, num_leaves = num_leaves_Target,
                         min_data_in_leaf = min_data_in_leaf_Target, reg_lambda = reg_lambda_new_list)
print(param_grid_Target)

#%%
##　資料預處理 & 切分資料
cdata_Train_result = []
cdata_Test_result = []

cdata_Train_MAE = []
cdata_Test_MAE = []

cdata_Test_RMSD = []
cdata_Train_RMSD = []

cdata_Test_R2 = []
cdata_Train_R2 = []

cdata_Test_score = []
cdata_Train_score = []

Best_Hyperparameter_Target_All = []

cdata_train = Train_dataset[0]
keys = cdata_train.keys()
cdata_feature_array = cdata_train[keys[:-1]]
cdata_target_array = cdata_train[keys[-1]].to_numpy()
cdata_target_array = cdata_target_array.astype('int')

# Split data into train and test subsets (train:0.75, test:0.25)
from sklearn.model_selection import train_test_split
train_data, test_data, train_target, test_target = train_test_split(
    cdata_feature_array, cdata_target_array, test_size = config['test_ratio'], shuffle = False)

print('Size of training dataset:', train_data.shape)
print('Size of testing dataset:', test_data.shape)



## 訓練模型: LightGBM
LGBM_Ori = lgb.LGBMRegressor(application='multiclass')

# instantiate the grid
grid = GridSearchCV(LGBM_Ori, param_grid_Target, cv=config['cv'], scoring='neg_root_mean_squared_error',return_train_score=True)

del LGBM_Ori
# fit the grid with data
grid.fit(train_data, train_target)
    
# view the results
SumAccuracy = pd.DataFrame(grid.cv_results_)[['mean_test_score', 'std_test_score','mean_train_score' ,'std_train_score',  'params']]
SumAccuracy = SumAccuracy.sort_values(['mean_test_score','std_test_score'],ascending=False)

TryP_Accuracy_df = pd.DataFrame.from_dict(SumAccuracy)
TryP_Accuracy_df.to_csv(Path2 + RunCondition + '_TryP_Accuracy_data' + str(1) + '.csv' ,header=True , index=False )

# SumAccuracy['params'][1]['boosting']
# examine the best model
print(grid.best_score_)
print(grid.best_params_)

# make predictions of grid search model
predicted = grid.predict(test_data)
accuracy_grid = metrics.mean_absolute_error(test_target, predicted)
print('\nAccuracy of grid search model: ', accuracy_grid)

# Using the best parameters in LGBM model
#LGBM = lgb.LGBMRegressor(application='multiclass', boosting='gbdt', learning_rate=0.1, max_depth=-5, feature_fraction=0.5, random_state=42)
LGBM = lgb.LGBMRegressor(application='multiclass', learning_rate=grid.best_params_['learning_rate'], max_depth=grid.best_params_['max_depth'], 
                         num_leaves = grid.best_params_['num_leaves'] ,min_data_in_leaf = grid.best_params_['min_data_in_leaf'], 
                         reg_lambda = grid.best_params_['reg_lambda'], colsample_bytree = None, subsample =None, min_child_samples = None, 
                         n_estimators = config['n_estimators'],
                         bagging_fraction = config['bagging_fraction'],
                         feature_fraction = config['feature_fraction'])

# Best Hyperparameter in Target
Best_Hyperparameter_Target = [grid.best_params_['learning_rate'],grid.best_params_['max_depth'],
                           grid.best_params_['num_leaves'],grid.best_params_['min_data_in_leaf'],
                           grid.best_params_['reg_lambda']]

LGBM.fit(train_data, train_target)

# make predictions of LGBM model
predicted_test = LGBM.predict(test_data)
predicted_train = LGBM.predict(train_data)

# Calculate accuracy
score_test = LGBM.score(test_data, test_target)
score_train = LGBM.score(train_data, train_target)
MAE_test = metrics.mean_absolute_error(test_target, predicted_test)
MAE_train = metrics.mean_absolute_error(train_target, predicted_train)
print('MAE of LGBM model: ', MAE_test)
RMSD_test = metrics.mean_squared_error(test_target, predicted_test)
RMSD_train = metrics.mean_squared_error(train_target, predicted_train)
print('RMSD of LGBM model: ', RMSD_test)

R2_test = metrics.r2_score(test_target, predicted_test)
R2_train = metrics.r2_score(train_target, predicted_train)
print('R2 of LGBM model: ', R2_test)


## Calculate accuracy
#accuracy = metrics.accuracy_score(test_target, predicted)
#print('Accuracy of testing dataset: ', accuracy)
## np.r_[cdata_feature_array,predicted]
import numpy as np
Output = np.c_[cdata_feature_array,cdata_target_array]
cdata_Train_result.append(Output)
cdata_Test_MAE.append(MAE_test)
cdata_Train_MAE.append(MAE_train)
cdata_Test_RMSD.append(RMSD_test)
cdata_Train_RMSD.append(RMSD_train)
cdata_Test_R2.append(R2_test)
cdata_Train_R2.append(R2_train)
cdata_Test_score.append(score_test)
cdata_Train_score.append(score_train)

## Record the best Hyperparameter
Best_Hyperparameter_Target_All.append(Best_Hyperparameter_Target)
#Best_Hyperparameter_NaN_All.append(Best_Hyperparameter_NaN)

train_df = pd.DataFrame.from_dict(Output)
train_df.to_csv(Path2 + RunCondition + '_Train_data' + str(1) + '.csv' ,header=False , index=False )

del Output,MAE_test,MAE_train


## 預測驗證集資料 
    
# 套用訓練模型
cdata_test = Test_dataset[0]
keys2 = cdata_test.keys()
"""
"""
cdata_feature_array_Target = cdata_test[keys2[:-1]]
predicted2 =  LGBM.predict(cdata_feature_array_Target)

# https://codertw.com/%E7%A8%8B%E5%BC%8F%E8%AA%9E%E8%A8%80/357982/
Output = np.c_[cdata_feature_array_Target,predicted2]
cdata_Test_result.append(Output)

## 匯出CSV檔
#　https://www.dotblogs.com.tw/CYLcode/2020/03/17/175255                
#def OutputCSV(): 
df_SAMPLE2 = pd.DataFrame.from_dict(Output)
df_SAMPLE2.to_csv(Path2 + RunCondition + '_Output_data' + str(1) + '.csv' ,header=False , index=False )
'''    
MAE_train_df = pd.DataFrame.from_dict(cdata_Train_MAE)
MAE_train_df.to_csv(Path2 + RunCondition + '_TrainAccuracyFin_data.csv' ,header=False , index=False )
MAE_test_df = pd.DataFrame.from_dict(cdata_Test_MAE)
MAE_test_df.to_csv(Path2 + RunCondition + '_TestAccuracyFin_data.csv' ,header=False , index=False )
'''
MAE_train_df = pd.DataFrame.from_dict(cdata_Train_MAE)
MAE_test_df = pd.DataFrame.from_dict(cdata_Test_MAE)
score_train_df = pd.DataFrame.from_dict(cdata_Train_score)
score_test_df = pd.DataFrame.from_dict(cdata_Train_score)


Accuracy_df = pd.DataFrame(list(zip(cdata_Train_MAE, cdata_Test_MAE,cdata_Train_RMSD,cdata_Test_RMSD,cdata_Train_R2,cdata_Test_R2,cdata_Train_score,cdata_Test_score)),
           columns =['Train_MAE', 'Val_MAE','Train_RMSD', 'Val_RMSD','Train_R2', 'Val_R2', 'Train_Score', 'Val_Score'], 
           #index =['data1', 'data2', 'data3', 'data4', 'data5', 'data6', 'data7'])    
           index =['data1'])    
Accuracy_df.to_csv(Path2 + RunCondition + '_AccuracyFinal.csv' ,header=True , index=True)

'''
Hyperparameter.df再加上註解
'''

#Best_Hyperparameter_NaN_All = pd.DataFrame(list(Best_Hyperparameter_NaN_All),
#                          columns =['learning_rate', 'max_depth', 'num_leaves', 'min_data_in_leaf'], 
#                          index =['data1', 'data2', 'data3', 'data4', 'data5', 'data6', 'data7'])
#Best_Hyperparameter_NaN_All.to_csv(Path2 + RunCondition + '_Best_Hyperparameter_NaN.csv' ,header=True , index=True)

Best_Hyperparameter_Target_All = pd.DataFrame(list(Best_Hyperparameter_Target_All),
                          columns =['learning_rate', 'max_depth', 'num_leaves', 'min_data_in_leaf','reg_lambda'],
                          index =['data1'])                          
                          # index =['data1', 'data2', 'data3', 'data4', 'data5', 'data6', 'data7'])
Best_Hyperparameter_Target_All.to_csv(Path2 + RunCondition + '_Best_Hyperparameter_Target.csv' ,header=True , index=True)



#%%
from sklearn.metrics import mean_squared_error
print("RMSE of the validation set:", np.sqrt(mean_squared_error(cdata_test[keys2[-1]],predicted2)))


#%%
cdata_testAns = TestAns_dataset[0]
keys2Ans = cdata_testAns.keys()

cdata_testAns_Target = cdata_testAns[keys2Ans[-1]]
cdata_testAns_Feature = cdata_testAns[keys[:-1]]

# Calculate accuracy

MAE_test2 = metrics.mean_absolute_error(cdata_testAns_Target, predicted2 )
print('MAE of LGBM model: ', MAE_test2)

RMSD_test2 = metrics.mean_squared_error(cdata_testAns_Target, predicted2 )
print('RMSD of LGBM model: ', RMSD_test2)

R2_test2 = metrics.r2_score(cdata_testAns_Target, predicted2 )
print('R2 of LGBM model: ', R2_test2)

score_test2 = LGBM.score(cdata_testAns_Feature, cdata_testAns_Target)
print('Accuracy of LGBM model: ', score_test2)

## Calculate accuracy
#accuracy = metrics.accuracy_score(test_target, predicted)
#print('Accuracy of testing dataset: ', accuracy)
## np.r_[cdata_feature_array,predicted]
import numpy as np
cdata_Test2_MAE = []
cdata_Test2_RMSD = []
cdata_Test2_R2 = []
cdata_Test2_score = []

cdata_Test2_MAE.append(MAE_test2)
cdata_Test2_RMSD.append(RMSD_test2)
cdata_Test2_R2.append(R2_test2)
cdata_Test2_score.append(score_test2)


Accuracy_df = pd.DataFrame(list(zip( cdata_Test2_MAE,cdata_Test2_RMSD,cdata_Test2_R2,cdata_Test2_score)),
           columns =[ 'Test_MAE','Test_RMSD', 'Test_R2', 'Test_Score'], 
           #index =['data1', 'data2', 'data3', 'data4', 'data5', 'data6', 'data7'])    
           index =['data1'])    
Accuracy_df.to_csv(Path2 + RunCondition + '_AccuracyFinal_Test.csv' ,header=True , index=True)


#%%
## 找重要特徵
import lightgbm as lgb
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(10, 7))
lgb.plot_importance(LGBM, max_num_features=30, ax=ax, importance_type ='gain')
plt.title("LightGBM - Feature Importance");


fig, ax = plt.subplots(figsize=(10, 7))
lgb.plot_importance(LGBM, max_num_features=30, ax=ax)
plt.title("LightGBM - Feature Importance");

## 畫出樹
from lightgbm import plot_tree

lgb.create_tree_digraph(LGBM, tree_index=0, show_info=None, precision=3, orientation='horizontal')
fig, ax = plt.subplots(figsize=(10, 7))
lgb.plot_tree(LGBM, ax=None, tree_index=0, figsize=None, dpi=None, show_info=None, precision=3, orientation='horizontal')
plt.title("LightGBM - PrintTree");


boost = LGBM.booster_
print('Feature names',boost.feature_name())
