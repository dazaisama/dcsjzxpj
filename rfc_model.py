# -*- coding: utf-8 -*-
"""
Created on Wed Nov  1 13:09:11 2023

@author: 86188
"""

#随机森林分类

import pandas as pd
from sklearn import metrics  # 分类结果评价函数
import pickle
import joblib




def dt_model21(data):
# 新的数据

    x1 = data[[ 'x1','x2','x3','x4','x5','x6','x7','x8','x9','x10','x11']]  # 特征
    y1 = data['y']  # 标签
    #调用模型
    loaded_rfc = joblib.load('app01/random_forest_model_r.pkl')
    predicted1 = loaded_rfc.predict(x1)
    accuracy1 = metrics.accuracy_score(y1, predicted1)  # 求精度
    print("level: %.2f%%" % ((1-accuracy1) * 100.0))
    level=(1-accuracy1) * 100.0
    return predicted1,level


def dt_model22(data):
# 新的数据
    x1 = data[[ 'x1','x2','x3','x4','x5','x6','x7','x8','x9','x10','x11']]  # 特征
    y1 = data['y']  # 标签
    
    #调用模型
    loaded_rfc = joblib.load('app01/random_forest_model_p.pkl')
    predicted1 = loaded_rfc.predict(x1)
    accuracy1 = metrics.accuracy_score(y1, predicted1)  # 求精度
    print("level: %.2f%%" % ((1-accuracy1) * 100.0))
    level=(1-accuracy1) * 100.0
    return predicted1,level