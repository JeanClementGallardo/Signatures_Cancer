'''
Jean-Clément GALLARDO
03/04/2020
Multinomial Logistic Regression parameters' seeker
'''

import numpy as np
import pandas as pd
import os

from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import LogisticRegressionCV
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV
from sklearn.preprocessing import LabelEncoder


## Function ##
# Import data and encode clinical label
def import_csv(path, encoder):
    ## Import
    script_dir = os.path.dirname("multinomial_logistic_regression")
    data_path = os.path.join(script_dir, path)
    data = pd.read_csv(data_path)
    data.drop(["patient_id"], axis = 1, inplace = True)
    data = data.dropna(subset=["subtype"])
    data = data.fillna(0)
    ## Encode
    categorical_feature_mask = data.dtypes==object
    categorical_cols = data.columns[categorical_feature_mask].tolist()
    data[categorical_cols] = data[categorical_cols].apply(lambda col: encoder.fit_transform(col))
    return data

# Search best parameter for model
def params_seeker(geneTrain, list_tol, list_C):
    logistic = LogisticRegression(solver = 'saga', penalty = 'l1', max_iter = 2000)
    parameters = {'tol' : list_tol, 'C' : list_C}
    logistic_cv = GridSearchCV(logistic, parameters)
    logistic_cv.fit(geneTrain.loc[:, geneTrain.columns != 'subtype'], geneTrain.subtype)
    sorted(logistic_cv.cv_results_.keys())
    print(logistic_cv.cv_results_)

if __name__ == "__main__":
    encoder = LabelEncoder()

    #Import data
    data = import_csv("../data/data_article.csv", encoder)

    #Split data : train and test
    geneTrain, geneTest= train_test_split(data, train_size = 120, random_state = 69780)

    #Found best parameters
    list_C = [0.005,0.008,0.01,0.025,0.05,0.1,0.5]
    list_tol = [0.001,0.01,0.025,0.05,0.1,0.2]

    params_seeker(geneTrain, list_tol, list_C)
