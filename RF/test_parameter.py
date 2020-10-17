'''
Jean-Clément GALLARDO
15/04/2020
Random Forest parameters' seeker
'''

import numpy as np
import pandas as pd
import os

from sklearn.ensemble import RandomForestClassifier
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
    data = data.dropna()
    ## Encode
    categorical_feature_mask = data.dtypes==object
    categorical_cols = data.columns[categorical_feature_mask].tolist()
    data[categorical_cols] = data[categorical_cols].apply(lambda col: encoder.fit_transform(col))
    return data

# Search best parameter for model
def params_seeker(geneTrain, list_estimator, list_criterion):
    random_forest = RandomForestClassifier()
    parameters = {"criterion" : list_criterion, "n_estimators" : list_estimator}
    RF_cv = GridSearchCV(random_forest, parameters)
    RF_cv.fit(geneTrain.loc[:, geneTrain.columns != 'subtype'], geneTrain.subtype)
    sorted(RF_cv.cv_results_.keys())
    print(RF_cv.cv_results_)

if __name__ == "__main__":
    encoder = LabelEncoder()

    #Import data
    data = import_csv("../data/data.csv", encoder)

    #Split data : train and test
    geneTrain, geneTest= train_test_split(data, train_size = 700, random_state = 69780)

    #Found best parameters
    list_criterion = ["gini", "entropy"]
    list_estimator = [20,40,60,80,100,200,500,800,1000]

    params_seeker(geneTrain, list_estimator, list_criterion)
