'''
Jean-Clément GALLARDO
12/05/2020
k-NN : Classifier
'''

import numpy as np
import pandas as pd
import os

from sklearn.neighbors import KNeighborsClassifier

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
def params_seeker(geneTrain, list_weight, list_nb_neighbors, list_p, list_leaf):
    model = KNeighborsClassifier()
    parameters = {'weights' : list_weight, 'n_neighbors' : list_nb_neighbors, 'p' : list_p, 'leaf_size' : list_leaf}
    model_cv = GridSearchCV(model, parameters)
    model_cv.fit(geneTrain.loc[:, geneTrain.columns != 'subtype'], geneTrain.subtype)
    sorted(model_cv.cv_results_.keys())
    print(model_cv.cv_results_)

if __name__ == "__main__":
    encoder = LabelEncoder()

    #Import data
    data = import_csv("../data/data.csv", encoder)

    #Split data : train and test
    geneTrain, geneTest= train_test_split(data, train_size = 700, random_state = 69780)

    #Found best parameters
    list_p = [1, 2]
    list_leaf = [10,20,30,50,100,250,500,1000]
    list_nb_neighbors = [1,3,5,10,20,30,50]
    list_weight = ['uniform','distance']

    params_seeker(geneTrain, list_weight, list_nb_neighbors, list_p, list_leaf)
