'''
Jean-Clément GALLARDO
12/05/2020
k-NN : Classifier
Best params :
{'leaf_size': 20, 'n_neighbors': 1, 'p': 2, 'weights': 'uniform'}
{'leaf_size': 1000, 'n_neighbors': 3, 'p': 1, 'weights': 'uniform'}
{'leaf_size': 20, 'n_neighbors': 3, 'p': 1, 'weights': 'uniform'}
'''
import numpy as np
import pandas as pd
import os
import csv

from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder

from sklearn.metrics import balanced_accuracy_score
from sklearn.metrics import confusion_matrix

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

# Create and fit logistic model
def create_model(geneTrain, list_params):
    model = KNeighborsClassifier(
        leaf_size = list_params[0], n_neighbors = list_params[1], p = list_params[2],
        weights = list_params[3]
    )
    model.fit(geneTrain.loc[:, geneTrain.columns != 'subtype'], geneTrain.subtype)
    return model

if __name__ == "__main__":
    encoder = LabelEncoder()

    ## Import data
    data = import_csv("../data/data.csv", encoder)

    ## Split data into training and test dataset
    geneTrain, geneTest= train_test_split(data, train_size = 700, random_state = 69780)
    print("data split")
    ## Create and fit models
    params_list1 = [20, 1, 2, 'uniform']
    params_list2 = [1000, 3, 1, 'uniform']
    params_list3 = [20, 3, 1, 'uniform']

    model1 = create_model(geneTrain, params_list1)
    model2 = create_model(geneTrain, params_list2)
    model3 = create_model(geneTrain, params_list3)
    subtypes = encoder.inverse_transform(model1.classes_)
    print("model created")
    print("results model params 1 :")

    score = model1.score(geneTest.loc[:, geneTest.columns != 'subtype'], geneTest.subtype)
    balanced_score = balanced_accuracy_score(geneTest.subtype, model1.predict(geneTest.loc[:, geneTest.columns != 'subtype']))
    matrix = confusion_matrix(encoder.inverse_transform(geneTest.subtype),
                              encoder.inverse_transform(model1.predict(geneTest.loc[:, geneTest.columns != 'subtype'])),
                              labels = subtypes)

    print("score : ",score)
    print("balanced score : ",balanced_score)
    print("confusion matrix : \n",matrix)

    print("results model params 2 :")

    score = model2.score(geneTest.loc[:, geneTest.columns != 'subtype'], geneTest.subtype)
    balanced_score = balanced_accuracy_score(geneTest.subtype, model2.predict(geneTest.loc[:, geneTest.columns != 'subtype']))
    matrix = confusion_matrix(encoder.inverse_transform(geneTest.subtype),
                              encoder.inverse_transform(model2.predict(geneTest.loc[:, geneTest.columns != 'subtype'])),
                              labels = subtypes)

    print("score : ",score)
    print("balanced score : ",balanced_score)
    print("confusion matrix : \n",matrix)

    print("results model params 3 :")

    score = model3.score(geneTest.loc[:, geneTest.columns != 'subtype'], geneTest.subtype)
    balanced_score = balanced_accuracy_score(geneTest.subtype, model3.predict(geneTest.loc[:, geneTest.columns != 'subtype']))
    matrix = confusion_matrix(encoder.inverse_transform(geneTest.subtype),
                              encoder.inverse_transform(model3.predict(geneTest.loc[:, geneTest.columns != 'subtype'])),
                              labels = subtypes)

    print("score : ",score)
    print("balanced score : ",balanced_score)
    print("confusion matrix : \n",matrix)
