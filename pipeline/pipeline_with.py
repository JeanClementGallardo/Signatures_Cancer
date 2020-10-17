'''
Jean-Clément GALLARDO
14/05/2020
Pipeline V2
'''


import numpy as np
import pandas as pd
import os
import csv

from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder

from sklearn.metrics import balanced_accuracy_score
from sklearn.metrics import confusion_matrix

import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

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

def import_file(path):
    script_dir = os.path.dirname("multinomial_logistic_regression")
    pam50_path = os.path.join(script_dir, path)
    f = open(pam50_path, "r")
    pam50 = f.read()
    f.close()
    list_pam50 = pam50.split("\n")
    return list_pam50

# Create and fit model
def create_rf(geneTrain):
    random_forest = RandomForestClassifier(
        criterion = 'entropy', n_estimators = 1000
    )
    random_forest.fit(geneTrain.loc[:, geneTrain.columns != 'subtype'], geneTrain.subtype)
    return random_forest

def create_all_model(geneTrain,number):
    list_model = []
    for idx in range(number):
        tmp = create_rf(geneTrain)
        list_model.append(tmp)
    return list_model

def create_logistic(newTrain):
    logistic = LogisticRegression(
        C = 0.5, tol = 0.001, penalty = 'l1', solver = "saga", max_iter = 2000
    )
    logistic.fit(newTrain.loc[:, newTrain.columns != 'subtype'], newTrain.subtype)
    return logistic

# Select features from models
def select_features(list_model,data, list_pam50):
    genes = data.columns[1:]
    features_selection = []
    list_genes = {}

    for idx in range(len(list_model)):
        tmp_selection = model_value(list_model[idx], data, list_pam50)
        for gene in range(len(tmp_selection)):
            if(tmp_selection[gene] in list_genes.keys()):
                list_genes[tmp_selection[gene]]+=1
            if(tmp_selection[gene] not in list_genes.keys()):
                list_genes[tmp_selection[gene]] = 1

    for gene in list_genes.keys():
        if (list_genes[gene] > (len(list_model)/2)):
            features_selection.append(gene)

    return features_selection

# Select features from the median of PAM50 coefs
def features_selection(genes_coef):
    features_selected = []

    pam50_coef = []

    for pam50 in genes_coef["PAM50"].keys():
        pam50_coef.append(genes_coef["PAM50"][pam50])
        features_selected.append(pam50)

    median = np.percentile(pam50_coef, 50)

    for others in genes_coef["OTHERS"].keys():
        if(genes_coef["OTHERS"][others] >= median):
            features_selected.append(others)

    print(len(features_selected)-50,"genes selected without PAM50 label")

    return features_selected

def model_value(model,data,list_pam50):
    genes_coef = {"PAM50" : {}, "OTHERS" : {}}
    coefficients = model.feature_importances_
    genes = data.columns[1:]

    for idx in range(len(coefficients)):
        if (genes[idx] in list_pam50 and coefficients[idx] != 0):
            genes_coef["PAM50"][genes[idx]] = coefficients[idx]
        if (genes[idx] not in list_pam50 and coefficients[idx] != 0):
            genes_coef["OTHERS"][genes[idx]] = coefficients[idx]

    tmp_list = features_selection(genes_coef)
    return tmp_list

# Create new dataset with selected features
def create_dataset(data,features_selected):
    features_selected.append("subtype")
    newData = data.loc[:, features_selected]

    return newData

# Export final results
def export_result(model, data, geneTest, filename):
    subtypes = encoder.inverse_transform(model.classes_)
    sparsity = np.mean(model.coef_ == 0) * 100
    score = model.score(geneTest.loc[:, geneTest.columns != 'subtype'], geneTest.subtype)
    balanced_score = balanced_accuracy_score(geneTest.subtype, model.predict(geneTest.loc[:, geneTest.columns != 'subtype']))
    matrix = confusion_matrix(encoder.inverse_transform(geneTest.subtype),
                              encoder.inverse_transform(model.predict(geneTest.loc[:, geneTest.columns != 'subtype'])),
                              labels = subtypes)
    coefficient = model.coef_
    genes = data.columns[1:]
    export_dict = {"sparsity" : sparsity, "score" : score,
                    "balanced_score" : balanced_score, "confusion_matrix" : matrix,
                    subtypes[0] : {}, subtypes[1] : {}, subtypes[2] : {},
                    subtypes[3] : {}, subtypes[4] : {} }

    for subtype in range(len(coefficient)):
        for idx in range(len(coefficient[subtype])):
            export_dict[subtypes[subtype]][genes[idx]] = coefficient[subtype][idx]

    with open(filename, 'w') as csvfile :
        writer = csv.DictWriter(csvfile, fieldnames = export_dict.keys())
        writer.writeheader()
        writer.writerow(export_dict)

def check_selection(features_selected, list_pam50):
    tmp = []
    print('number of genes selected : ',len(features_selected))
    for i in range(len(features_selected)):
        if(features_selected[i] in list_pam50):
            tmp.append(features_selected[i])
    print("with : ",len(tmp),'from PAM50 list')

if __name__ == "__main__":
    encoder = LabelEncoder()

    ## Import data
    data = import_csv("../data/data.csv", encoder)
    list_pam50 = import_file("../data/PAM50_geneIDs.txt")
    print(len(list_pam50))

    ## Split data into training and test dataset
    geneTrain, geneTest= train_test_split(data, train_size = 700, random_state = 69780)
    print("data split")
    ## Create and fit models
    list_rf = create_all_model(geneTrain, 30)
    print("model created")
    ## Show results from Random Forest model
    sparsity = 0
    score = 0
    balanced_score = 0

    for idx in range(len(list_rf)):
        sparsity+=np.mean(list_rf[idx].feature_importances_ == 0) * 100
        score+=list_rf[idx].score(geneTest.loc[:, geneTest.columns != 'subtype'], geneTest.subtype)
        balanced_score+=balanced_accuracy_score(geneTest.subtype, list_rf[idx].predict(geneTest.loc[:, geneTest.columns != 'subtype']))

    print("Results from Random Forest model")
    print("score : ",score/len(list_rf))
    print("balanced score : ",balanced_score/len(list_rf))
    print("sparsity : ",sparsity/len(list_rf))

    ## Features selection
    features_selected = select_features(list_rf, data, list_pam50)
    print("features selected")
    check_selection(features_selected,list_pam50)
    ## Create new dataset from selected features
    newData = create_dataset(data,features_selected)
    print("New dataset created")
    ## Split new dataset
    newTrain, newTest = train_test_split(newData, train_size = 700, random_state = 69780)
    ## Create logistic model with selected features
    logistic = create_logistic(newTrain)
    print("logistic model created")
    ## Export final result
    export_result(logistic, newData, newTest, "with_30model.csv")
    print("result exported")
