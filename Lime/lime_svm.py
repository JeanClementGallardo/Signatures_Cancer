'''
Jean-Clement GALLARDO
19/06/2020
test lime package with svm
'''


import numpy as np
import pandas as pd
import os
import csv

import lime
import lime.lime_tabular
import matplotlib.pyplot as plt
import seaborn as sbn

from sklearn.svm import SVC
from sklearn.model_selection import train_test_split
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

# Create and fit logistic model
def create_model(X_train, Y_train, C, tol):
    model = SVC(
        C = C, tol = tol, kernel = "linear", probability = True
    )
    model.fit(X_train, Y_train)
    return model

## Create Lime explainer
def lime_explainer(model, training, test, sample, features_names, subtype, sub):
    explainer = lime.lime_tabular.LimeTabularExplainer(training, feature_names = features_names, class_names = subtype)
    exp = explainer.explain_instance(test[sample], model.predict_proba, num_features = 50, labels = [0,1,2,3,4])
    return exp.as_list(label = sub)

if __name__ == "__main__":
    encoder = LabelEncoder()
    data = import_csv("../data/data_pam50.csv", encoder)
    print("Data imported")

    geneTrain, geneTest = train_test_split(data, train_size = 700, random_state = 69780)
    X_train = geneTrain.loc[:, geneTrain.columns != 'subtype']
    Y_train = geneTrain.subtype
    X_test = geneTest.loc[:, geneTest.columns != 'subtype']
    Y_test = geneTest.subtype
    print("Data split")

    svm = create_model(X_train,Y_train, 0.5,0.001)
    print("Model created")

    training = X_train.to_numpy()
    test = X_test.to_numpy()
    features_names = data.columns[1:]
    subtype = encoder.inverse_transform(svm.classes_)
    sub = Y_test.to_numpy()
    sub_decode = encoder.inverse_transform(sub)
    mat_result = {}

    for sample in range(len(test)):
        tmp_result = lime_explainer(svm, training, test, sample, features_names, subtype, sub[sample])
        mat_result[sample] = {}
        for idx in range(len(tmp_result)):
            split_row = tmp_result[idx][0].split(" ")
            value = tmp_result[idx][1]
            j = 0
            while(j < len(split_row)):
                if (split_row[j].isalnum() == False):
                    split_row.pop(j)
                else :
                    j+=1
            gene = str(split_row)
            mat_result[sample][gene] = float(value)

    y_label = mat_result.keys()
    x_label = sorted(mat_result[1].keys())

    mat = []
    maxi = 0
    for sample in mat_result.keys():
        gene = mat_result[sample].keys()
        gene = sorted(gene)
        tmp = []
        for idx in range(len(gene)):
            tmp.append(mat_result[sample][gene[idx]])
            if (maxi < mat_result[sample][gene[idx]]):
                maxi = mat_result[sample][gene[idx]]
        mat.append(tmp)

    mini = -maxi

    i = 0
    while(i < len(mat)):
        if (mat[i] == []):
            mat.pop(i)
        else :
            i+=1

    print("begin compute heatmap")
    sbn.set(color_codes = True)

    dict_map = {"Her2" : "blue", "Basal" : "red", "Normal" : "green", "LumA" : "yellow", "LumB" : "black"}
    row_colors = []
    for idx in range(len(test)):
        row_colors.append(dict_map[sub_decode[idx]])

    fg = sbn.clustermap(mat, row_colors = row_colors, xticklabels = x_label,
        cmap = "twilight_shifted", vmin=mini, vmax=maxi)

    plt.show()
