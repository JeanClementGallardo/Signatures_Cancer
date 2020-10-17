'''
JeanClement GALLARDO
09/07/2020
Comparaison des résultats explicatifs de LIME sur deux modèles différents
'''

import numpy as np
import pandas as pd
import os
import csv

import lime
import lime.lime_tabular
import matplotlib.pyplot as plt
import seaborn as sbn

from sklearn.linear_model import LogisticRegression
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
def create_logistic(X_train, Y_train, C_value, tol_value):
    logistic_model = LogisticRegression(
        C = C_value, penalty='l1', solver='saga', tol=tol_value, max_iter = 2000
    )
    logistic_model.fit(X_train, Y_train)
    return logistic_model

def create_svm(X_train, Y_train, C, tol):
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

    logistic = create_logistic(X_train,Y_train, 0.5,0.001)
    svm = create_svm(X_train,Y_train,0.5,0.001)
    print("Model created")

    training = X_train.to_numpy()
    test = X_test.to_numpy()
    features_names = data.columns[1:]
    subtype = encoder.inverse_transform(logistic.classes_)
    sub = Y_test.to_numpy()
    sub_decode = encoder.inverse_transform(sub)
    mat_lrm = {}
    mat_svm = {}

    for sample in range(len(test)):
        tmp_logistic = lime_explainer(logistic, training, test, sample, features_names, subtype, sub[sample])
        tmp_svm = lime_explainer(svm, training, test, sample, features_names, subtype, sub[sample])
        mat_lrm[sample] = {}
        mat_svm[sample] = {}
        for idx in range(len(tmp_logistic)):
            split_row = tmp_logistic[idx][0].split(" ")
            value = tmp_logistic[idx][1]
            j = 0
            while(j < len(split_row)):
                if (split_row[j].isalnum() == False):
                    split_row.pop(j)
                else :
                    j+=1
            gene = str(split_row)
            mat_lrm[sample][gene] = float(value)
        for idx in range(len(tmp_svm)):
            split_row = tmp_svm[idx][0].split(" ")
            value = tmp_svm[idx][1]
            j = 0
            while(j < len(split_row)):
                if (split_row[j].isalnum() == False):
                    split_row.pop(j)
                else :
                    j+=1
            gene = str(split_row)
            mat_svm[sample][gene] = float(value)


    x_label = sorted(mat_lrm[1].keys())

    lrm = []
    svm = []
    maxi = 0
    for sample in mat_lrm.keys():
        gene = mat_lrm[sample].keys()
        gene = sorted(gene)
        tmp_lrm = []
        tmp_svm = []
        for idx in range(len(gene)):
            tmp_lrm.append(mat_lrm[sample][gene[idx]])
            tmp_svm.append(mat_svm[sample][gene[idx]])
            if (maxi < mat_lrm[sample][gene[idx]]):
                maxi = mat_lrm[sample][gene[idx]]
        lrm.append(tmp_lrm)
        svm.append(tmp_svm)

    mini = -maxi

    i = 0
    while(i < len(lrm)):
        if (lrm[i] == []):
            lrm.pop(i)
        else :
            i+=1

    i = 0
    while(i < len(svm)):
        if (svm[i] == []):
            svm.pop(i)
        else :
            i+=1

    mat_final = lrm

    for idx in range(len(svm)):
        mat_final.append(svm[idx])

    print("begin compute heatmap")
    sbn.set(color_codes = True)

    dict_map = {"Her2" : "blue", "Basal" : "red", "Normal" : "green", "LumA" : "yellow", "LumB" : "black"}
    row_colors = []
    for idx in range(len(sub_decode)):
        row_colors.append(dict_map[sub_decode[idx]])

    test_sous_type = row_colors*2

    fg = sbn.clustermap(mat_final, row_colors = test_sous_type, xticklabels = x_label,
        cmap = "twilight_shifted", vmin=mini, vmax=maxi)

    plt.show()
