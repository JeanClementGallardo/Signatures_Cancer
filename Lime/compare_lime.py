'''
Jean-Clement GALLARDO
13/06/2020
test lime package with logistic regression
Comparison between PAM50 and All genome
'''


import numpy as np
import pandas as pd
import os
import csv

import lime
import lime.lime_tabular
import matplotlib.pyplot as plt

from sklearn.linear_model import LogisticRegression
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
def create_model(geneTrain, C_value, tol_value):
    logistic_model = LogisticRegression(
        C = C_value, penalty='l1', solver='saga', tol=tol_value, max_iter = 2000
    )
    logistic_model.fit(geneTrain.loc[:, geneTrain.columns != 'subtype'], geneTrain.subtype)
    return logistic_model

#Lime process
def lime_histo(numpy_data, test, features_names, subtype, encoder, model, sample, sub, num):
    explainer = lime.lime_tabular.LimeTabularExplainer(numpy_data,
                feature_names = features_names, class_names = subtype)
    exp = explainer.explain_instance(test[sample], model.predict_proba, labels=[0,1,2,3,4], num_features=num)
    return exp.as_list(label=sub)

def compare_lime(explainer_pam50, explainer_genome, logistic_pam50,
            logistic_genome, numpy_subtype_pam50, numpy_data_pam50,
            numpy_data_genome, encoder_pam50, encoder_genome, subtype):
    compare_result = {}
    for sample in range(150):
        exp_pam50 = lime_histo(numpy_data_pam50, numpy_test_pam50, features_names_pam50, subtype, encoder_pam50, logistic_pam50, sample, numpy_subtype_pam50[sample], 50)
        exp_genome = lime_histo(numpy_data_genome, numpy_test_genome, features_names_genome, subtype, encoder_genome, logistic_genome, sample, numpy_subtype_pam50[sample], 2000)
        compare_result[sample] = {"pam50" : [], "genome" : []}
        for i in range(len(exp_pam50)):
            pam50 = exp_pam50[i][0].split(" ")
            j = 0
            while(j < len(pam50)):
                if (pam50[j].isalnum() == False):
                    pam50.pop(j)
                else :
                    j+=1
            compare_result[sample]["pam50"].append(pam50)

        for i in range(len(exp_genome)):
            genome = exp_genome[i][0].split(" ")
            j = 0
            while(j < len(genome)):
                if (genome[j].isalnum() == False):
                    genome.pop(j)
                else :
                    j+=1
            compare_result[sample]["genome"].append(genome)

    return compare_result

def pourcent_corres(compare_result):
    pourcent_score = []
    for sample in compare_result.keys():
        tmp_score = 50
        for idx in range(len(compare_result[sample]["pam50"])):
            if (compare_result[sample]["pam50"][idx] not in compare_result[sample]["genome"]):
                tmp_score-=1
        pourcent_score.append(tmp_score)
    return pourcent_score


if __name__ == "__main__":
    encoder_pam50 = LabelEncoder()
    encoder_genome = LabelEncoder()

    ## Import data
    pam50 = import_csv("../data/data_pam50.csv", encoder_pam50)
    genome = import_csv("../data/data.csv", encoder_genome)

    ## Split data into training and test dataset
    geneTrain_pam50, geneTest_pam50 = train_test_split(pam50, train_size = 700, random_state = 69780)
    X_train_pam50 = geneTrain_pam50.loc[:, geneTrain_pam50.columns != 'subtype']
    Y_train_pam50 = geneTrain_pam50.subtype
    X_test_pam50 = geneTest_pam50.loc[:, geneTest_pam50.columns != 'subtype']
    Y_test_pam50 = geneTest_pam50.subtype

    geneTrain_genome, geneTest_genome = train_test_split(genome, train_size = 700, random_state = 69780)
    X_train_genome = geneTrain_genome.loc[:, geneTrain_genome.columns != 'subtype']
    Y_train_genome = geneTrain_genome.subtype
    X_test_genome = geneTest_genome.loc[:, geneTest_genome.columns != 'subtype']
    Y_test_genome = geneTest_genome.subtype

    print("data split")


    ## Create and fit models
    params = [0.5,0.001]
    logistic_pam50 = create_model(geneTrain_pam50, params[0], params[1])
    logistic_genome = create_model(geneTrain_genome, params[0], params[1])
    print("model created")


    ## Lime results
    numpy_data_pam50 = X_train_pam50.to_numpy()
    numpy_subtype_pam50 = Y_train_pam50.to_numpy()
    numpy_test_pam50 = X_test_pam50.to_numpy()
    features_names_pam50 = pam50.columns[1:]
    subtype_pam50 = logistic_pam50.classes_
    explainer_pam50 = lime.lime_tabular.LimeTabularExplainer(numpy_data_pam50,
                feature_names = features_names_pam50, class_names = subtype_pam50)

    numpy_data_genome = X_train_genome.to_numpy()
    numpy_subtype_genome = Y_train_genome.to_numpy()
    numpy_test_genome = X_test_genome.to_numpy()
    features_names_genome = genome.columns[1:]
    subtype_genome = logistic_genome.classes_
    explainer_genome = lime.lime_tabular.LimeTabularExplainer(numpy_data_genome,
                feature_names = features_names_genome, class_names = subtype_genome)

    result = compare_lime(explainer_pam50, explainer_genome, logistic_pam50,
                logistic_genome, numpy_subtype_pam50, numpy_data_pam50,
                numpy_data_genome, encoder_pam50, encoder_genome, subtype_pam50)

    pourcent_score = pourcent_corres(result)
    plt.hist(pourcent_score, range=(0, 50), color = 'red', edgecolor = 'black')
    plt.xlabel('Nombres de features sélectionnées par Lime qui sont identiques entre les deux jeux de données')
    plt.ylabel('Nombres de samples')
    plt.show()

    print(result)

    print("Process finished")
