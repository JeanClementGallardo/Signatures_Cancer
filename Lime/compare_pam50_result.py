import numpy as np
import pandas as pd
import os
import csv

from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import balanced_accuracy_score

import matplotlib.pyplot as plt
import seaborn as sbn

## Function ##

# Import data and encode clinical label
def import_csv(path, encoder):
    ## Import
    data = pd.read_csv(path)
    data.drop(["patient_id"], axis = 1, inplace = True)
    data = data.dropna()
    ## Encode
    categorical_feature_mask = data.dtypes==object
    categorical_cols = data.columns[categorical_feature_mask].tolist()
    data[categorical_cols] = data[categorical_cols].apply(lambda col: encoder.fit_transform(col))
    return data

# Create and fit model
def create_svm(geneTrain, C, tol):
    model = SVC(
        C = C, tol = tol, kernel = "linear"
    )
    model.fit(geneTrain.loc[:, geneTrain.columns != 'subtype'], geneTrain.subtype)
    return model

def create_logistic(geneTrain, C_value, tol_value):
    logistic_model = LogisticRegression(
        C = C_value, penalty='l1', solver='saga', tol=tol_value, max_iter = 2000
    )
    logistic_model.fit(geneTrain.loc[:, geneTrain.columns != 'subtype'], geneTrain.subtype)
    return logistic_model


if __name__ == "__main__":
    encoder = LabelEncoder()

    ## Import data
    data = import_csv("../data/data_pam50.csv", encoder)
    remove_svm = ["MKI67","KIF2C","TYMS","UBE2C","MYC","CCNE1",
                   "ANLN","RRM2","BLVRA","MDM2","GRB7","ACTR3B"]

    remove_logistic = ["GRB7","BLVRA","FGFR4","MDM2","CCNE1","RRM2", "RRM2",
                   "CEP55","ERBB2","KIF2C","ACTR3B","TYMS","MKI67"]

    data_logistic = data.drop(remove_logistic,axis = 1)
    data_svm = data.drop(remove_svm,axis = 1)

    ## Split data into training and test dataset
    logisticTrain, logisticTest= train_test_split(data_logistic, train_size = 700, random_state = 69780)
    svmTrain, svmTest= train_test_split(data_svm, train_size = 700, random_state = 69780)
    geneTrain, geneTest= train_test_split(data, train_size = 700, random_state = 69780)
    print("data split")
    ## Create and fit models
    svm_remove = create_svm(svmTrain, 0.005, 0.001)
    logistic_remove = create_logistic(logisticTrain, 0.5, 0.001)
    svm = create_svm(geneTrain, 0.005, 0.001)
    logistic = create_logistic(geneTrain, 0.5, 0.001)
    print("model created")

    score_svm = svm.score(geneTest.loc[:, geneTest.columns != 'subtype'], geneTest.subtype)
    balanced_svm = balanced_accuracy_score(geneTest.subtype, svm.predict(geneTest.loc[:, geneTest.columns != 'subtype']))

    score_logistic = logistic.score(geneTest.loc[:, geneTest.columns != 'subtype'], geneTest.subtype)
    balanced_logistic = balanced_accuracy_score(geneTest.subtype, logistic.predict(geneTest.loc[:, geneTest.columns != 'subtype']))

    score_svm40 = svm_remove.score(svmTest.loc[:, svmTest.columns != 'subtype'], svmTest.subtype)
    balanced_svm40 = balanced_accuracy_score(svmTest.subtype, svm_remove.predict(svmTest.loc[:, svmTest.columns != 'subtype']))

    score_logistic40 = logistic_remove.score(logisticTest.loc[:, logisticTest.columns != 'subtype'], logisticTest.subtype)
    balanced_logistic40 = balanced_accuracy_score(logisticTest.subtype, logistic_remove.predict(logisticTest.loc[:, logisticTest.columns != 'subtype']))

    barWidth = 0.4
    y1 = [score_svm, score_svm40, score_logistic, score_logistic40]
    y2 = [balanced_svm, balanced_svm40, balanced_logistic, balanced_logistic40]
    print(y1)
    print(y2)
    r1 = range(len(y1))
    r2 = [x + barWidth for x in r1]

    ax = plt.gca()
    ax.set_ylim(0,1)

    plt.bar(r1, y1, width = barWidth, color = ['red' for i in y1])
    plt.bar(r2, y2, width = barWidth, color = ['black' for i in y1])
    plt.xticks([r + barWidth / 2 for r in range(len(y1))], ["SVM", "SVM40", "LRM", "LRM40"])
    plt.show()

    print("Process finished")
