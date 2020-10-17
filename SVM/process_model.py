'''
Jean-Clément GALLARD
12/05/2020
Support Vector Machine : Classifier
Best params :
{'C': 0.005, 'kernel': 'linear', 'tol': 0.001}
{'C': 0.005, 'kernel': 'linear', 'tol': 0.2}
{'C': 1, 'kernel': 'linear', 'tol': 0.2}
'''

import numpy as np
import pandas as pd
import os
import csv

from sklearn.svm import SVC
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

# Create and fit model
def create_model(geneTrain, list_params):
    model = SVC(
        C = list_params[0], tol=list_params[1], kernel = "linear"
    )
    model.fit(geneTrain.loc[:, geneTrain.columns != 'subtype'], geneTrain.subtype)
    return model

# Create one list of coef
def unique_coef(list_coef):
    unique_coef = []
    for one in range(len(list_coef)):
        if(one == 0):
            for two in range(len(list_coef[one])):
                unique_coef.append(list_coef[one][two])
        if(one != 0):
            for two in range(len(list_coef[one])):
                unique_coef[two] = unique_coef[two] + list_coef[one][two]

    for idx in range(len(unique_coef)):
        unique_coef[idx] = unique_coef[idx]/len(list_coef)

    return unique_coef

# Calculate sparsity
def sparse(coef):
    cpt = 0
    for i in range(len(coef)):
        if (coef[i] == 0):
            cpt+=1
    cpt = (cpt/len(coef)) * 100

    return cpt

# Export results
def export_result(model, data, geneTest,filename):
    score = model.score(geneTest.loc[:, geneTest.columns != 'subtype'], geneTest.subtype)
    coefficient = model.coef_
    coef = unique_coef(coefficient)
    sparsity = sparse(coef)
    genes = data.columns[1:]
    balanced_score = balanced_accuracy_score(geneTest.subtype, model.predict(geneTest.loc[:, geneTest.columns != 'subtype']))
    matrix = confusion_matrix(encoder.inverse_transform(geneTest.subtype),
                              encoder.inverse_transform(model.predict(geneTest.loc[:, geneTest.columns != 'subtype'])))
    export_dict = {"sparsity" : sparsity, "score" : score,
                    "coefficient" : {}, "confusion_matrix": matrix, "balanced_score": balanced_score }

    for i in range(len(coef)):
        export_dict["coefficient"][genes[i]] = coef[i]

    with open(filename, 'w') as csvfile :
        writer = csv.DictWriter(csvfile, fieldnames = export_dict.keys())
        writer.writeheader()
        writer.writerow(export_dict)



if __name__ == "__main__":
    encoder = LabelEncoder()

    ## Import data
    data = import_csv("../data/data.csv", encoder)

    ## Split data into training and test dataset
    geneTrain, geneTest= train_test_split(data, train_size = 700, random_state = 69780)
    print("data split")
    ## Create and fit models
    params_list1 = [0.005,0.001]
    params_list2 = [0.005,0.2]
    params_list3 = [1,0.2]

    model1 = create_model(geneTrain, params_list1)
    model2 = create_model(geneTrain, params_list2)
    model3 = create_model(geneTrain, params_list3)
    print("model created")
    ## Export results
    export_result(model1, data, geneTest, 'result1.csv')
    export_result(model2, data, geneTest, 'result2.csv')
    export_result(model3, data, geneTest, 'result3.csv')

    print("Process finished")
