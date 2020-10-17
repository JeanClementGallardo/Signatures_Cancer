'''
Jean-Clément GALLARDO
09/04/2020
Logistic regression model creation
'''


import numpy as np
import pandas as pd
import os
import csv

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

# Export results
def export_result(model, data, geneTest,filename):
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



if __name__ == "__main__":
    encoder = LabelEncoder()

    ## Import data
    data = import_csv("../data/data.csv", encoder)
    print(data)
    # selected_gene = ["subtype","PHGDH","KIF2C","CDC6","UBE2T","MMP11","NUF2","MLPH","CENPF","NDC80","CCNE1","BAG1"]
    # new_data = data[selected_gene]

    ## Split data into training and test dataset
    geneTrain, geneTest= train_test_split(data, train_size = 116, random_state = 69780)
    print("data split")
    ## Create and fit models
    params_rank1 = [0.5,0.001]
    # params_rank2 = [0.5,0.01]
    # params_rank3 = [0.1,0.001]
    # params_rank4 = [0.1,0.01]
    # params_rank5 = [0.5,0.025]
    # params_rank6 = [0.5,0.05]

    logistic_rank1 = create_model(geneTrain, params_rank1[0], params_rank1[1])
    # logistic_rank2 = create_model(geneTrain, params_rank2[0], params_rank2[1])
    # logistic_rank3 = create_model(geneTrain, params_rank3[0], params_rank3[1])
    # logistic_rank4 = create_model(geneTrain, params_rank4[0], params_rank4[1])
    # logistic_rank5 = create_model(geneTrain, params_rank5[0], params_rank5[1])
    # logistic_rank6 = create_model(geneTrain, params_rank6[0], params_rank6[1])
    print("model created")

    ## Export results
    export_result(logistic_rank1, data, geneTest, 'article_result.csv')
    # export_result(logistic_rank2, data, geneTest, 'result2.csv')
    # export_result(logistic_rank3, data, geneTest, 'result3.csv')
    # export_result(logistic_rank4, data, geneTest, 'result4.csv')
    # export_result(logistic_rank5, data, geneTest, 'result5.csv')
    # export_result(logistic_rank6, data, geneTest, 'result6.csv')
    print("Process finished")
