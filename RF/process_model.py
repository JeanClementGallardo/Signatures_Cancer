'''
Jean-Clément GALLARDO
15/04/2020
Logistic regression model creation
'''


import numpy as np
import pandas as pd
import os
import csv

from sklearn.ensemble import RandomForestClassifier
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

# Create and fit Random Forest model
def create_model(geneTrain):
    random_forest = RandomForestClassifier(
        criterion = 'entropy', n_estimators = 1000
    )
    random_forest.fit(geneTrain.loc[:, geneTrain.columns != 'subtype'], geneTrain.subtype)
    return random_forest

# Export results
def export_result(model, data, geneTest):
    sparsity = np.mean(model.feature_importances_ == 0) * 100
    score = model.score(geneTest.loc[:, geneTest.columns != 'subtype'], geneTest.subtype)
    subtypes = encoder.inverse_transform(model.classes_)
    balanced_score = balanced_accuracy_score(geneTest.subtype, model.predict(geneTest.loc[:, geneTest.columns != 'subtype']))
    matrix = confusion_matrix(encoder.inverse_transform(geneTest.subtype),
                              encoder.inverse_transform(model.predict(geneTest.loc[:, geneTest.columns != 'subtype'])),
                              labels = subtypes)
    coefficient = model.feature_importances_
    # decision = model.decision_path(geneTest.loc[:, geneTest.columns != 'subtype'])
    # leaves = model.apply(geneTest.loc[:, geneTest.columns != 'subtype'])
    genes = data.columns[1:]
    export_dict = {"sparsity" : sparsity, "score" : score,
                   "balanced_score" : balanced_score, "confusion_matrix" : matrix ,
                   "coefficient" : {}}

    for idx in range(len(coefficient)):
        export_dict["coefficient"][genes[idx]] = coefficient[idx]

    with open('result.csv', 'w') as csvfile :
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
    random_rank1 = create_model(geneTrain)
    print("model created")

    ## Export results
    export_result(random_rank1, data, geneTest)
    print("Process finished")
