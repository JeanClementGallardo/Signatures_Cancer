import numpy as np
import pandas as pd
import os
import csv

from sklearn.linear_model import Lasso
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
def create_model(geneTrain, alpha_value):
    lasso_model = Lasso(
        alpha = alpha_value, fit_intercept = True, normalize = True, max_iter = 5000
    )
    lasso_model.fit(geneTrain.loc[:, geneTrain.columns != 'subtype'], geneTrain.subtype)
    return lasso_model

# Export results
def export_result(model, data, geneTest,filename):
    sparsity = np.mean(model.coef_ == 0) * 100
    score = model.score(geneTest.loc[:, geneTest.columns != 'subtype'], geneTest.subtype)
    coefficient = model.coef_
    genes = data.columns[1:]
    export_dict = {"sparsity" : sparsity, "score" : score,
                    "coefficient" : {}}

    for i in range(len(coefficient)):
        export_dict["coefficient"][genes[i]] = coefficient[i]

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
    lasso_model = create_model(geneTrain, 0.001)
    print("model created")

    ## Export results
    export_result(lasso_model, data, geneTest, 'result.csv')
    print("Process finished")
