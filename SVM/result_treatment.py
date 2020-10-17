'''
Jean-Clément GALLARDO
13/05/2020
Result treatment for Support Vector Machine Model
'''
import numpy as np
import pandas as pd
import os
import ast
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

#Function
##Import
def import_csv(path):
    ## Import
    script_dir = os.path.dirname("multinomial_logistic_regression")
    data_path = os.path.join(script_dir, path)
    file = pd.read_csv(data_path)
    if("result" in path):
        data = post_treatment_data(file)
    if("data" in path):
        file.drop(["patient_id"], axis = 1, inplace = True)
        data = file.dropna()
    return data

def post_treatment_data(file):
    data = {"sparsity" : float(file.sparsity), "score" : float(file.score),
            "coefficient" : ast.literal_eval(file.coefficient[0]), "balanced_score" : float(file.balanced_score),
            "confusion_matrix" : list(file.confusion_matrix)}
    return data

def import_file(path):
    script_dir = os.path.dirname("multinomial_logistic_regression")
    pam50_path = os.path.join(script_dir, path)
    f = open(pam50_path, "r")
    pam50 = f.read()
    f.close()
    list_pam50 = pam50.split("\n")
    return list_pam50

#Show results
def min_range(list1, list2):
    min1 = min(list1)
    min2 = min(list2)
    if(min1 < min2):
        min_range = min1
    if(min1 > min2):
        min_range = min2
    return min_range

def max_range(list1, list2):
    max1 = max(list1)
    max2 = max(list2)
    if (max1 < max2):
        max_range = max2
    if (max1 > max2):
        max_range = max1
    return max_range

def histo_force(force_coef):
    pam50 = []
    others = []

    for selector in force_coef.keys():
        if (selector == "Pam50"):
            for gene in force_coef[selector].keys():
                pam50.append(force_coef[selector][gene])
        if (selector == "others"):
            for gene in force_coef[selector].keys():
                others.append(force_coef[selector][gene])

    lim_histo = [min_range(pam50, others), max_range(pam50, others)]




    plt.figure(figsize = (15,15))
    plt.subplot(1,2,1)
    plt.hist(pam50, color = "red", range = (lim_histo[0], lim_histo[1]))
    plt.gca().set_yticklabels(['{:.0f}%'.format(x*100/len(pam50)) for x in plt.gca().get_yticks()])
    plt.xlabel('Pam50 coefs')
    plt.subplot(1,2,2)
    plt.hist(others, color = "blue", range = (lim_histo[0], lim_histo[1]))
    plt.gca().set_yticklabels(['{:.0f}%'.format(x*100/len(others)) for x in plt.gca().get_yticks()])
    plt.xlabel('others coefs')


    plt.show()

def boxplt(force_coef):
    pam50 = []
    others = []

    for selector in force_coef.keys():
        if (selector == "Pam50"):
            for gene in force_coef[selector].keys():
                pam50.append(force_coef[selector][gene])
        if (selector == "others"):
            for gene in force_coef[selector].keys():
                others.append(force_coef[selector][gene])


    plt.figure(figsize = (15,15))
    plt.boxplot([pam50, others])
    plt.gca().xaxis.set_ticklabels(['Pam50', 'others'])

    plt.show()



def ratio_subtype(data):
    filter_her2 = data["subtype"]=="Her2"
    Her2 = data.where(filter_her2)
    Her2 = Her2.dropna()

    filter_basal = data["subtype"]=="Basal"
    Basal = data.where(filter_basal)
    Basal = Basal.dropna()

    filter_normal = data["subtype"]=="Normal"
    Normal = data.where(filter_normal)
    Normal = Normal.dropna()

    filter_luma = data["subtype"]=="LumA"
    LumA = data.where(filter_luma)
    LumA = LumA.dropna()

    filter_lumb = data["subtype"]=="LumB"
    LumB = data.where(filter_lumb)
    LumB = LumB.dropna()

    print("Nb subtype")
    print("Her2 : ", len(Her2))
    print("Basal : ", len(Basal))
    print("Normal : ", len(Normal))
    print("LumA : ", len(LumA))
    print("LumB : ", len(LumB))

def show_results(list_pam50, result, data):
    sparsity = result["sparsity"]
    score = result["score"]
    genes = result["coefficient"].keys()
    matrix = result["confusion_matrix"]
    balanced_score = result["balanced_score"]
    force_coef = {"Pam50" : {}, "others" : {}}

    for selector in result.keys():
        if (selector == "coefficient"):
            for gene in result[selector].keys():
                if (gene in list_pam50 and result[selector][gene] != 0):
                    force_coef['Pam50'][gene] = result[selector][gene]
                if (gene not in list_pam50 and result[selector][gene] != 0):
                    force_coef['others'][gene] = result[selector][gene]
        else :
            continue


    print("Sparsity : %.2f%%" % sparsity)
    print("Test score : %.4f" % score)
    print("balanced score : %.4f" % balanced_score)
    print("Confusion matrix : \n",matrix[0])

    number = len(force_coef["Pam50"]) + len(force_coef["others"])
    print("Number of genes used for subtyping : ",number)
    print("With : ",len(force_coef["Pam50"])," genes from pam50 list")
    print("\n")

    ratio_subtype(data)
    histo_force(force_coef)
    boxplt(force_coef)


if __name__ == "__main__":
    result = import_csv("result3.csv")
    data = import_csv("../data/data.csv")
    list_pam50 = import_file("../data/PAM50_geneIDs.txt")
    show_results(list_pam50, result, data)
