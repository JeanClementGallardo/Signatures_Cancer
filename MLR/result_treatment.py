'''
Jean-Clément GALLARDO
09/04/2020
Result treatment
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
    if("test" in path):
        data = post_treatment_data(file)
    if("data" in path):
        file.drop(["patient_id"], axis = 1, inplace = True)
        data = file.dropna()
    return data

def post_treatment_data(file):
    data = {"sparsity" : float(file.sparsity), "score" : float(file.score),
            "balanced_score" : float(file.balanced_score), "confusion_matrix" : list(file.confusion_matrix),
            "Normal" : ast.literal_eval(file.Normal[0]),
            "Basal" : ast.literal_eval(file.Basal[0]),
            "Her2" : ast.literal_eval(file.Her2[0]),
            "LumA" : ast.literal_eval(file.LumA[0]),
            "LumB" : ast.literal_eval(file.LumB[0])}
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

def compare_lim(list_min, list_max):
    list_lim = []
    tmp_min = 0
    tmp_max = 0
    for idx in range(len(list_min)):
        if (tmp_min > list_min[idx]):
            tmp_min = list_min[idx]
    list_lim.append(tmp_min)
    for idx in range(len(list_max)):
        if (tmp_max < list_max[idx]):
            tmp_max = list_max[idx]
    list_lim.append(tmp_max)

    return list_lim

def all_pam50(genes,force_coef,list_pam50):
    pam50_found = []

    for subtype in force_coef.keys():
        potential_gene = []
        for pam in force_coef[subtype]['Pam50']:
            if(force_coef[subtype]['Pam50'][pam] != 0):
                potential_gene.append(pam)
        for gene in list_pam50:
            if(gene in potential_gene):
                if (gene not in pam50_found):
                    pam50_found.append(gene)
    print(pam50_found, len(pam50_found))

def compare_genes(force_coef):
    list_gene = []

    for subtype in force_coef.keys():
        for genes in force_coef[subtype]["Pam50"].keys():
            if (genes not in list_gene):
                list_gene.append(genes)
        for genes in force_coef[subtype]["others"].keys():
            if (genes not in list_gene):
                list_gene.append(genes)

    print("Total of different genes used in prediction : ", len(list_gene))

def histo_force(force_coef):
    Basal_pam50 = []
    Basal_others = []
    Her2_pam50 = []
    Her2_others = []
    LumA_pam50 = []
    LumA_others = []
    LumB_pam50 = []
    LumB_others = []
    Normal_pam50 = []
    Normal_others = []

    for subtype in force_coef.keys():
        if (subtype == "Basal"):
            for pam in force_coef[subtype]["Pam50"].keys():
                Basal_pam50.append(force_coef[subtype]["Pam50"][pam])
            for pam in force_coef[subtype]["others"].keys():
                Basal_others.append(force_coef[subtype]["others"][pam])
        if (subtype == "Normal"):
            for pam in force_coef[subtype]["Pam50"].keys():
                Normal_pam50.append(force_coef[subtype]["Pam50"][pam])
            for pam in force_coef[subtype]["others"].keys():
                Normal_others.append(force_coef[subtype]["others"][pam])
        if (subtype == "Her2"):
            for pam in force_coef[subtype]["Pam50"].keys():
                Her2_pam50.append(force_coef[subtype]["Pam50"][pam])
            for pam in force_coef[subtype]["others"].keys():
                Her2_others.append(force_coef[subtype]["others"][pam])
        if (subtype == "LumA"):
            for pam in force_coef[subtype]["Pam50"].keys():
                LumA_pam50.append(force_coef[subtype]["Pam50"][pam])
            for pam in force_coef[subtype]["others"].keys():
                LumA_others.append(force_coef[subtype]["others"][pam])
        if (subtype == "LumB"):
            for pam in force_coef[subtype]["Pam50"].keys():
                LumB_pam50.append(force_coef[subtype]["Pam50"][pam])
            for pam in force_coef[subtype]["others"].keys():
                LumB_others.append(force_coef[subtype]["others"][pam])

    list_max = []
    list_min = []

    list_max.append(max_range(Basal_pam50, Basal_others))
    list_min.append(min_range(Basal_pam50, Basal_others))

    list_max.append(max_range(Normal_pam50, Normal_others))
    list_min.append(min_range(Normal_pam50, Normal_others))

    list_max.append(max_range(Her2_pam50, Her2_others))
    list_min.append(min_range(Her2_pam50, Her2_others))

    list_max.append(max_range(LumA_pam50, LumA_others))
    list_min.append(min_range(LumA_pam50, LumA_others))

    list_max.append(max_range(LumB_pam50, LumB_others))
    list_min.append(min_range(LumB_pam50, LumB_others))

    lim_histo = compare_lim(list_min, list_max)




    plt.figure(figsize = (15,15))
    plt.subplot(2,5,1)
    plt.hist(Normal_pam50, color = "red", range = (lim_histo[0], lim_histo[1]))
    plt.gca().set_yticklabels(['{:.0f}%'.format(x*100/len(Normal_pam50)) for x in plt.gca().get_yticks()])
    plt.xlabel('Normal subtype : Pam50 coefs')
    plt.subplot(2,5,2)
    plt.hist(Her2_pam50, color = "red", range = (lim_histo[0], lim_histo[1]))
    plt.gca().set_yticklabels(['{:.0f}%'.format(x*100/len(Her2_pam50)) for x in plt.gca().get_yticks()])
    plt.xlabel('Her2 subtype : Pam50 coefs')
    plt.subplot(2,5,3)
    plt.hist(Basal_pam50, color = "red", range = (lim_histo[0], lim_histo[1]))
    plt.gca().set_yticklabels(['{:.0f}%'.format(x*100/len(Basal_pam50)) for x in plt.gca().get_yticks()])
    plt.xlabel('Basal subtype : Pam50 coefs')
    plt.subplot(2,5,4)
    plt.hist(LumA_pam50, color = "red", range = (lim_histo[0], lim_histo[1]))
    plt.gca().set_yticklabels(['{:.0f}%'.format(x*100/len(LumA_pam50)) for x in plt.gca().get_yticks()])
    plt.xlabel('LumA subtype : Pam50 coefs')
    plt.subplot(2,5,5)
    plt.hist(LumB_pam50, color = "red", range = (lim_histo[0], lim_histo[1]))
    plt.gca().set_yticklabels(['{:.0f}%'.format(x*100/len(LumB_pam50)) for x in plt.gca().get_yticks()])
    plt.xlabel('LumB subtype : Pam50 coefs')
    plt.subplot(2,5,6)
    plt.hist(Normal_others, color = "blue", range = (lim_histo[0], lim_histo[1]))
    plt.gca().set_yticklabels(['{:.0f}%'.format(x*100/len(Normal_others)) for x in plt.gca().get_yticks()])
    plt.xlabel('Normal subtype : others coefs')
    plt.subplot(2,5,7)
    plt.hist(Her2_others, color = "blue", range = (lim_histo[0], lim_histo[1]))
    plt.gca().set_yticklabels(['{:.0f}%'.format(x*100/len(Her2_others)) for x in plt.gca().get_yticks()])
    plt.xlabel('Her2 subtype : others coefs')
    plt.subplot(2,5,8)
    plt.hist(Basal_others, color = "blue", range = (lim_histo[0], lim_histo[1]))
    plt.gca().set_yticklabels(['{:.0f}%'.format(x*100/len(Basal_others)) for x in plt.gca().get_yticks()])
    plt.xlabel('Basal subtype : others coefs')
    plt.subplot(2,5,9)
    plt.hist(LumA_others, color = "blue", range = (lim_histo[0], lim_histo[1]))
    plt.gca().set_yticklabels(['{:.0f}%'.format(x*100/len(LumA_others)) for x in plt.gca().get_yticks()])
    plt.xlabel('LumA subtype : others coefs')
    plt.subplot(2,5,10)
    plt.hist(LumB_others, color = "blue", range = (lim_histo[0], lim_histo[1]))
    plt.gca().set_yticklabels(['{:.0f}%'.format(x*100/len(LumB_others)) for x in plt.gca().get_yticks()])
    plt.xlabel('LumB subtype : others coefs')

    plt.show()

def boxplt(force_coef):
    Basal_pam50 = []
    Basal_others = []
    Her2_pam50 = []
    Her2_others = []
    LumA_pam50 = []
    LumA_others = []
    LumB_pam50 = []
    LumB_others = []
    Normal_pam50 = []
    Normal_others = []

    for subtype in force_coef.keys():
        if (subtype == "Basal"):
            for pam in force_coef[subtype]["Pam50"].keys():
                Basal_pam50.append(force_coef[subtype]["Pam50"][pam])
            for pam in force_coef[subtype]["others"].keys():
                Basal_others.append(force_coef[subtype]["others"][pam])
        if (subtype == "Normal"):
            for pam in force_coef[subtype]["Pam50"].keys():
                Normal_pam50.append(force_coef[subtype]["Pam50"][pam])
            for pam in force_coef[subtype]["others"].keys():
                Normal_others.append(force_coef[subtype]["others"][pam])
        if (subtype == "Her2"):
            for pam in force_coef[subtype]["Pam50"].keys():
                Her2_pam50.append(force_coef[subtype]["Pam50"][pam])
            for pam in force_coef[subtype]["others"].keys():
                Her2_others.append(force_coef[subtype]["others"][pam])
        if (subtype == "LumA"):
            for pam in force_coef[subtype]["Pam50"].keys():
                LumA_pam50.append(force_coef[subtype]["Pam50"][pam])
            for pam in force_coef[subtype]["others"].keys():
                LumA_others.append(force_coef[subtype]["others"][pam])
        if (subtype == "LumB"):
            for pam in force_coef[subtype]["Pam50"].keys():
                LumB_pam50.append(force_coef[subtype]["Pam50"][pam])
            for pam in force_coef[subtype]["others"].keys():
                LumB_others.append(force_coef[subtype]["others"][pam])

    plt.figure(figsize = (15,15))
    plt.subplot(2,5,1)
    plt.boxplot([Normal_pam50, Normal_others])
    plt.gca().xaxis.set_ticklabels(['Pam50', 'others'])
    plt.xlabel('Normal subtype')
    plt.subplot(2,5,2)
    plt.boxplot([Her2_pam50, Her2_others])
    plt.gca().xaxis.set_ticklabels(['Pam50', 'others'])
    plt.xlabel('Her2 subtype')
    plt.subplot(2,5,3)
    plt.boxplot([Basal_pam50, Basal_others])
    plt.gca().xaxis.set_ticklabels(['Pam50', 'others'])
    plt.xlabel('Basal subtype')
    plt.subplot(2,5,4)
    plt.boxplot([LumA_pam50, LumA_others])
    plt.gca().xaxis.set_ticklabels(['Pam50', 'others'])
    plt.xlabel('LumA subtype')
    plt.subplot(2,5,5)
    plt.boxplot([LumB_pam50, LumB_others])
    plt.gca().xaxis.set_ticklabels(['Pam50', 'others'])
    plt.xlabel('LumB subtype')

    plt.show()

def high_force(force_coef, genes,list_pam50, value_pos, value_neg):
    high_coef = {  'Basal':{'Pam50' : {}, 'others' : {}},
                    'Her2':{'Pam50' : {}, 'others' : {}},
                    'LumA':{'Pam50' : {}, 'others' : {}},
                    'LumB':{'Pam50' : {}, 'others' : {}},
                    'Normal':{'Pam50' : {}, 'others' : {}}}

    for subtype in force_coef.keys():
        for gene in force_coef[subtype]["Pam50"].keys():
            if (force_coef[subtype]['Pam50'][gene] >= value_pos or force_coef[subtype]['Pam50'][gene] <= value_neg):
                high_coef[subtype]["Pam50"][gene] = force_coef[subtype]['Pam50'][gene]
        for gene in force_coef[subtype]["others"].keys():
            if (force_coef[subtype]['others'][gene] >= value_pos or force_coef[subtype]['others'][gene] <= value_neg):
                high_coef[subtype]["others"][gene] = force_coef[subtype]['others'][gene]

    return high_coef

def print_dict(high_coef):
    for subtype in high_coef.keys():
        print("\n")
        print("High genes coef for subtype %s" % subtype)
        for label in high_coef[subtype].keys():
            print(label)
            for genes in high_coef[subtype][label].keys():
                print(genes," : ",high_coef[subtype][label][genes])
    print("\n")

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
    balanced_score = result["balanced_score"]
    matrix = result["confusion_matrix"]
    genes = result["Normal"].keys()
    force_coef = {  'Basal':{'Pam50' : {}, 'others' : {}},
                    'Her2':{'Pam50' : {}, 'others' : {}},
                    'LumA':{'Pam50' : {}, 'others' : {}},
                    'LumB':{'Pam50' : {}, 'others' : {}},
                    'Normal':{'Pam50' : {}, 'others' : {}}}

    list_subtype = ["Her2","Basal","Normal","LumA","LumB"]

    for selector in result.keys():
        if (selector not in list_subtype):
            continue
        else :
            for gene in result[selector].keys():
                if (gene in list_pam50 and result[selector][gene] != 0):
                    force_coef[selector]['Pam50'][gene] = result[selector][gene]
                if (gene not in list_pam50 and result[selector][gene] != 0):
                    force_coef[selector]['others'][gene] = result[selector][gene]


    print("Sparsity with L1 penalty: %.2f%%" % sparsity)
    print("Test score with L1 penalty: %.4f" % score)
    print("Test balanced score with L1 penalty: %.4f" % balanced_score)
    print("Confusion matrix : \n",matrix[0])

    for subtype in force_coef.keys():
        number = len(force_coef[subtype]["Pam50"]) + len(force_coef[subtype]["others"])
        print("Number of genes used for subtyping ",subtype," : ",number)
        print("With : ",len(force_coef[subtype]["Pam50"])," genes from pam50 list")
        print("\n")

    print("Genes found with pam50 label for all subtypes")
    all_pam50(genes,force_coef,list_pam50)
    compare_genes(force_coef)
    high_coef = high_force(force_coef, genes, list_pam50, 0.004, -0.004)
    # print_dict(high_coef)
    # ratio_subtype(data)
    histo_force(force_coef)
    boxplt(force_coef)


if __name__ == "__main__":
    result = import_csv("test.csv")
    data = import_csv("../data/data_pam50.csv")
    list_pam50 = import_file("../data/PAM50_geneIDs.txt")
    show_results(list_pam50, result, data)
