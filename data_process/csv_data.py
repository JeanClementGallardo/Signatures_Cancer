'''
Jean-Clement GALLARDO
16/03/2020
Create csv file to contain the data
'''

#Import package
import csv

#Open file
def open_file(FileName, split_sep):
    f = open(FileName,"r")
    contents = f.read()
    data = contents.split(split_sep)
    f.close()
    return data

#Pre-process Genomics and miRNA data
def process_expr_data(data):
    list_sample = []
    expr_dict = {}

    for block in range(len(data)):
        rows = data[block].split("\t")
        if(block == 0):
            for row in range(len(rows)):
                if (rows[row] == "sample"):
                    continue
                else :
                    list_sample.append(rows[row])
                    expr_dict[rows[row]] = {}
        if(block > 0):
            for row in range(len(rows)):
                if (row == 0):
                    continue
                if (rows[row] == "NA"):
                    expr_dict[list_sample[row-1]][rows[0]] = rows[row]
                else :
                    expr_dict[list_sample[row-1]][rows[0]] = float(rows[row])
    return expr_dict

#Pre-process clinical data
def process_clinical_data(data):
    list_identifier = []
    clinical_data = {}

    for block in range(len(data)):
        rows = data[block].split("\t")
        if (block == 0):
            for row in range(len(rows)):
                list_identifier.append(rows[row])
        if (block > 0):
            for row in range(len(rows)):
                if (row == 0):
                    clinical_data[rows[row]] = {}
                else :
                    clinical_data[rows[0]][list_identifier[row]] = rows[row]

    return clinical_data



#Write in csv file
def create_csv(value_list, FileName):
    with open(FileName, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerows(value_list)


#Create fieldnames list
def create_fieldnames(gene_data):
    idx = 0
    list_fieldnames = []
    list_fieldnames.append('patient_id')
    list_fieldnames.append('subtype')
    for patient in gene_data.keys():
        for gene in gene_data[patient].keys():
            if (idx == 0):
                list_fieldnames.append(gene)
            if (idx == 1):
                break
        idx+=1

    return list_fieldnames

#Create value list index on fieldnames list
def create_value_list(gene_data,clinical_data):
    value_list = []
    value_list.append(create_fieldnames(gene_data))
    tmp_list = []
    for patient in gene_data.keys():
        tmp_list.append(patient)
        tmp_list.append(clinical_data[patient]["PAM50Call_RNAseq"])
        for gene in gene_data[patient].keys():
            tmp_list.append(gene_data[patient][gene])
        value_list.append(tmp_list)
        tmp_list = []

    return value_list


if __name__ == "__main__":
    expr_gene = open_file("../HiSeqV2.txt","\n")
    gene_data = process_expr_data(expr_gene)

    clinical = open_file("../BRCA_clinicalMatrix.txt","\n")
    clinical_data = process_clinical_data(clinical)

    value_list = create_value_list(gene_data, clinical_data)

    create_csv(value_list, '../data.csv')

    print("Csv ready to use")
