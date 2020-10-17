'''
Jean-Clement GALLARDO
27/04/2020
Create csv file to contain the clinical data
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

def import_file(path):
    f = open(path, "r")
    pam50 = f.read()
    f.close()
    list_pam50 = pam50.split("\n")
    return list_pam50

#Pre-process Genomics and miRNA data
def process_expr_data(data,list_pam50):
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
                if (rows[0] in list_pam50):
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
def create_fieldnames(clinical_data, clinical_id):
    idx = 0
    list_fieldnames = []

    for idx in range(len(clinical_id)):
            list_fieldnames.append(clinical_id[idx])

    return list_fieldnames

#Create value list index on fieldnames list
def create_value_list(gene_data,clinical_data, clinical_id):
    value_list = []
    value_list.append(create_fieldnames(clinical_data, clinical_id))
    tmp_list = []
    for patient in gene_data.keys():
        for clinical in clinical_data[patient].keys():
            if (clinical in clinical_id):
                tmp_list.append(clinical_data[patient][clinical])
        value_list.append(tmp_list)
        tmp_list = []

    return value_list


if __name__ == "__main__":
    list_pam50 = import_file("../PAM50_geneIDs.txt")

    expr_gene = open_file("../HiSeqV2.txt","\n")
    gene_data = process_expr_data(expr_gene,list_pam50)

    clinical = open_file("../BRCA_clinicalMatrix.txt","\n")
    clinical_data = process_clinical_data(clinical)

    clinical_id = ["PAM50Call_RNAseq",
                "breast_cancer_surgery_margin_status",
                "breast_carcinoma_estrogen_receptor_status",
                "breast_carcinoma_immunohistochemistry_er_pos_finding_scale",
                "breast_carcinoma_immunohistochemistry_pos_cell_score",
                "breast_carcinoma_immunohistochemistry_prgstrn_rcptr_ps_fndng_scl",
                "breast_carcinoma_primary_surgical_procedure_name",
                "breast_carcinoma_progesterone_receptor_status",
                "breast_carcinoma_surgical_procedure_name",
                "cytokeratin_immunohistochemistry_staining_method_mcrmtstss_ndctr",
                "distant_metastasis_present_ind2",
                "her2_and_centromere_17_positive_finding_other_measuremnt_scl_txt",
                "her2_erbb_method_calculation_method_text",
                "her2_erbb_pos_finding_cell_percent_category",
                "her2_erbb_pos_finding_fluorescence_n_st_hybrdztn_clcltn_mthd_txt",
                "her2_immunohistochemistry_level_result",
                "her2_neu_and_centromere_17_copy_number_analysis_npt_ttl_nmbr_cnt",
                "her2_neu_breast_carcinoma_copy_analysis_input_total_number",
                "her2_neu_chromosone_17_signal_ratio_value",
                "her2_neu_metastatic_breast_carcinoma_copy_analysis_inpt_ttl_nmbr",
                "histological_type",
                "immunohistochemistry_positive_cell_score",
                "metastatic_breast_carcinm_ps_fndng_prgstrn_rcptr_thr_msr_scl_txt",
                "metastatic_breast_carcinom_lb_prc_hr2_n_mmnhstchmstry_rcptr_stts",
                "metastatic_breast_carcinoma_erbb2_immunohistochemistry_levl_rslt",
                "metastatic_breast_carcinoma_estrogen_receptor_detection_mthd_txt",
                "metastatic_breast_carcinoma_estrogen_receptor_status",
                "metastatic_breast_carcinoma_estrogen_receptr_lvl_cll_prcnt_ctgry",
                "metastatic_breast_carcinoma_her2_erbb_method_calculatin_mthd_txt",
                "metastatic_breast_carcinoma_her2_erbb_pos_findng_cll_prcnt_ctgry",
                "metastatic_breast_carcinoma_her2_neu_chromosone_17_signal_rat_vl",
                "metastatic_breast_carcinoma_immunhstchmstry_r_pstv_fndng_scl_typ",
                "metastatic_breast_carcinoma_immunohistochemistry_er_pos_cell_scr",
                "metastatic_breast_carcinoma_immunohistochemistry_pr_pos_cell_scr",
                "metastatic_breast_carcinoma_lab_proc_hr2_n_n_st_hybrdztn_tcm_typ",
                "metastatic_breast_carcinoma_pos_finding_hr2_rbb2_thr_msr_scl_txt",
                "metastatic_breast_carcinoma_progestern_rcptr_lvl_cll_prcnt_ctgry",
                "metastatic_breast_carcinoma_progesterone_receptor_dtctn_mthd_txt",
                "metastatic_breast_carcinoma_progesterone_receptor_status",
                "metastatic_site_at_diagnosis",
                "menopause_status",
                "pathologic_M",
                "pathologic_N",
                "pathologic_T",
                "pathologic_stage",
                "pos_finding_her2_erbb2_other_measurement_scale_text",
                "pos_finding_metastatic_brst_crcnm_strgn_rcptr_thr_msrmnt_scl_txt",
                "pos_finding_progesterone_receptor_other_measurement_scale_text",
                "positive_finding_estrogen_receptor_other_measurement_scale_text",
                "postoperative_rx_tx",
                "progesterone_receptor_level_cell_percent_category"]

    value_list = create_value_list(gene_data, clinical_data, clinical_id)

    create_csv(value_list, '../test.csv')

    print("Csv ready to use")
