#!/usr/bin/python
from __future__ import print_function

from RMSbufferingMethods import dataAnalysisMethods


import numpy
import mysql.connector
import itertools
import collections


####################################
# GLOBAL VARIABLES

host = ""
user = ""
password = ""
database = ""
directory = ""
organism = ""
dataType = ""
enzymeClassification = ""
height = ""


def set_global_variables(h, u, p, da, di, org, dt, enzC, ht):
    global host
    global user
    global password
    global database
    global directory
    global organism
    global dataType
    global enzymeClassification
    global height
    host = h
    user = u
    password = p
    database = da
    directory = di
    organism = org
    dataType = dt
    enzymeClassification = enzC
    height = ht

    return

######################################################################


def retrieve_pair_enzyme_data_for_tests_pombe(*args, **kwargs):
    # connection to the DB
    global host
    global user
    global password
    global database
    group_type = kwargs.get('enzymeClassification', "RMS")
    with_backup_only = kwargs.get('with_backup_only', False)

    pair_var_cor = {}
    # pairs_no_positive_interaction = {}
    # pairs_no_negative_interaction = {}
    # enzyme_list = []
    enzyme_sd = {}  # key = enzyme id, value = [ 0:SD_N, 1:SD_H, 2:delta_SD, 3: mean of the fold change, 4: SD of FC]

    if group_type == "RMS":
        height = kwargs.get('height', None)
        # Retrieve counts data
        # NORMAL CONDITION
        enzyme_counts_N = retrieve_count_data(group_type, "N", height=height)
        # OXY CONDITION
        enzyme_counts_H = retrieve_count_data(group_type, "H", height=height)

        # classification retrieval
        g_enz_cpd, enz_g_cpd = retrieve_enzyme_group_cpd(group_type, height)
        # retrieve enzymes to study (with RMS H2 OR EC number AND not alone in their group)
        enzyme_list = retrieve_enzymes_to_study(group_type, height, with_backup_only)

    elif group_type == "EC":
        height = -1
        # Retrieve counts data
        # NORMAL CONDITION
        enzyme_counts_N = retrieve_count_data(group_type, "N")
        # OXY CONDITION
        enzyme_counts_H = retrieve_count_data(group_type, "H")

        # classification retrieval
        g_enz_cpd, enz_g_cpd = retrieve_enzyme_group_cpd("EC", height)
        # retrieve enzymes to study (with RMS H2 OR EC number AND not alone in their group)
        enzyme_list = retrieve_enzymes_to_study(group_type, height, with_backup_only)

    elif group_type == "all":
        height = kwargs.get('height', 2)
        # Retrieve counts data

        # NORMAL CONDITION
        enzyme_counts_N = retrieve_count_data(group_type, "N")
        # OXY CONDITION
        enzyme_counts_H = retrieve_count_data(group_type, "H")

        # classification retrieval
        g_enz_cpd, enz_g_cpd = retrieve_enzyme_group_cpd(group_type, height)

        # retrieve enzymes to study (with RMS H2 OR EC number AND not alone in their group)
        enzyme_list = retrieve_enzymes_to_study(group_type, height, with_backup_only)

    else:
        print("Unrecognized enzyme classification :/")
        return

    # negative_gen_pairs = get_genetic_interactions("Negative")
    # positive_gen_pairs = get_genetic_interactions("Positive")

    # analysis for individual enzymes
    for enz in enzyme_list:
        if enz in enzyme_counts_H.keys() and enz in enzyme_counts_H.keys():
            enzyme_sd[enz] = [0, 0, 0, 0, 0]  # 0:SD_N, 1:SD_H, 2:delta_SD, 3: mean of the fold change, 4: SD of FC
            enzyme_sd[enz][0] = numpy.std(numpy.log2(list(enzyme_counts_N[enz].values())), dtype=numpy.float64)  # SD in N condition
            enzyme_sd[enz][1] = numpy.std(numpy.log2(list(enzyme_counts_H[enz].values())), dtype=numpy.float64)  # SD in OXY condition
            enzyme_sd[enz][2] = numpy.absolute(
                enzyme_sd[enz][0] - enzyme_sd[enz][1])  # DELTA variance between two conditions

            fcs = {}  # fold change for each strain
            for strain in enzyme_counts_H[enz].keys():
                if strain in enzyme_counts_N[enz].keys() and enzyme_counts_N[enz][strain] != 0:
                    fcs[strain] = numpy.log2(enzyme_counts_H[enz][strain]) / numpy.log2(enzyme_counts_N[enz][strain])

            enzyme_sd[enz][3] = numpy.mean([fcs[e] for e in fcs.keys() if fcs[e] != 0], dtype=numpy.float64)  # mean of the fold change
            enzyme_sd[enz][4] = numpy.std([fcs[e] for e in fcs.keys() if fcs[e] != 0], dtype=numpy.float64)  # SD of the fold change
        else:
            # remove enzyme from enzymeList
            enzyme_list.remove(enz)

    # analysis for enzyme pairs
    for pair in itertools.combinations(enzyme_list, 2):
        gmean_sum_H, pair_sd_H, pair_correlation_H = dataAnalysisMethods.analyse_pair(pair, enzyme_counts_H)
        gmean_sum_N, pair_sd_N, pair_correlation_N = dataAnalysisMethods.analyse_pair(pair, enzyme_counts_N)

        delta_sd = pair_sd_N - pair_sd_H

        sd_fc_pair, sd_reduction_fc_min, mean_fc_pair = dataAnalysisMethods.analyse_pair_fold_change(pair, enzyme_counts_N, enzyme_counts_H)

        sd_reduction_score_N_min = pair_sd_N - min([enzyme_sd[pair[0]][0], enzyme_sd[pair[1]][0]])
        sd_reduction_score_H_min = pair_sd_H - min([enzyme_sd[pair[0]][1], enzyme_sd[pair[1]][1]])

        same_group = dataAnalysisMethods.is_in_same_group(pair, enz_g_cpd)

        pair_var_cor[pair] = [gmean_sum_N, pair_sd_N, pair_correlation_N, gmean_sum_H, pair_sd_H, pair_correlation_H,
                              delta_sd, sd_reduction_score_N_min, sd_reduction_score_H_min, sd_fc_pair,
                              sd_reduction_fc_min, mean_fc_pair, same_group, "RNA"]

    return pair_var_cor  # , pairs_no_positive_interaction, pairs_no_negative_interaction

######################################################################


def retrieve_pair_enzyme_data_for_tests_pombe_prot(*args, **kwargs):
    # connection to the DB
    global host
    global user
    global password
    global database
    group_type = kwargs.get('enzymeClassification', "RMS")

    pair_var_cor = {}
    pairs_no_positive_interaction = {}
    pairs_no_negative_interaction = {}
    enzyme_list = []
    enzyme_sd = {}  # key = enzyme id, value = [ 0:SD_N, 1:SD_H, 2:delta_SD, 3: mean of the fold change, 4: SD of FC]

    if group_type == "RMS":
        height = kwargs.get('height', None)
        # Retrieve counts data
        # NORMAL CONDITION
        enzyme_counts_N = retrieve_protein_count_data(group_type, "N", height=height)
        # OXY CONDITION
        enzyme_counts_H = retrieve_protein_count_data(group_type, "H", height=height)

        # classification retrieval
        g_enz_cpd, enz_g_cpd = retrieve_enzyme_group_cpd_protein(group_type, height)
        # retrieve enzymes to study (with RMS H2 OR EC number AND not alone in their group)
        enzyme_list = retrieve_enzymes_to_study_protein(group_type, height)

    elif group_type == "EC":
        height = -1
        # Retrieve counts data
        # NORMAL CONDITION
        enzyme_counts_N = retrieve_protein_count_data(group_type, "N")
        # OXY CONDITION
        enzyme_counts_H = retrieve_protein_count_data(group_type, "H")

        # classification retrieval
        g_enz_cpd, enz_g_cpd = retrieve_enzyme_group_cpd_protein("EC", height)
        # retrieve enzymes to study (with RMS H2 OR EC number AND not alone in their group)
        enzyme_list = retrieve_enzymes_to_study_protein(group_type, height)

    elif group_type == "all":
        height = kwargs.get('height', 2)
        # Retrieve counts data

        # NORMAL CONDITION
        enzyme_counts_N = retrieve_protein_count_data(group_type, "N")
        # OXY CONDITION
        enzyme_counts_H = retrieve_protein_count_data(group_type, "H")

        # classification retrieval
        g_enz_cpd, enz_g_cpd = retrieve_enzyme_group_cpd_protein(group_type, height)

        # retrieve enzymes to study (with RMS H2 OR EC number AND not alone in their group)
        enzyme_list = retrieve_enzymes_to_study_protein(group_type, height)

    else:
        print("Unrecognized enzyme classification :/")
        return

    # negative_gen_pairs = get_genetic_interactions("Negative")
    # positive_gen_pairs = get_genetic_interactions("Positive")

    # analysis for individual enzymes
    for enz in enzyme_list:
        if enz in enzyme_counts_H.keys() and enz in enzyme_counts_H.keys():
            enzyme_sd[enz] = [0, 0, 0, 0, 0]  # 0:SD_N, 1:SD_H, 2:delta_SD, 3: mean of the fold change, 4: SD of FC
            enzyme_sd[enz][0] = numpy.std(numpy.log2(list(enzyme_counts_N[enz].values())), dtype=numpy.float64)  # SD in N condition
            enzyme_sd[enz][1] = numpy.std(numpy.log2(list(enzyme_counts_H[enz].values())), dtype=numpy.float64)  # SD in OXY condition
            enzyme_sd[enz][2] = numpy.absolute(
                enzyme_sd[enz][0] - enzyme_sd[enz][1])  # DELTA variance between two conditions

            fcs = {}  # fold change for each strain
            for strain in enzyme_counts_H[enz].keys():
                if strain in enzyme_counts_N[enz].keys() and enzyme_counts_N[enz][strain] != 0:
                    fcs[strain] = numpy.log2(enzyme_counts_H[enz][strain]) / numpy.log2(enzyme_counts_N[enz][strain])

            enzyme_sd[enz][3] = numpy.mean([fcs[e] for e in fcs.keys() if fcs[e] != 0], dtype=numpy.float64)  # mean of the fold change
            enzyme_sd[enz][4] = numpy.std([fcs[e] for e in fcs.keys() if fcs[e] != 0], dtype=numpy.float64)  # SD of the fold change
        else:
            # remove enzyme from enzymeList
            enzyme_list.remove(enz)

    # analysis for enzyme pairs
    for pair in itertools.combinations(enzyme_list, 2):
        gmean_sum_H, pair_sd_H, pair_correlation_H = dataAnalysisMethods.analyse_pair(pair, enzyme_counts_H)
        gmean_sum_N, pair_sd_N, pair_correlation_N = dataAnalysisMethods.analyse_pair(pair, enzyme_counts_N)

        delta_sd = pair_sd_N - pair_sd_H

        sd_fc_pair, sd_reduction_fc_min, mean_fc_pair = dataAnalysisMethods.analyse_pair_fold_change(pair, enzyme_counts_N, enzyme_counts_H)

        sd_reduction_score_N_min = pair_sd_N - min([enzyme_sd[pair[0]][0], enzyme_sd[pair[1]][0]])
        sd_reduction_score_H_min = pair_sd_H - min([enzyme_sd[pair[0]][1], enzyme_sd[pair[1]][1]])

        same_group = dataAnalysisMethods.is_in_same_group(pair, enz_g_cpd)

        pair_var_cor[pair] = [gmean_sum_N, pair_sd_N, pair_correlation_N, gmean_sum_H, pair_sd_H, pair_correlation_H,
                              delta_sd, sd_reduction_score_N_min, sd_reduction_score_H_min, sd_fc_pair,
                              sd_reduction_fc_min, mean_fc_pair, same_group, "protein"]




    return pair_var_cor

######################################################################


def retrieve_pair_enzyme_data_for_tests_cerevisiae(**kwargs):
    # connection to the DB
    global host
    global user
    global password
    global database
    group_type = kwargs.get('enzymeClassification', "RMS")
    # with_backup_only = kwargs.get('withBackUpOnly', False)

    pair_var_cor = {}
    pairs_no_positive_interaction = {}
    pairs_no_negative_interaction = {}
    enzyme_list = []
    enzyme_sd = {}  # key = enzyme id, value = SD

    if group_type == "RMS":
        height = kwargs.get('height', None)
        # Retrieve counts data

        enzyme_counts = retrieve_count_data_cer(group_type, height=height)

        # classification retrieval
        g_enz_cpd, enz_g_cpd = retrieve_enzyme_group_cpd(group_type, height)
        # retrieve enzymes to study (with RMS H2 OR EC number AND not alone in their group)
        enzyme_list = retrieve_enzymes_to_study_cer(group_type, height)

    elif group_type == "EC":
        height = -1
        # Retrieve counts data
        enzyme_counts = retrieve_count_data_cer(group_type)

        # classification retrieval
        g_enz_cpd, enz_g_cpd = retrieve_enzyme_group_cpd("EC", height)
        # retrieve enzymes to study (with RMS H2 OR EC number AND not alone in their group)
        enzyme_list = retrieve_enzymes_to_study_cer(group_type, height)

    elif group_type == "all":
        height = kwargs.get('height', 2)
        # Retrieve counts data
        enzyme_counts = retrieve_count_data_cer(group_type)

        # classification retrieval
        g_enz_cpd, enz_g_cpd = retrieve_enzyme_group_cpd(group_type, height)

        # retrieve enzymes to study (with RMS H2 OR EC number AND not alone in their group)
        enzyme_list = retrieve_enzymes_to_study_cer(group_type, height)

    else:
        print("Unrecognized enzyme classification :/")
        return

    # negative_gen_pairs = get_genetic_interactions("Negative")
    # positive_gen_pairs = get_genetic_interactions("Positive")

    # analysis for individual enzymes
    for enz in enzyme_list:
        if enz in enzyme_counts.keys():
            enzyme_sd[enz] = 0
            enzyme_sd[enz] = numpy.std(numpy.log2(list(enzyme_counts[enz].values())), dtype=numpy.float64)
        else:
            # remove enzyme from enzymeList
            enzyme_list.remove(enz)

    # analysis for enzyme pairs
    for pair in itertools.combinations(enzyme_list, 2):
        gmean_sum, pair_sd, pair_correlation = dataAnalysisMethods.analyse_pair(pair, enzyme_counts)

        sd_reduction_score_min = pair_sd - min([enzyme_sd[pair[0]], enzyme_sd[pair[1]]])

        same_group = dataAnalysisMethods.is_in_same_group(pair, enz_g_cpd)

        pair_var_cor[pair] = [gmean_sum, pair_sd, pair_correlation, sd_reduction_score_min, same_group, "RNA"]

    return pair_var_cor
######################################################################


def retrieve_pair_enzyme_data_for_tests_cerevisiae_prot(**kwargs):
    # connection to the DB
    global host
    global user
    global password
    global database
    group_type = kwargs.get('enzymeClassification', "RMS")
    # with_backup_only = kwargs.get('withBackUpOnly', False)

    pair_var_cor = {}

    enzyme_list = []
    enzyme_sd = {}  # key = enzyme id, value =SD

    if group_type == "RMS":
        height = kwargs.get('height', None)
        # Retrieve counts data

        enzyme_counts = retrieve_count_data_prot_cer(group_type, height=height)

        # classification retrieval
        g_enz_cpd, enz_g_cpd = retrieve_enzyme_group_cpd_protein(group_type, height)
        # retrieve enzymes to study (with RMS H2 OR EC number AND not alone in their group)
        enzyme_list = retrieve_enzymes_to_study_cer_prot(group_type, height)

    elif group_type == "EC":
        height = -1
        # Retrieve counts data
        enzyme_counts = retrieve_count_data_prot_cer(group_type)

        # classification retrieval
        g_enz_cpd, enz_g_cpd = retrieve_enzyme_group_cpd_protein("EC", height)
        # retrieve enzymes to study (with RMS H2 OR EC number AND not alone in their group)
        enzyme_list = retrieve_enzymes_to_study_cer_prot(group_type, height)

    elif group_type == "all":
        height = kwargs.get('height', 2)
        # Retrieve counts data
        enzyme_counts = retrieve_enzymes_to_study_cer_prot(group_type)

        # classification retrieval
        g_enz_cpd, enz_g_cpd = retrieve_enzyme_group_cpd_protein(group_type, height)

        # retrieve enzymes to study (with RMS H2 OR EC number AND not alone in their group)
        enzyme_list = retrieve_enzymes_to_study_cer(group_type, height)

    else:
        print("Unrecognized enzyme classification :/")
        return

    # negative_gen_pairs = get_genetic_interactions("Negative")
    # positive_gen_pairs = get_genetic_interactions("Positive")

    # analysis for individual enzymes
    for enz in enzyme_list:
        if enz in enzyme_counts.keys():
            enzyme_sd[enz] = 0
            enzyme_sd[enz] = numpy.std(numpy.log2(list(enzyme_counts[enz].values())), dtype=numpy.float64)
        else:
            # remove enzyme from enzymeList
            enzyme_list.remove(enz)

    # analysis for enzyme pairs
    for pair in itertools.combinations(enzyme_list, 2):
        gmean_sum, pair_sd, pair_correlation = dataAnalysisMethods.analyse_pair(pair, enzyme_counts)

        sd_reduction_score_min = pair_sd - min([enzyme_sd[pair[0]], enzyme_sd[pair[1]]])

        same_group = dataAnalysisMethods.is_in_same_group(pair, enz_g_cpd)

        pair_var_cor[pair] = [gmean_sum, pair_sd, pair_correlation, sd_reduction_score_min, same_group, "protein"]

    return pair_var_cor

######################################################################


def retrieve_count_data(classification_type, condition, **kwargs):
    ###################################################
    # connection to the DB
    global host
    global user
    global password
    global database

    enzyme_counts = collections.OrderedDict()
    conn = mysql.connector.connect(host=host, user=user, password=password, database=database)
    cursor = conn.cursor()
    if classification_type == "all":
        cursor.execute("SELECT EnzymeNode_id, strain, count "
                       "FROM RNAcounts "
                       "WHERE exp_condition='" + condition + "' "
                                                             "ORDER BY strain;")
    elif classification_type == "RMS":
        height = kwargs.get('height', 2)
        cursor.execute("SELECT EnzymeNode_id, strain, count FROM RNAcounts "
                       "INNER JOIN Enzyme_RMS_CPD USING(EnzymeNode_id) "
                       "INNER JOIN hasEnzymeBackUp USING(EnzymeNode_id) "
                       "WHERE exp_condition='"+condition+"' "
                                                         "AND Enzyme_RMS_CPD.height="+str(height)+" AND RMSid IN (SELECT RMSid FROM Enzyme_RMS_CPD GROUP BY RMSid ) AND hasEnzymeBackUp.height="+str(height)+";")
    elif classification_type == "EC":
        cursor.execute("SELECT EnzymeNode_id, strain, count "
                       "FROM RNAcounts "
                       "INNER JOIN EnzymeNode USING(EnzymeNode_id) "
                       "INNER JOIN YeastMetaBase.UniProt_ProteinAnnotation USING(ORF_id) "
                       "INNER JOIN YeastMetaBase.UP_ECnumber_CPD USING(UP_id) "
                       "WHERE exp_condition='"+condition+"' AND isECcomplete=1 ;")

    counts = cursor.fetchall()
    for c in counts:
        # 0: EnzymeNode_id, 1=strain , 2 =count

        if c[0] in enzyme_counts.keys():
            if c[2] > 0:
                # enzymeCounts[c[0]][c[1]]=numpy.log2(c[2])
                enzyme_counts[c[0]][c[1]] = c[2]
        else:
            if c[2] > 0:
                enzyme_counts[c[0]] = collections.OrderedDict()
                # enzymeCounts[c[0]][c[1]]=numpy.log2(c[2])
                enzyme_counts[c[0]][c[1]] = c[2]

    conn.close()
    return enzyme_counts


########################################################################

def retrieve_enzyme_group_cpd(classification_type, height):
    ct_enzyme_cpd = {}
    enzyme_ct_cpd = {}

    global host
    global user
    global password
    global database

    conn = mysql.connector.connect(host=host, user=user, password=password, database=database)
    cursor = conn.cursor()

    if classification_type == "RMS":
        cursor.execute("SELECT EnzymeNode_id, RMSid FROM RNAcounts INNER JOIN Enzyme_RMS_CPD USING(EnzymeNode_id) "
                       "WHERE height="+str(height)+" AND RMSid IN "
                                                   "(SELECT RMSid "
                                                   "FROM Enzyme_RMS_CPD GROUP BY RMSid ) "
                                                   "GROUP BY EnzymeNode_id, RMSid;")
        cpds = cursor.fetchall()
        for cpd in cpds:
            if cpd[1] in ct_enzyme_cpd.keys():
                ct_enzyme_cpd[cpd[1]].append(cpd[0])
            else:
                ct_enzyme_cpd[cpd[1]] = []
                ct_enzyme_cpd[cpd[1]].append(cpd[0])

            if cpd[0] in enzyme_ct_cpd.keys():
                enzyme_ct_cpd[cpd[0]].append(cpd[1])
            else:
                enzyme_ct_cpd[cpd[0]] = []
                enzyme_ct_cpd[cpd[0]].append(cpd[1])

    elif classification_type == "EC":
        cursor.execute(
            "SELECT EnzymeNode_id, EC_id FROM RNAcounts INNER JOIN EnzymeNode USING(EnzymeNode_id) "
            "INNER JOIN YeastMetaBase.UniProt_ProteinAnnotation USING(ORF_id) "
            "INNER JOIN YeastMetaBase.UP_ECnumber_CPD USING(UP_id)  GROUP BY EnzymeNode_id, EC_id;")
        cpds = cursor.fetchall()
        for cpd in cpds:
            if cpd[1] in ct_enzyme_cpd.keys():
                ct_enzyme_cpd[cpd[1]].append(cpd[0])
            else:
                ct_enzyme_cpd[cpd[1]] = []
                ct_enzyme_cpd[cpd[1]].append(cpd[0])

            if cpd[0] in enzyme_ct_cpd.keys():
                enzyme_ct_cpd[cpd[0]].append(cpd[1])
            else:
                enzyme_ct_cpd[cpd[0]] = []
                enzyme_ct_cpd[cpd[0]].append(cpd[1])

    elif classification_type == "all":
        # first RMS
        cursor.execute(
            "SELECT EnzymeNode_id, RMSid FROM RNAcounts INNER JOIN Enzyme_RMS_CPD USING(EnzymeNode_id) WHERE height=" + str(
                height) + " AND RMSid IN (SELECT RMSid FROM Enzyme_RMS_CPD GROUP BY RMSid ) GROUP BY EnzymeNode_id, RMSid;")
        cpds = cursor.fetchall()
        for cpd in cpds:
            if cpd[1] in ct_enzyme_cpd.keys():
                ct_enzyme_cpd[cpd[1]].append(cpd[0])
            else:
                ct_enzyme_cpd[cpd[1]] = []
                ct_enzyme_cpd[cpd[1]].append(cpd[0])

            if cpd[0] in enzyme_ct_cpd.keys():
                enzyme_ct_cpd[cpd[0]].append(cpd[1])
            else:
                enzyme_ct_cpd[cpd[0]] = []
                enzyme_ct_cpd[cpd[0]].append(cpd[1])

        # then EC numbers
        cursor.execute("SELECT EnzymeNode_id, EC_id FROM RNAcounts INNER JOIN EnzymeNode USING(EnzymeNode_id) "
                       "INNER JOIN YeastMetaBase.UniProt_ProteinAnnotation USING(ORF_id) "
                       "INNER JOIN YeastMetaBase.UP_ECnumber_CPD USING(UP_id) WHERE isECcomplete=1  GROUP BY EnzymeNode_id, EC_id;")
        cpds = cursor.fetchall()
        for cpd in cpds:
            if cpd[1] in ct_enzyme_cpd.keys():
                ct_enzyme_cpd[cpd[1]].append(cpd[0])
            else:
                ct_enzyme_cpd[cpd[1]] = []
                ct_enzyme_cpd[cpd[1]].append(cpd[0])

            if cpd[0] in enzyme_ct_cpd.keys():
                enzyme_ct_cpd[cpd[0]].append(cpd[1])
            else:
                enzyme_ct_cpd[cpd[0]] = []
                enzyme_ct_cpd[cpd[0]].append(cpd[1])

    conn.close()
    return ct_enzyme_cpd, enzyme_ct_cpd


########################################################################################################################

def get_genetic_interactions(gen_int_type):
    ###################################################
    # connection to the DB
    global host
    global user
    global password
    global database

    pairs = []

    conn = mysql.connector.connect(host=host, user=user, password=password, database=database)
    cursor = conn.cursor()

    cursor.execute("SELECT EnzymeNode_id_A, EnzymeNode_id_B FROM BIOGRID "
                   "WHERE experimentalSystem='"+gen_int_type+" Genetic' ;")
    counts = cursor.fetchall()
    for c in counts:
        # enzA enzB
        pairs.append([c[0], c[1]])

    conn.close()

    return pairs

#######################################################################################################################


def retrieve_enzymes_to_study(group_type, height, with_backup_only):
    enzyme_list = []
    global host
    global user
    global password
    global database

    conn = mysql.connector.connect(host=host, user=user, password=password, database=database)
    cursor = conn.cursor()

    if with_backup_only:
        hasBackup = "b.hasBackup=1 AND"
    else:
        hasBackup = ""

    if group_type=="all":
        cursor.execute(
            "SELECT DISTINCT(b.EnzymeNode_id) FROM hasEnzymeBackUp b "
            "INNER JOIN RNAcounts r1 USING(EnzymeNode_id) "
            "INNER JOIN RNAcounts r2 USING(EnzymeNode_id)  WHERE "+hasBackup+" (b.height=" + str(
                height) + " OR b.height=-1) and r1.exp_condition='N' AND r2.exp_condition='H';")

    elif group_type == "RMS":
        cursor.execute(
            "SELECT DISTINCT(b.EnzymeNode_id) FROM hasEnzymeBackUp b "
            "INNER JOIN RNAcounts r1 USING(EnzymeNode_id) "
            "INNER JOIN RNAcounts r2 USING(EnzymeNode_id)  WHERE "+hasBackup+" b.height=" + str(
                height) + " and r1.exp_condition='N' AND r2.exp_condition='H';")

    elif group_type == "EC":
        cursor.execute(
            "SELECT DISTINCT(b.EnzymeNode_id) FROM hasEnzymeBackUp b "
            "INNER JOIN RNAcounts r1 USING(EnzymeNode_id) "
            "INNER JOIN RNAcounts r2 USING(EnzymeNode_id)  "
            "WHERE "+hasBackup+" b.height=-1 and r1.exp_condition='N' AND r2.exp_condition='H';")

    enz = cursor.fetchall()
    for e in enz:
        enzyme_list.append(e[0])
    conn.close()

    return enzyme_list

######################################################################


def retrieve_protein_count_data(classification_type, condition, **kwargs):
    ###################################################
    # connection to the DB
    global host
    global user
    global password
    global database

    enzyme_counts = collections.OrderedDict()
    conn = mysql.connector.connect(host=host, user=user, password=password, database=database)
    cursor = conn.cursor()
    if classification_type == "all":
        cursor.execute("SELECT EnzymeNode_id, strain, count "
                       "FROM ProtCounts "
                       "WHERE exp_condition='" + condition + "' "
                                                             "ORDER BY strain;")
    elif classification_type == "RMS":
        height = kwargs.get('height', 2)
        cursor.execute("SELECT EnzymeNode_id, strain, count FROM ProtCounts "
                       "INNER JOIN Enzyme_RMS_CPD USING(EnzymeNode_id) "
                       "INNER JOIN hasEnzymeBackUp USING(EnzymeNode_id) "
                       "WHERE exp_condition='"+condition+"' "
                                                         "AND Enzyme_RMS_CPD.height="+str(height)+" AND RMSid IN (SELECT RMSid FROM Enzyme_RMS_CPD GROUP BY RMSid ) AND hasEnzymeBackUp.height="+str(height)+";")
    elif classification_type == "EC":
        cursor.execute("SELECT EnzymeNode_id, strain, count "
                       "FROM ProtCounts "
                       "INNER JOIN EnzymeNode USING(EnzymeNode_id) "
                       "INNER JOIN YeastMetaBase.UniProt_ProteinAnnotation USING(ORF_id) "
                       "INNER JOIN YeastMetaBase.UP_ECnumber_CPD USING(UP_id) "
                       "WHERE exp_condition='"+condition+"' AND isECcomplete=1 ;")

    counts = cursor.fetchall()
    for c in counts:
        # 0: EnzymeNode_id, 1=strain , 2 =count

        if c[0] in enzyme_counts.keys():
            if c[2] > 0:
                # enzymeCounts[c[0]][c[1]]=numpy.log2(c[2])
                enzyme_counts[c[0]][c[1]] = c[2]
        else:
            if c[2] > 0:
                enzyme_counts[c[0]] = collections.OrderedDict()
                # enzymeCounts[c[0]][c[1]]=numpy.log2(c[2])
                enzyme_counts[c[0]][c[1]] = c[2]

    conn.close()
    return enzyme_counts


########################################################################


def retrieve_enzyme_group_cpd_protein(classification_type, height):
    ct_enzyme_cpd = {}
    enzyme_ct_cpd = {}

    global host
    global user
    global password
    global database

    conn = mysql.connector.connect(host=host, user=user, password=password, database=database)
    cursor = conn.cursor()

    if classification_type == "RMS":
        cursor.execute("SELECT EnzymeNode_id, RMSid FROM ProtCounts INNER JOIN Enzyme_RMS_CPD USING(EnzymeNode_id) "
                       "WHERE height="+str(height)+" AND RMSid IN "
                                                   "(SELECT RMSid "
                                                   "FROM Enzyme_RMS_CPD GROUP BY RMSid ) "
                                                   "GROUP BY EnzymeNode_id, RMSid;")
        cpds = cursor.fetchall()
        for cpd in cpds:
            if cpd[1] in ct_enzyme_cpd.keys():
                ct_enzyme_cpd[cpd[1]].append(cpd[0])
            else:
                ct_enzyme_cpd[cpd[1]] = []
                ct_enzyme_cpd[cpd[1]].append(cpd[0])

            if cpd[0] in enzyme_ct_cpd.keys():
                enzyme_ct_cpd[cpd[0]].append(cpd[1])
            else:
                enzyme_ct_cpd[cpd[0]] = []
                enzyme_ct_cpd[cpd[0]].append(cpd[1])

    elif classification_type == "EC":
        cursor.execute(
            "SELECT EnzymeNode_id, EC_id FROM ProtCounts INNER JOIN EnzymeNode USING(EnzymeNode_id) "
            "INNER JOIN YeastMetaBase.UniProt_ProteinAnnotation USING(ORF_id) "
            "INNER JOIN YeastMetaBase.UP_ECnumber_CPD USING(UP_id)  GROUP BY EnzymeNode_id, EC_id;")
        cpds = cursor.fetchall()
        for cpd in cpds:
            if cpd[1] in ct_enzyme_cpd.keys():
                ct_enzyme_cpd[cpd[1]].append(cpd[0])
            else:
                ct_enzyme_cpd[cpd[1]] = []
                ct_enzyme_cpd[cpd[1]].append(cpd[0])

            if cpd[0] in enzyme_ct_cpd.keys():
                enzyme_ct_cpd[cpd[0]].append(cpd[1])
            else:
                enzyme_ct_cpd[cpd[0]] = []
                enzyme_ct_cpd[cpd[0]].append(cpd[1])

    elif classification_type == "all":
        # first RMS
        cursor.execute(
            "SELECT EnzymeNode_id, RMSid FROM ProtCounts INNER JOIN Enzyme_RMS_CPD USING(EnzymeNode_id) WHERE height=" + str(
                height) + " AND RMSid IN (SELECT RMSid FROM Enzyme_RMS_CPD GROUP BY RMSid ) GROUP BY EnzymeNode_id, RMSid;")
        cpds = cursor.fetchall()
        for cpd in cpds:
            if cpd[1] in ct_enzyme_cpd.keys():
                ct_enzyme_cpd[cpd[1]].append(cpd[0])
            else:
                ct_enzyme_cpd[cpd[1]] = []
                ct_enzyme_cpd[cpd[1]].append(cpd[0])

            if cpd[0] in enzyme_ct_cpd.keys():
                enzyme_ct_cpd[cpd[0]].append(cpd[1])
            else:
                enzyme_ct_cpd[cpd[0]] = []
                enzyme_ct_cpd[cpd[0]].append(cpd[1])

        # then EC numbers
        cursor.execute("SELECT EnzymeNode_id, EC_id FROM ProtCounts INNER JOIN EnzymeNode USING(EnzymeNode_id) "
                       "INNER JOIN YeastMetaBase.UniProt_ProteinAnnotation USING(ORF_id) "
                       "INNER JOIN YeastMetaBase.UP_ECnumber_CPD USING(UP_id) WHERE isECcomplete=1  GROUP BY EnzymeNode_id, EC_id;")
        cpds = cursor.fetchall()
        for cpd in cpds:
            if cpd[1] in ct_enzyme_cpd.keys():
                ct_enzyme_cpd[cpd[1]].append(cpd[0])
            else:
                ct_enzyme_cpd[cpd[1]] = []
                ct_enzyme_cpd[cpd[1]].append(cpd[0])

            if cpd[0] in enzyme_ct_cpd.keys():
                enzyme_ct_cpd[cpd[0]].append(cpd[1])
            else:
                enzyme_ct_cpd[cpd[0]] = []
                enzyme_ct_cpd[cpd[0]].append(cpd[1])

    conn.close()
    return ct_enzyme_cpd, enzyme_ct_cpd


########################################################################################################################


def retrieve_enzymes_to_study_protein(group_type, height):
    enzyme_list = []
    global host
    global user
    global password
    global database

    conn = mysql.connector.connect(host=host, user=user, password=password, database=database)
    cursor = conn.cursor()

    if group_type == "all":
        cursor.execute(
            "SELECT DISTINCT(b.EnzymeNode_id) FROM hasEnzymeBackUp b "
            "INNER JOIN ProtCounts r1 USING(EnzymeNode_id) "
            "INNER JOIN ProtCounts r2 USING(EnzymeNode_id)  WHERE b.hasBackup=1 AND (b.height=" + str(
                height) + " OR b.height=-1) and r1.exp_condition='N' AND r2.exp_condition='H';")

    elif group_type == "RMS":
        cursor.execute(
            "SELECT DISTINCT(b.EnzymeNode_id) FROM hasEnzymeBackUp b "
            "INNER JOIN ProtCounts r1 USING(EnzymeNode_id) "
            "INNER JOIN ProtCounts r2 USING(EnzymeNode_id)  WHERE b.hasBackup=1 AND b.height=" + str(
                height) + " and r1.exp_condition='N' AND r2.exp_condition='H';")

    elif group_type == "EC":
        cursor.execute(
            "SELECT DISTINCT(b.EnzymeNode_id) FROM hasEnzymeBackUp b "
            "INNER JOIN ProtCounts r1 USING(EnzymeNode_id) "
            "INNER JOIN ProtCounts r2 USING(EnzymeNode_id)  "
            "WHERE b.hasBackup=1 AND b.height=-1 and r1.exp_condition='N' AND r2.exp_condition='H';")

    enz = cursor.fetchall()
    for e in enz:
        enzyme_list.append(e[0])
    conn.close()

    return enzyme_list

######################################################################

def retrieve_count_data_cer(classification_type, **kwargs):
    ###################################################
    # connection to the DB
    global host
    global user
    global password
    global database

    enzyme_counts = collections.OrderedDict()
    conn = mysql.connector.connect(host=host, user=user, password=password, database=database)
    cursor = conn.cursor()
    if classification_type == "all":
        cursor.execute("SELECT EnzymeNode_id, strain, count "
                       "FROM RNAcounts ORDER BY strain;")
    elif classification_type == "RMS":
        height = kwargs.get('height', 2)
        cursor.execute("SELECT EnzymeNode_id, strain, count FROM RNAcounts "
                       "INNER JOIN Enzyme_RMS_CPD USING(EnzymeNode_id) "
                       "INNER JOIN hasEnzymeBackUp USING(EnzymeNode_id) "
                       "WHERE Enzyme_RMS_CPD.height="+str(height)+" AND RMSid IN (SELECT RMSid FROM Enzyme_RMS_CPD GROUP BY RMSid ) AND hasEnzymeBackUp.height="+str(height)+";")
    elif classification_type == "EC":
        cursor.execute("SELECT EnzymeNode_id, strain, count "
                       "FROM RNAcounts "
                       "INNER JOIN EnzymeNode USING(EnzymeNode_id) "
                       "INNER JOIN YeastMetaBase.UniProt_ProteinAnnotation USING(ORF_id) "
                       "INNER JOIN YeastMetaBase.UP_ECnumber_CPD USING(UP_id) "
                       "WHERE isECcomplete=1 ;")

    counts = cursor.fetchall()
    for c in counts:
        # 0: EnzymeNode_id, 1=strain , 2 =count

        if c[0] in enzyme_counts.keys():
            if c[2] > 0:
                # enzymeCounts[c[0]][c[1]]=numpy.log2(c[2])
                enzyme_counts[c[0]][c[1]] = c[2]
        else:
            if c[2] > 0:
                enzyme_counts[c[0]] = collections.OrderedDict()
                # enzymeCounts[c[0]][c[1]]=numpy.log2(c[2])
                enzyme_counts[c[0]][c[1]] = c[2]

    conn.close()
    return enzyme_counts

########################################################################


def retrieve_count_data_prot_cer(classification_type, **kwargs):
    ###################################################
    # connection to the DB
    global host
    global user
    global password
    global database

    enzyme_counts = collections.OrderedDict()
    conn = mysql.connector.connect(host=host, user=user, password=password, database=database)
    cursor = conn.cursor()
    if classification_type == "all":
        cursor.execute("SELECT EnzymeNode_id, strain, count "
                       "FROM ProtCounts ORDER BY strain;")

    elif classification_type == "RMS":
        height = kwargs.get('height', 2)
        cursor.execute("SELECT EnzymeNode_id, strain, count FROM ProtCounts "
                       "INNER JOIN Enzyme_RMS_CPD USING(EnzymeNode_id) "
                       "INNER JOIN hasEnzymeBackUp USING(EnzymeNode_id) "
                       "WHERE Enzyme_RMS_CPD.height="+str(height)+" AND RMSid IN (SELECT RMSid FROM Enzyme_RMS_CPD GROUP BY RMSid ) AND hasEnzymeBackUp.height="+str(height)+";")
    elif classification_type == "EC":
        cursor.execute("SELECT EnzymeNode_id, strain, count "
                       "FROM ProtCounts "
                       "INNER JOIN EnzymeNode USING(EnzymeNode_id) "
                       "INNER JOIN YeastMetaBase.UniProt_ProteinAnnotation USING(ORF_id) "
                       "INNER JOIN YeastMetaBase.UP_ECnumber_CPD USING(UP_id) "
                       "WHERE isECcomplete=1 ;")

    counts = cursor.fetchall()
    for c in counts:
        # 0: EnzymeNode_id, 1=strain , 2 =count

        if c[0] in enzyme_counts.keys():
            if c[2] > 0:
                # enzymeCounts[c[0]][c[1]]=numpy.log2(c[2])
                enzyme_counts[c[0]][c[1]] = c[2]
        else:
            if c[2] > 0:
                enzyme_counts[c[0]] = collections.OrderedDict()
                # enzymeCounts[c[0]][c[1]]=numpy.log2(c[2])
                enzyme_counts[c[0]][c[1]] = c[2]

    conn.close()
    return enzyme_counts

########################################################################


def retrieve_enzymes_to_study_cer(group_type, height):
    enzyme_list = []
    global host
    global user
    global password
    global database

    conn = mysql.connector.connect(host=host, user=user, password=password, database=database)
    cursor = conn.cursor()

    if group_type == "all":
        cursor.execute(
            "SELECT DISTINCT(b.EnzymeNode_id) FROM hasEnzymeBackUp b "
            "INNER JOIN RNAcounts r1 USING(EnzymeNode_id) "
            "WHERE b.hasBackup=1 AND (b.height=" + str(height) + " OR b.height=-1) ;")

    elif group_type == "RMS":
        cursor.execute(
            "SELECT DISTINCT(b.EnzymeNode_id) FROM hasEnzymeBackUp b "
            "INNER JOIN RNAcounts r1 USING(EnzymeNode_id) "
            "WHERE b.hasBackup=1 AND b.height=" + str(height) + " ;")

    elif group_type == "EC":
        cursor.execute(
            "SELECT DISTINCT(b.EnzymeNode_id) FROM hasEnzymeBackUp b "
            "INNER JOIN RNAcounts r1 USING(EnzymeNode_id) WHERE b.hasBackup=1 AND b.height=-1 ;")

    enz = cursor.fetchall()
    for e in enz:
        enzyme_list.append(e[0])
    conn.close()

    return enzyme_list

######################################################################


def retrieve_enzymes_to_study_cer_prot(group_type, height):
    enzyme_list = []
    global host
    global user
    global password
    global database

    conn = mysql.connector.connect(host=host, user=user, password=password, database=database)
    cursor = conn.cursor()

    if group_type == "all":
        cursor.execute(
            "SELECT DISTINCT(b.EnzymeNode_id) FROM hasEnzymeBackUp b "
            "INNER JOIN ProtCounts r1 USING(EnzymeNode_id) "
            "WHERE b.hasBackup=1 AND (b.height=" + str(height) + " OR b.height=-1) ;")

    elif group_type == "RMS":
        cursor.execute(
            "SELECT DISTINCT(b.EnzymeNode_id) FROM hasEnzymeBackUp b "
            "INNER JOIN ProtCounts r1 USING(EnzymeNode_id) "
            "WHERE b.hasBackup=1 AND b.height=" + str(height) + " ;")

    elif group_type == "EC":
        cursor.execute(
            "SELECT DISTINCT(b.EnzymeNode_id) FROM hasEnzymeBackUp b "
            "INNER JOIN ProtCounts r1 USING(EnzymeNode_id) WHERE b.hasBackup=1 AND b.height=-1 ;")

    enz = cursor.fetchall()
    for e in enz:
        enzyme_list.append(e[0])
    conn.close()

    return enzyme_list

######################################################################


# this function allows to retrieve pairs of enzymes (and their general statistics) that are in same QTL (exp_condition!)
# and according to precised conditions (on their correlation or sdReduction) if precised
def retrieve_pairs_with_qtl(**kwargs):

    connector = kwargs.get('connector', None)
    data_type = kwargs.get('data_type', 'RNA')

    condition = kwargs.get('exp_condition', 'both')

    if connector is None:
        print("Can't connect to the database!")
        return


    #retrieving scores if precised

    correlationN = kwargs.get("correlationN", None)  # shoud be a string like ">=0.5"
    correlationH = kwargs.get("correlationH", None)
    sdReductionN = kwargs.get("sdReductionN", None)  # should be a string like "<=0"
    sdReductionH = kwargs.get("sdReductionH", None)

    # construct query

    condition1 = ""
    condition2 = ""
    condition3 = ""
    condition4 = ""

    if correlationN is not None:
        condition1 = " AND pt.correlationN"+str(correlationN)+" "
    if correlationH is not None:
        condition2 = " AND pt.correlationH"+str(correlationH)+" "
    if sdReductionN is not None:
        condition3 = " AND pt.sdReductionN"+str(sdReductionN)+" "
    if sdReductionH is not None:
        condition4 = " AND pt.sdReductionH"+str(sdReductionH)+" "

    # change the * : EnzymeNode_id_A, EnzymeNode_id_B, GROUP_CONCAT(DISTINCT(qtl1.marker), SEPARATOR '$') .... GROUP BY EnzymeNode_id_A, EnzymeNode_id_B
    query_qtl_inN = "SELECT pt.EnzymeNode_id_A, pt.EnzymeNode_id_B, GROUP_CONCAT(DISTINCT(qtl1.marker) SEPARATOR '$') " \
                    "FROM PairData pt " \
                    "INNER JOIN EnzymeNode en1 ON(pt.EnzymeNode_id_A=en1.EnzymeNode_id) " \
                    "INNER JOIN YeastMetaBase.eQTLmappings qtl1 ON(en1.ORF_id=qtl1.ORF_id) " \
                    "INNER JOIN EnzymeNode en2 ON(pt.EnzymeNode_id_B=en2.EnzymeNode_id) " \
                    "INNER JOIN YeastMetaBase.eQTLmappings qtl2 ON(en2.ORF_id=qtl2.ORF_id) " \
                    "INNER JOIN Enzyme_RMS_CPD trms ON (pt.EnzymeNode_id_A = trms.EnzymeNode_id)" \
                    "WHERE pt.dataType ='" + data_type + "' " \
                    "AND pt.inSameGroup='True' " + condition1 + condition2 + condition3 + condition4 + " " \
                    "AND qtl1.fdr<=0.05 AND qtl2.fdr<=0.05 AND qtl1.exp_condition='N' AND qtl2.exp_condition='N' " \
                    "AND qtl1.marker=qtl2.marker " \
                    "AND trms.RMSid NOT LIKE 'RMS-0.1833%' AND trms.RMSid NOT LIKE 'RMS-0.2613%' " \
                    "GROUP BY pt.EnzymeNode_id_A, pt.EnzymeNode_id_B;"

                    ##"AND trms.RMSid NOT LIKE 'RMS-0.1048%'

    query_qtl_inH = "SELECT pt.EnzymeNode_id_A, pt.EnzymeNode_id_B, GROUP_CONCAT(DISTINCT(qtl1.marker) SEPARATOR '$') " \
                    "FROM PairData pt " \
                    "INNER JOIN EnzymeNode en1 ON(pt.EnzymeNode_id_A=en1.EnzymeNode_id) " \
                    "INNER JOIN YeastMetaBase.eQTLmappings qtl1 ON(en1.ORF_id=qtl1.ORF_id) " \
                    "INNER JOIN EnzymeNode en2 ON(pt.EnzymeNode_id_B=en2.EnzymeNode_id) " \
                    "INNER JOIN YeastMetaBase.eQTLmappings qtl2 ON(en2.ORF_id=qtl2.ORF_id) " \
                    "INNER JOIN Enzyme_RMS_CPD trms ON (pt.EnzymeNode_id_A = trms.EnzymeNode_id)" \
                    "WHERE pt.dataType ='" + data_type + "' " \
                    "AND pt.inSameGroup='True' " + condition1 + condition2 + condition3 + condition4 + " " \
                    "AND qtl1.fdr<=0.05 AND qtl2.fdr<=0.05 AND qtl1.exp_condition='H' AND qtl2.exp_condition='H' " \
                    "AND qtl1.marker=qtl2.marker " \
                    "AND trms.RMSid NOT LIKE 'RMS-0.1833%' AND trms.RMSid NOT LIKE 'RMS-0.2613%' " \
                    "GROUP BY pt.EnzymeNode_id_A, pt.EnzymeNode_id_B;"

                    #"AND trms.RMSid NOT LIKE 'RMS-0.1048%'

    cursor = connector.connect_to_db()

    biglist = {}

    if condition == "N":
        # run just query_qtl_inN
        print(query_qtl_inN)
        cursor.execute(query_qtl_inN)
        pairs = cursor.fetchall()
        for line in pairs:
            biglist[(line[0], line[1])] = line[2].split("$")
            if line[0]=="idn1-SPAC4G9.12" and line[1]=="SPCC162.11c":
                print("FOUND!")
    elif condition == "H":
        # run just query_qtl_inH
        print(query_qtl_inH)
        cursor.execute(query_qtl_inH)
        pairs = cursor.fetchall()
        for line in pairs:
            biglist[(line[0], line[1])] = line[2].split("$")
    elif condition == "both":

        biglistN = {}
        print(query_qtl_inN)
        cursor.execute(query_qtl_inN)
        pairs = cursor.fetchall()
        for line in pairs:
            biglistN[(line[0], line[1])] = line[2].split("$")

        biglistH = {}
        print(query_qtl_inH)
        cursor.execute(query_qtl_inH)
        pairs = cursor.fetchall()
        for line in pairs:
            biglistH[(line[0], line[1])] = line[2].split("$")
        # + do filtering - select pairs of enzymes that are in both lists

        intersect = []
        for pair in biglistN.keys():
            if pair in biglistH.keys():
                intersect.append(pair)
        for pair in intersect:
            biglist[pair] = biglistN[pair]



    # in biglist, there are now pairs of enzymes with a list of eQTLs that they share
    # need to sort the markers, to keep only those that don't segregate strains in same groups
    simplified_biglist = sort_qtl_markers(biglist, cursor)

    return simplified_biglist



def sort_qtl_markers(biglist, cursor):

    simplified_biglist = {}

    for pair in biglist.keys():
        markers = biglist[pair]
        simplified_biglist[pair]=[]

        if len(markers) == 1:

            #nothing to do
            simplified_biglist[pair] = markers
        else:
            # need to sort
            strains_by_marker = {}  # key = marker, value = [ set0, set1 ] sets are made of strains
            for marker in markers:
                strains_by_marker[marker] = [[], []]
                query = "SELECT marker, strain, allele from YeastMetaBase.GenotypeReference WHERE marker="+marker+";"
                cursor.execute(query)
                res = cursor.fetchall()
                flag = 0
                for line in res:
                    strain = line[1]
                    allele = line[2]
                    #print(allele)
                    if allele == "1":

                        strains_by_marker[marker][1].append(strain)
                    elif allele == "0":
                        strains_by_marker[marker][0].append(strain)
                    elif allele == "NA":
                        flag=1

                for omarker in strains_by_marker.keys():
                    if omarker != marker:
                        # compare their strain sets
                        if set(strains_by_marker[marker][0]) == set(strains_by_marker[omarker][0]) or set(strains_by_marker[marker][0]) == set(strains_by_marker[omarker][1]):
                            flag=1

                if flag==1:
                    strains_by_marker.pop(marker)

            for m in strains_by_marker.keys():
                simplified_biglist[pair].append(m)

    return simplified_biglist


########################################################################################################################


def network_pair_qtl_analysis(**kwargs):
    connector = kwargs.get('connector', None)
    data_type = kwargs.get('data_type', 'RNA')
    condition = kwargs.get('exp_condition', 'both')
    pair = kwargs.get('pair', None)

    if connector is None:
        print("Can't connect to the database!")
        return

    if pair is None:
        print("No pair input, return")
        return


    return




