#!/usr/bin/python
from __future__ import print_function

from RMSbufferingMethods import dataImportMethods
from RMSbufferingMethods import dataOutputMethods
from RMSbufferingMethods import MySQLConnector

import sys
import re
import os
import time
import random
import numpy
import scipy
import mysql.connector
import itertools
from itertools import chain, combinations
import collections
from collections import defaultdict
from scipy.stats.stats import pearsonr
from scipy.stats.mstats import *
from scipy.stats import ttest_ind
import operator
import matplotlib.pyplot as plt
import matplotlib.patches as mpatch

########################################################################################################################

host = ""
user = ""
password = ""
database = ""
directory = ""
organism = ""
dataType = ""
enzymeClassification = ""
height = ""


def analyse_pair(pair, enzyme_counts):
    sump = enzyme_counts[pair[0]].copy()
    for strain in enzyme_counts[pair[1]].keys():
            if strain in sump.keys():
                sump[strain]=sump[strain]+enzyme_counts[pair[1]][strain]

    pair_sd = numpy.std([numpy.log2(x) for x in list(sump.values()) if x > 0], dtype=numpy.float64)

    pair_gmean = gmean([numpy.log2(x) for x in list(sump.values()) if x > 0], dtype=numpy.float64)

    pair_correlation=0
    if len(list(enzyme_counts[pair[0]].values())) == len(list(enzyme_counts[pair[1]].values())):
        try:
            pair_correlation = scipy.stats.pearsonr(list(enzyme_counts[pair[0]].values()), list(enzyme_counts[pair[1]].values()))
        except FloatingPointError as err:
            print(str(pair[0])+"\t"+str(pair[1])+"\t"+str(err))

    return pair_gmean, pair_sd, pair_correlation

########################################################################################################################


def analyse_pair_fold_change(pair, enz_n, enz_h):
    numpy.seterr(all='ignore', divide='ignore', over='ignore', under='ignore', invalid='ignore')

    fcpair = {}  # fold change of the sum for each strain
    fcenz1 = {}  # fold change for the enzyme 1
    fcenz2 = {}  # fold change for the enzyme 2
    for strain in enz_h[pair[0]].keys():
        if strain in enz_n[pair[0]] and strain in enz_n[pair[1]] and strain in enz_h[pair[1]] and enz_n[pair[0]][strain] != 0 and enz_n[pair[1]][strain] != 0:
                fcenz1[strain] = numpy.log2(enz_h[pair[0]][strain]) / numpy.log2(enz_n[pair[0]][strain])
                fcenz2[strain] = numpy.log2(enz_h[pair[1]][strain]) / numpy.log2(enz_n[pair[1]][strain])

                fcpair[strain] = (numpy.log2(enz_h[pair[0]][strain] + enz_h[pair[1]][strain])) / (numpy.log2(enz_n[pair[0]][strain] + enz_n[pair[1]][strain]))

    mean_fc_pair = numpy.mean(list(fcpair.values()), dtype=numpy.float64)  # mean of the fold change of the pair
    sd_fc_pair = numpy.std(list(fcpair.values()), dtype=numpy.float64)  # SD of the fold change
    # varFCpair = numpy.var(list(fcpair.values()))

    # meanFCenz1 = mean( list(fcenz1.values()) )
    sd_fc_enz1 = numpy.std(list(fcenz1.values()))
    # varFCenz1 = numpy.var(list(fcenz1.values()))

    # meanFCenz2 = mean( list(fcenz2.values()) )
    sd_fc_enz2 = numpy.std(list(fcenz2.values()))
    # varFCenz2 = numpy.var(list(fcenz2.values()))

    sd_reduction_fc_min = sd_fc_pair - min([sd_fc_enz1, sd_fc_enz2])
    # varReductionFC_min = varFCpair - min([varFCenz1, varFCenz2])

    return sd_fc_pair, sd_reduction_fc_min, mean_fc_pair


#######################################################################################################################


def is_in_same_group(pair, enz_g_cpd):
    anw = False
    if len(pair) == 2:
        if not set(enz_g_cpd[pair[0]] ).isdisjoint(enz_g_cpd[pair[1]]):
            anw = True
    return anw

#######################################################################################################################


def analyse_in_rms_not_in_ec(**kwargs):
    height = kwargs.get('height', 2)
    connector = kwargs.get('connector', None)
    data_type = kwargs.get('data_type', 'RNA')

    if connector is None:
        print("Can't connect to the database!")
        return

    enzymes_without_ec = []

    cursor = connector.connect_to_db()

    #fetching enzymes that are not in EC
    query = "SELECT DISTINCT(EnzymeNode_id) " \
            "FROM  Enzyme_RMS_CPD WHERE EnzymeNode_id NOT IN(" \
            "select EnzymeNode_id from EnzymeNode INNER JOIN YeastMetaBase.UniProt_ProteinAnnotation USING(ORF_id) " \
            "INNER JOIN YeastMetaBase.UP_ECnumber_CPD USING(UP_id) GROUP BY EnzymeNode_id, EC_id);"
    cursor.execute(query)
    enzymes = cursor.fetchall()
    for enz in enzymes:
        enzymes_without_ec.append(enz[0])

    #fetching pair data

    query = "SELECT * FROM PairData WHERE dataType='"+data_type+"';"

    # +-----------------+-------------+------+-----+---------+-------+
    # | Field | Type | Null | Key | Default | Extra |
    # +-----------------+-------------+------+-----+---------+-------+
    # | EnzymeNode_id_A | varchar(20) | YES | MUL | NULL | |
    # | EnzymeNode_id_B | varchar(20) | YES | MUL | NULL | |
    # |2 gmeanSumN | float | YES | | NULL | |
    # |3 sdSumN | float | YES | | NULL | |
    # |4 correlationN | float | YES | | NULL | |
    # |5 gmeanSumH | float | YES | | NULL | |
    # |6 sdSumH | float | YES | | NULL | |
    # |7 correlationH | float | YES | | NULL | |
    # | deltaSD | float | YES | | NULL | |
    # | sdReductionN | float | YES | | NULL | |
    # | sdReductionH | float | YES | | NULL | |
    # | sdFC | float | YES | | NULL | |
    # | sdReductionFC | float | YES | | NULL | |
    # | meanFC | float | YES | | NULL | |
    # | inSameGroup | varchar(6) | YES | | NULL | |
    # | dataType | varchar(8) | YES | | NULL | |
    # +-----------------+-------------+------+-----+---------+-------+

    cursor.execute(query)
    pairs = cursor.fetchall()

    pairs_without_ec_in_same_rms = {}
    pairs_without_ec_not_in_same_rms = {}
    pairs_with_ec_in_same_rms={}
    pairs_with_ec_not_in_same_rms = {}

    for line in pairs:
        if line[0] in enzymes_without_ec or line[1] in enzymes_without_ec:
            # at least one of the two enzymes doesn't have an EC
            if line[14] == 'True':  # pair in same RMS
                pairs_without_ec_in_same_rms[(line[0], line[1])] = [line[2], line[3], line[4], line[5], line[6], line[7], line[8], line[9], line[10], line[11], line[12], line[13]]
            elif line[14] == 'False':  # pair not in same RMS
                pairs_without_ec_not_in_same_rms[(line[0], line[1])] = [line[2], line[3], line[4], line[5], line[6], line[7], line[8], line[9], line[10], line[11], line[12], line[13]]
        else:
            # both enzymes in the pair have an EC
            if line[14] == 'True':  # pair in same RMS
                pairs_with_ec_in_same_rms[(line[0], line[1])] = [line[2], line[3], line[4], line[5], line[6], line[7], line[8], line[9], line[10], line[11], line[12], line[13]]
            elif line[14] == 'False':
                pairs_with_ec_not_in_same_rms[(line[0], line[1])] = [line[2], line[3], line[4], line[5], line[6], line[7], line[8], line[9], line[10], line[11], line[12], line[13]]

    # stat analyses

    #############################
    # plots
    #############################

    #############################

    prepare_data_for_plots_without_ec(0, pairs_without_ec_in_same_rms, pairs_with_ec_in_same_rms)
    prepare_data_for_plots_without_ec(1, pairs_without_ec_in_same_rms, pairs_with_ec_in_same_rms)
    prepare_data_for_plots_without_ec(2, pairs_without_ec_in_same_rms, pairs_with_ec_in_same_rms)
    prepare_data_for_plots_without_ec(3, pairs_without_ec_in_same_rms, pairs_with_ec_in_same_rms)
    prepare_data_for_plots_without_ec(4, pairs_without_ec_in_same_rms, pairs_with_ec_in_same_rms)
    prepare_data_for_plots_without_ec(5, pairs_without_ec_in_same_rms, pairs_with_ec_in_same_rms)
    prepare_data_for_plots_without_ec(9, pairs_without_ec_in_same_rms, pairs_with_ec_in_same_rms)
    prepare_data_for_plots_without_ec(10, pairs_without_ec_in_same_rms, pairs_with_ec_in_same_rms)
    prepare_data_for_plots_without_ec(11, pairs_without_ec_in_same_rms, pairs_with_ec_in_same_rms)




    # Correlation(enz1,enz2)

    connector.disconnect_from_db()

    return


#######################################################################################################################

def prepare_data_for_plots_without_ec(measure, list_without_ec, list_with_ec):

    if measure==0:
        title = "gmean N | pairs in same RMS | 1or2"
        xlab = "gmean(log2(enz1+enz2))"
    elif measure == 1:
        title = "SD N | pairs in same RMS | 1or2"
        xlab = "SD(log2(enz1+enz2))"
    elif measure == 2:
        title = "Correlation N | pairs in same RMS | 1or2"
        xlab = "Correlation(enz1,enz2)"
    elif measure == 3:
        title = "gmean H | pairs in same RMS | 1or2"
        xlab = "gmean(log2(enz1+enz2))"
    elif measure == 4:
        title = "SD H | pairs in same RMS | 1or2"
        xlab = "SD(log2(enz1+enz2))"
    elif measure == 5:
        title = "Correlation H | pairs in same RMS | 1or2"
        xlab = "Correlation(enz1,enz2)"
    elif measure == 9:
        title = "SD Fold Change | pairs in same RMS | 1or2"
        xlab = "FC(log2(enz1+enz2))"
    elif measure == 10:
        title = "SD reduction Fold Change | pairs in same RMS | 1or2"
        xlab = "FC(log2(enz1+enz2)) - min[FC(log2(enz1)), FC(log2(enz2))]"
    elif measure == 11:
        title = "mean Fold Change | pairs in same RMS | 1or2"
        xlab = "mean(log2(enz1+enz2))"

    list1 = [list_without_ec[pair][measure] for pair in list_without_ec.keys()]
    list2 = [list_with_ec[pair][measure] for pair in list_with_ec.keys()]
    label1 = "Pairs without EC number: " + str(len(list1))
    label2 = "Pairs with EC number: " + str(len(list2))
    dataOutputMethods.plot_hists_ec_analysis(list1, list2, label1, label2, title, xlab)

    return list1, list2, label1, label2, title, xlab


######################################################################################################################

def single_qtl_analysis(**kwargs):
    connector = kwargs.get('connector', None)
    data_type = kwargs.get('data_type', 'RNA')
    exp_condition = kwargs.get('exp_condition', 'N')
    pair = kwargs.get('pair', None)
    marker = kwargs.get("marker", None)

    enz1allele0 = {}
    enz1allele1 = {}
    enz2allele0 = {}
    enz2allele1 = {}

    # for the pair of enzymes, need to retrieve data in total 4 lists, by combination enzyme/allele
    # and them, compute scores - 2 sets of scores, 1 by allele, instead of one
    # + don't forget to do sums!!!!!

    if connector is None or pair is None or marker is None:
        print("Can't don anything without database or data")
        return
    cursor = connector.connect_to_db()
    if data_type == "RNA":
        queryEnz1 = "SELECT EnzymeNode_id, strain, count, allele FROM RNAcounts INNER JOIN YeastMetaBase.GenotypeReference USING(strain) WHERE EnzymeNode_id = '"+pair[0]+"' AND marker="+marker+" AND exp_condition='"+exp_condition+"';"
        queryEnz2 = "SELECT EnzymeNode_id, strain, count, allele FROM RNAcounts INNER JOIN YeastMetaBase.GenotypeReference USING(strain) WHERE EnzymeNode_id = '"+pair[1]+"' AND marker="+marker+" AND exp_condition='"+exp_condition+"';"

        cursor.execute(queryEnz1)

        enz1data = cursor.fetchall()
        for line in enz1data:
            strain = line[1]
            count = line[2]
            allele = line[3]
            if allele == "0":
                #add count to :
                enz1allele0[strain] = count
            else :
                #add count to :
                enz1allele1[strain] = count

        # for the second enzyme
        cursor.execute(queryEnz2)
        enz2data = cursor.fetchall()
        for line in enz2data:
            strain = line[1]
            count = line[2]
            allele = line[3]
            if allele == "0":
                # add count to :
                enz2allele0[strain] = count
            else :
                # add count to :
                enz2allele1[strain] = count

    # there are four tables, for each enzyme and each strain, there are counts

    # compare enz1allele0 to enz2alle0 AND enz1allele1 to enz2allele1
    allele0_pair_gmean, allele0_pair_sd, allele0_pair_correlation = analyse_splitted_pair(enz1allele0, enz2allele0)
    allele1_pair_gmean, allele1_pair_sd, allele1_pair_correlation = analyse_splitted_pair(enz1allele1, enz2allele1)

    allele0_sdReduction = allele0_pair_sd - min( numpy.std([numpy.log2(x) for x in list(enz1allele0.values()) if x > 0], dtype=numpy.float64), numpy.std([numpy.log2(x) for x in list(enz2allele0.values()) if x > 0], dtype=numpy.float64) )
    allele1_sdReduction = allele1_pair_sd - min( numpy.std([numpy.log2(x) for x in list(enz1allele1.values()) if x > 0], dtype=numpy.float64), numpy.std([numpy.log2(x) for x in list(enz2allele1.values()) if x > 0], dtype=numpy.float64) )



    # if isinstance(allele0_pair_correlation, tuple) and isinstance(allele1_pair_correlation, tuple) and numpy.absolute(allele0_pair_correlation[0] - allele1_pair_correlation[0]) > 0.1:
    #     print("real correlation difference in H condition")
    # else:
    #     print("no")

    if isinstance(allele0_pair_correlation, tuple) and isinstance(allele1_pair_correlation, tuple): #and numpy.sign(allele0_pair_correlation[0]) != numpy.sign(allele1_pair_correlation[0]) :

        #if numpy.absolute(allele0_pair_correlation[0] - allele1_pair_correlation[0]) >= 0.5:
        print("\n"+str(pair) + " " + marker)
        print("correlation allele0 " + str(allele0_pair_correlation))
        print("correlation allele1 " + str(allele1_pair_correlation))
        print("allele0 sd reduction " + str(allele0_sdReduction))
        print("allele1 sd reduction " + str(allele1_sdReduction))
        print("SD sum allele 0 "+str(allele0_pair_sd))
        print("SD sum allele 1 "+str(allele1_pair_sd))
        print("SD enzyme 1 allele 0 "+str(  numpy.std([numpy.log2(x) for x in list(enz1allele0.values()) if x > 0], dtype=numpy.float64) ))
        print("SD enzyme 1 allele 1 "+str( numpy.std([numpy.log2(x) for x in list(enz1allele1.values()) if x > 0], dtype=numpy.float64)  ))
        print("SD enzyme 2 allele 0 "+str( numpy.std([numpy.log2(x) for x in list(enz2allele0.values()) if x > 0], dtype=numpy.float64)  ))
        print("SD enzyme 2 allele 1 "+str(  numpy.std([numpy.log2(x) for x in list(enz2allele1.values()) if x > 0], dtype=numpy.float64) ))

    return #TODO


def enzyme_deletion_viable(**kwargs):
    connector = kwargs.get('connector', None)
    pair = kwargs.get('pair', None)

    allViable = False

    # connect to the DB and check in the PomEnzNet.Viability table if deletion of both enzymes is viable, or at least one is inviable
    # 'condition-dependent' will be considered as viable
    if connector is None or pair is None:
        print("Can't don anything without database or data on viability")
        return

    cursor = connector.connect_to_db()

    query = ""
    cursor.execute(query)
    via = cursor.fetchall()
    #for ...



    return allViable


def analyse_splitted_pair(enzyme1counts, enzyme2counts):
    sump = enzyme1counts.copy()
    for strain in enzyme2counts.keys():
            if strain in sump.keys():
                sump[strain]=sump[strain]+enzyme2counts[strain]

    pair_sd = numpy.std([numpy.log2(x) for x in list(sump.values()) if x > 0], dtype=numpy.float64)

    pair_gmean = gmean([numpy.log2(x) for x in list(sump.values()) if x > 0], dtype=numpy.float64)

    pair_correlation=0
    if len(list(enzyme1counts.values())) == len(list(enzyme2counts.values())):
        try:
            pair_correlation = scipy.stats.pearsonr(list(enzyme1counts.values()), list(enzyme2counts.values()))
        except FloatingPointError as err:
            print("Error in calculating pearson correlation")

    return pair_gmean, pair_sd, pair_correlation


def multiple_qtl_analysis():

    # TODO
    return

