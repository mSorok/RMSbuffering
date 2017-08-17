#!/usr/bin/python
from __future__ import print_function
from RMSbufferingMethods import dataAnalysisMethods
from RMSbufferingMethods import dataImportMethods

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
from scipy.stats.mstats import gmean
from scipy.stats import ttest_ind
import operator
import matplotlib.pyplot as plt
import matplotlib.patches as mpatch

import warnings


def print_all_pairs_data(biglist, outfile):
    # biglist[pair]
    # 0: pair_sd_N
    # 1: skew_N
    # 2: kurtosis_N
    # 3 : correlation_N

    # 4 : sd_sum_H
    # 5 : skew_H
    # 6 : kurtosis_H
    # 7 : correlation_H

    # 8 : delta_sd

    # 9 : sd_reduction_score_min_N
    # 10 : sd_reduction_score_min_H

    # 11 : sd_fc_pair
    # 12 : sd_reduction_score_fc
    # 13 : mean_fc_pair

    # 14 : skew_reduction_score_min_N
    # 15 : skew_reduction_score_min_H
    # 16 : kurtosis_reduction_score_min_N
    # 17 : kurtosis_reduction_score_min_H

    # 18 : sdreduction_skew_N
    # 19 : sdreduction_skew_H
    # 20 : sdReduction_kurtosis_N
    # 21 : sdReduction_kurosis_H

    # 22 : p-val pitman_morgan N
    # 23 : p-val pitman_morgan H
    # 24 : in same group
    # 25 : data type ("RNA" or "protein")

    fo = open(outfile, "w")

    for pair in biglist.keys():
        corr_n = 0
        corr_h = 0
        if isinstance(biglist[pair][3], tuple):# and biglist[pair][2][1] <= 0.05:
            corr_n = biglist[pair][3][0]

        if isinstance(biglist[pair][7], tuple): # and biglist[pair][5][1] <= 0.05:
            corr_h = biglist[pair][7][0]

        outstring = pair[0] + "\t" + pair[1] + "\t" + str(biglist[pair][0]) + "\t" + str(
            biglist[pair][1]) + "\t" + str(biglist[pair][2]) + "\t" + str(corr_n) + "\t" + str(biglist[pair][4]) + "\t" + str(
            biglist[pair][5]) + "\t" + str(biglist[pair][6]) + "\t" + str(corr_h) + "\t" + str(biglist[pair][8]) + "\t" + str(
            biglist[pair][9]) + "\t" + str(biglist[pair][10]) + "\t" + str(
            biglist[pair][11]) + "\t" + str(biglist[pair][12]) + "\t" + str(
            biglist[pair][13]) + "\t" + str(biglist[pair][14]) + "\t" + str(biglist[pair][15]) + "\t" + str(
            biglist[pair][16]) + "\t" + str(biglist[pair][17]) + "\t" + str(biglist[pair][18]) + "\t" + str(
            biglist[pair][19]) + "\t" + str(biglist[pair][20]) + "\t" + str(biglist[pair][21]) + "\t" + str(
            biglist[pair][22]) + "\t" + str(biglist[pair][23]) + "\t" + str(
            biglist[pair][24]) + "\t" + str(biglist[pair][25]) + "\n"
        fo.write(outstring)

    fo.close()
    return


def print_all_pairs_data_cer(biglist, outfile):
    # biglist[pair]
    # 0 : pair_sd
    # 1 : correlation
    # 2 : sd_reduction_score_min
    # 3 : skew_reduction_score
    # 4 : kurtosis_reduction_score
    # 5 : sd_reduction_skew corrected
    # 6 : sd_reduction_kurtosis corrected
    # 7 : pair_skew
    # 8 : pair_kurtosis
    # 9 : pvalue pitman morgan
    # 10 : in same group
    # 11 : data type ("RNA" or "protein")

    fo = open(outfile, "w")

    for pair in biglist.keys():
        corr = 0
        if isinstance(biglist[pair][1], tuple):  # and biglist[pair][2][1] <= 0.05:
            corr = biglist[pair][1][0]

        outstring = pair[0] + "\t" + pair[1] + "\t" + str(biglist[pair][0]) + "\t" + str(corr) + "\t" + str(biglist[pair][2]) + "\t" + str(biglist[pair][3]) + "\t" + str(
            biglist[pair][4]) + "\t" + str(biglist[pair][5]) + "\t" + str(biglist[pair][6])+ "\t" + str(
            biglist[pair][7]) + "\t" + str(biglist[pair][8]) + "\t" + str(biglist[pair][9])+ "\t" + str(biglist[pair][10])+ "\t" + str(biglist[pair][11]) + "\n"
        fo.write(outstring)

    fo.close()
    return




def plot_hists_ec_analysis(list1, list2, label1, label2, title, xlab):
    plt.clf()
    numpy.seterr(all='ignore', divide='ignore', over='ignore', under='ignore', invalid='ignore')
    directory = "/data/user/msorokin/data/DB/PomEnzNet/RMS_buffering/withoutEC/"

    plt.hist(list1, bins=30, normed=True, color='b', alpha=0.5, label=label1)
    plt.hist(list2, bins=30, normed=True, color='r', alpha=0.5, label=label2)
    axes = plt.gca()
    plt.title(title)
    plt.xlabel(xlab)
    plt.ylabel("Frequency")

    # t-test
    res = ttest_ind(list1, list2)

    ymin, ymax = axes.get_ylim()
    xmin, xmax = axes.get_xlim()
    textstr = '$t-test$\n$statistic=%.2f$\n$p-value=%.2e$' % (res[0], res[1])
    plt.text(xmin, ymax, textstr, bbox={'facecolor': 'white', 'alpha': 1, 'pad': 10}, ha='center',
             va='bottom')

    plt.legend(fontsize='x-small', loc='upper right')
    fname = title.replace(" ", "_")
    fname = directory + fname + ".png"
    plt.savefig(filename=fname, format="png", orientation='landscape', dpi=150)
    #plt.show()
    return


########################################################################################################################


def prepare_data_for_pairQTL_pombe(**kwargs):
    outdir = kwargs.get('outdir', "~/data/")
    connector = kwargs.get('connector', None)

    # retrieve counts
    # sum enzyme pairs with backup
    # arrange them by strain
    # print out the right way in the file

    if connector is None:
        print("Can't connect to the database!")
        return

    cursor = connector.connect_to_db()

    # First for N condition
    countsN = {}
    strainsN = set()

    queryN = "SELECT CONCAT(EnzymeNode_id_A, '$',EnzymeNode_id_B), r1.strain, r1.count AS count_A, r2.count AS count_B " \
             "FROM PairData " \
             "INNER JOIN RNAcounts r1 ON(EnzymeNode_id_A=r1.EnzymeNode_id) " \
             "INNER JOIN RNAcounts r2 ON(EnzymeNode_id_B=r2.EnzymeNode_id)  " \
             "WHERE inSameGroup='True' AND dataType='RNA' AND r1.strain=r2.strain AND r1.exp_condition='N' AND r2.exp_condition='N';"

    cursor.execute(queryN)
    counts = cursor.fetchall()
    for c in counts:
        # c0=trait id, c1, strain, c2 count1, c3 count2
        strainsN.add(c[1])
        if c[0] in countsN.keys():
            countsN[c[0]][c[1]] = c[2]+c[3]
        else:
            countsN[c[0]] = {}
            countsN[c[0]][c[1]] = c[2] + c[3]

    outfileN = "pair_phenotype_RNA_N.tsv"
    foutn = open(outdir+outfileN, 'w')

    liststrains = list(strainsN)
    firstline = "\t".join(liststrains)
    foutn.write(firstline+"\n")

    for trait in countsN.keys():
        line = str(trait)
        for strain in liststrains:
            if strain in countsN[trait] :
                line = line+"\t"+str(countsN[trait][strain])
            else:
                line = line+"\tNA"
        line = line+"\n"
        foutn.write(line)

    foutn.close()

    # Next for H condition

    countsH = {}
    strainsH = set()

    queryH = "SELECT CONCAT(EnzymeNode_id_A, '$',EnzymeNode_id_B), r1.strain, r1.count AS count_A, r2.count AS count_B " \
             "FROM PairData " \
             "INNER JOIN RNAcounts r1 ON(EnzymeNode_id_A=r1.EnzymeNode_id) " \
             "INNER JOIN RNAcounts r2 ON(EnzymeNode_id_B=r2.EnzymeNode_id)  " \
             "WHERE inSameGroup='True' AND dataType='RNA' AND r1.strain=r2.strain AND r1.exp_condition='H' AND r2.exp_condition='H';"


    cursor.execute(queryH)
    counts = cursor.fetchall()
    for c in counts:
        # c0=trait id, c1, strain, c2 count1, c3 count2
        strainsH.add(c[1])
        if c[0] in countsH.keys():
            countsH[c[0]][c[1]] = c[2] + c[3]
        else:
            countsH[c[0]] = {}
            countsH[c[0]][c[1]] = c[2] + c[3]

    outfileH = "pair_phenotype_RNA_H.tsv"
    foutn = open(outdir + outfileH, 'w')

    liststrains = list(strainsH)
    firstline = "\t".join(liststrains)
    foutn.write(firstline + "\n")

    for trait in countsH.keys():
        line = str(trait)
        for strain in liststrains:
            if strain in countsH[trait]:
                line = line + "\t" + str(countsH[trait][strain])
            else:
                line = line + "\tNA"
        line = line + "\n"
        foutn.write(line)
    foutn.close()
    connector.disconnect_from_db()
    return
