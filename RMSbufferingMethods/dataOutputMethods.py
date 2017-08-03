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

import rpy2
from rpy2 import robjects
from rpy2.robjects import Formula, Environment
from rpy2.robjects.vectors import IntVector, FloatVector
from rpy2.robjects.lib import grid
from rpy2.rinterface import RRuntimeError
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
import warnings


def print_all_pairs_data(biglist, outfile):
    # biglist[pair]
    # 0: gmeain_sum_n
    # 1: pais_sd_n
    # 2 : correlation_n
    # 3 : gmean_sum_h
    # 4 : sd_sum_h
    # 5 : correlation_h
    # 6 : delat_sd
    # 7 : sd_reduction_score_min_n
    # 8 : sd_reduction_score_min_h
    # 9 : sd_fc_pair
    # 10 : sd_reduction_score_fc
    # 11 : mean_fc_pair
    # 12 : in same group
    # 13 : data type ("RNA" or "protein")

    fo = open(outfile, "w")

    for pair in biglist.keys():
        corr_n = 0
        corr_h = 0
        if isinstance(biglist[pair][2], tuple):# and biglist[pair][2][1] <= 0.05:
            corr_n = biglist[pair][2][0]

        if isinstance(biglist[pair][5], tuple): # and biglist[pair][5][1] <= 0.05:
            corr_h = biglist[pair][5][0]

        outstring = pair[0] + "\t" + pair[1] + "\t" + str(biglist[pair][0]) + "\t" + str(
            biglist[pair][1]) + "\t" + str(corr_n) + "\t" + str(biglist[pair][3]) + "\t" + str(
            biglist[pair][4]) + "\t" + str(corr_h) + "\t" + str(biglist[pair][6]) + "\t" + str(
            biglist[pair][7]) + "\t" + str(biglist[pair][8]) + "\t" + str(
            biglist[pair][9]) + "\t" + str(biglist[pair][10]) + "\t" + str(
            biglist[pair][11]) + "\t" + str(biglist[pair][12]) + "\t" + str(
            biglist[pair][13]) + "\n"
        fo.write(outstring)

    fo.close()
    return


def print_all_pairs_data_cer(biglist, outfile):
    # biglist[pair]
    # 0: gmeain_sum
    # 1: pais_sd
    # 2 : correlation
    # 3 : sd_reduction_score_min
    # 4 : in same group
    # 5 : data type ("RNA" or "protein")

    fo = open(outfile, "w")

    for pair in biglist.keys():
        corr = 0
        if isinstance(biglist[pair][2], tuple):  # and biglist[pair][2][1] <= 0.05:
            corr = biglist[pair][2][0]

        outstring = pair[0] + "\t" + pair[1] + "\t" + str(biglist[pair][0]) + "\t" + str(
            biglist[pair][1]) + "\t" + str(corr) + "\t" + str(biglist[pair][3]) + "\t" + str(
            biglist[pair][4]) + "\t" + str(biglist[pair][5]) + "\n"
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
