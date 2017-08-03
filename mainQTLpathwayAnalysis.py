#!/usr/bin/python
from __future__ import print_function

from RMSbufferingMethods.dataImportMethods import *
from RMSbufferingMethods.dataAnalysisMethods import *
from RMSbufferingMethods.dataOutputMethods import *
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
from collections import defaultdict
from scipy.stats.stats import pearsonr
import operator
import matplotlib.pyplot as plt
from scipy.stats import *

##########################
# DB info

host = "localhost"
user = "readonly"
password = "R38d!0nly"
database = "PomEnzNet"


def main():
    numpy.seterr(all='raise')

    organism = sys.argv[1]
    data_type = sys.argv[2]

    correlationN = "<=0"
    correlationH = ">=0.3"
    sdReductionN = "<=0.1"
    sdReductionH = "<=0.1"
    exp_condition = "N"

    global_connector = MySQLConnector.MySQLConnector(host, user, password, database)

    # retrieve RMS pairs in same RMS and the list of QTL markers they are associated to
    pairs_with_same_qtl = retrieve_pairs_with_qtl(connector=global_connector, data_type=data_type, #correlationN=correlationN,
                                                   #correlationH=correlationH,
                                                  #sdReductionN=sdReductionN, sdReductionH=sdReductionH,
                                                  exp_condition=exp_condition)
    print(pairs_with_same_qtl)
    # in this table, the key is the pair of enzymes, and the associated list is the qtls
    # for each pair, redo the analysis (correlations, etc) by separating the strains on QTL markers
    # if only one marker - separate only on this marker
    # if several markers - OR separate on combination of markers, OR do the same analysis on each marker

    pairs_of_interest = {}

    for pair in pairs_with_same_qtl.keys():
        if len(pairs_with_same_qtl[pair]) >= 1:
            #print(pair)
            single_qtl_analysis(exp_condition=exp_condition, connector=global_connector, data_type=data_type, pair=pair, marker=pairs_with_same_qtl[pair][0])

            #network_pair_qtl_analysis(connector=global_connector, data_type=data_type, exp_condition=exp_condition, pair=pair)
        if len(pairs_with_same_qtl[pair]) > 1:
            single_qtl_analysis(exp_condition=exp_condition, connector=global_connector, data_type=data_type, pair=pair,
                                marker=pairs_with_same_qtl[pair][1])
            multiple_qtl_analysis()



if __name__ == '__main__':
    main()
