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

##########################
# arguments
organism = "pombe"  # argv 1 (or cerevisiae)
data_type = "RNA"  # arv 2 (or protein)
enzyme_classification = "RMS"  # argv 3 (or "EC" or "all")
height = 2  # argv 4 ( -1 if EC and default 2 if "all")



def main():
    numpy.seterr(all='raise')

    with_backup_only = False

    organism = sys.argv[1]
    data_type = sys.argv[2]
    enzyme_classification = sys.argv[3]
    height = sys.argv[4]
    directory = ""

    if organism=="pombe":
        database = "PomEnzNet"
        if data_type == "RNA":
            directory = "/data/user/msorokin/data/DB/PomEnzNet/RMS_buffering/plots/pairEnz/"
        elif data_type == "protein":
            directory = "/data/user/msorokin/data/DB/PomEnzNet/RMS_buffering_prot/plots/pairEnz/"
    elif organism=="cerevisiae":
        database = "CerEnzNet"
        if data_type == "RNA":
            directory = "/data/user/msorokin/data/DB/CerEnzNet/RMS_buffering/plots/pairEnz/"
        elif data_type == "protein":
            directory = "/data/user/msorokin/data/DB/CerEnzNet/RMS_buffering_prot/plots/pairEnz/"

    set_global_variables(host, user, password, database, directory, organism, data_type, enzyme_classification, height)
    global_connector = MySQLConnector.MySQLConnector(host, user, password, database)

    ####################################################################################################################

    if organism == "pombe":
        # do for S. pombe
        if data_type == "RNA":
            # do RNA
            all_pairs_scores = retrieve_pair_enzyme_data_for_tests_pombe(enzymeClassification=enzyme_classification, height=height, with_backup_only=with_backup_only)

            outfile = "/data/user/msorokin/data/DB/PomEnzNet/data/allEnzPairDataRNA.DB"
            print_all_pairs_data(all_pairs_scores, outfile)

        elif data_type == "protein":
            # do protein
            all_pairs_scores = retrieve_pair_enzyme_data_for_tests_pombe_prot(enzymeClassification=enzyme_classification, height=height)
            outfile = "/data/user/msorokin/data/DB/PomEnzNet/data/allEnzPairDataProtein.DB"
            print_all_pairs_data(all_pairs_scores, outfile)

    elif organism == "cerevisiae":
        # do for S. cerevisiae
        if data_type == "RNA":
            print(database)
            # do RNA
            all_pairs_scores = retrieve_pair_enzyme_data_for_tests_cerevisiae(enzymeClassification=enzyme_classification, height=height)
            outfile = "/data/user/msorokin/data/DB/CerEnzNet/data/allEnzPairDataRNA.DB"
            print_all_pairs_data_cer(all_pairs_scores, outfile)
        elif data_type == "protein":
            # do protein
            all_pairs_scores = retrieve_pair_enzyme_data_for_tests_cerevisiae_prot(enzymeClassification=enzyme_classification, height=height)
            outfile = "/data/user/msorokin/data/DB/CerEnzNet/data/allEnzPairDataProtein.DB"
            print_all_pairs_data_cer(all_pairs_scores, outfile)
    else:
        print("Wrong organism: "+organism)



if __name__ == '__main__':
    main()
