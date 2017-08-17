#!/usr/bin/python
from __future__ import print_function

from RMSbufferingMethods.dataImportMethods import *
from RMSbufferingMethods.dataAnalysisMethods import *
from RMSbufferingMethods.dataOutputMethods import *
from RMSbufferingMethods import MySQLConnector

import sys
import numpy

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

    if organism=="pombe":
        database = "PomEnzNet"
        if data_type == "RNA":
            outdir = "/data/user/msorokin/data/QTLmappings/data/indata/pom/"
        elif data_type == "protein":
            outfile = "~/data/..."
    elif organism == "cerevisiae":
        database = "CerEnzNet"
        if data_type == "RNA":
            outfile = "/data/user/msorokin/data/QTLmappings/data/indata/cer/pair_phenotype_RNA.tsv"
        elif data_type == "protein":
            outfile = "~/data/..."
    set_global_variables(host, user, password, database, directory, organism, data_type, enzyme_classification, height)
    global_connector = MySQLConnector.MySQLConnector(host, user, password, database)

    ####################################################################################################################

    if organism == "pombe":
        # do for S. pombe
        if data_type == "RNA":
            print(database)
            #file name = pair_phenotype_RNA_N.tsv
            prepare_data_for_pairQTL_pombe(connector=global_connector, outdir=outdir)
        elif data_type == "protein":
            print(database)
    if organism == "cerevisiae":

        if data_type == "RNA":
            print(database)
        elif data_type == "protein":
            print(database)



if __name__ == '__main__':
    main()

