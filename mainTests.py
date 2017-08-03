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
    global_connector = MySQLConnector.MySQLConnector(host, user, password, database)

    analyse_in_rms_not_in_ec(connector=global_connector, height=2, data_type='RNA')





if __name__ == '__main__':
    main()


