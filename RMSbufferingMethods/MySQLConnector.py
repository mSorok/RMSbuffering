#!/usr/bin/python
from __future__ import print_function


import os
import sys
import mysql.connector


class MySQLConnector:
    def __init__(self, host, user, password, database):
        self.host = host
        self.user = user
        self.password = password
        self.database = database
        self.conn = None

    def connect_to_db(self):
        self.conn = mysql.connector.connect(host=self.host, user=self.user, password=self.password, database=self.database)
        cursor = self.conn.cursor()
        return cursor

    def disconnect_from_db(self):
        self.conn.close()
        return






