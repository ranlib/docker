#!/usr/bin/env python
"""
Retrieves data from JGI production databases.
"""

import sys
import base64  # encrypt/decrypt mysql pwd
import cx_Oracle as oracle
import MySQLdb as mySql
import db_access

def connectDb(dbName=None):
    """
    Connects to a specified db and returns a connection object.

    @param dbName string 'RQC' or 'ITS'
    @return obj Databse connection object
    """

    dbName = dbName.upper()
    #print("<I> connectDb: db name = " + dbName )
    if dbName == "ITS":
        db = db_access.jgi_connect_db("dw-pg")
    elif dbName == "RQC":
        #print("<I> connectDb: db name = " + dbName )
        db = db_access.jgi_connect_db("rqc")
    else:
        print("<E> connectDb: unknown db name" + dbName)
        sys.exit(1)
    return db


def queryDb(dbName=None, sql=None, sqlParams=None):
    """
    Given a query and optional parameters, queries specified database,
    then returns an iterable for retrieved data.

    @param dbName string 'RQC' or 'ITS'
    @param sql string SQL statement. See cx_Oracle or MySQLdb for interpolation style
    @param sqlParams string Optionally provide SQL interpolation params
    @return iter Returns a tuple of tuples containing data
    """
    #print("<I> queryDB: db = " + dbName)
    #db_access.jgi_connect_db("list")

    connection = connectDb(dbName=dbName)
    cursor = connection.cursor()

    if sqlParams:
        cursor.execute(sql, sqlParams)
    else:
        cursor.execute(
            sql,
        )
    return cursor


# test sql connection
def testSql():

    sql = """
SELECT DISTINCT
 datawh.seq_unit_name,
 datawh.sample_id,
 datawh.sample_name,
 ls.index_name,
 ls.index_sequence,
 datawh.sample_tube_plate_label,
 datawh.plate_location,
 datawh.target_fragment_size_bp,
 datawh.actual_run_type,
 datawh.amplified,
 datawh.lib_name
FROM dw.all_inclusive_report datawh
INNER JOIN dw.library_stock ls ON datawh.lib_name = ls.library_name
WHERE
  datawh.lib_name IN (:lib1)"""

    params = ["NAZS"]

    # Bryce - 2016-01-27: the field names changed and I didn't want to re-look them all up
    # print "ITS query with params"
    # for row in queryDb(dbName = 'ITS', sql = sql, sqlParams=params):
    #    print row[0], row[1]

    sql = "select acct_scientific_program, lib_name from DW.all_inclusive_report where lib_name in ('HOOO', 'NAZS')"

    print("ITS query without params:")
    for row in queryDb(dbName="ITS", sql=sql):
        print(row[0], row[1])

    print("RQC query without params:")
    sql = "select library_name, seq_proj_name from library_info li where li.library_name in ('PNWN','PNWB')"

    for row in queryDb(dbName="RQC", sql=sql):
        print(row[0], row[1])

    print("RQC query with params:")

    params = ("HOOO", "NPTS")
    sql = "select l.library_name, l.seq_proj_name from library_info l where l.library_name in (%s, %s)"

    for row in queryDb(dbName="RQC", sql=sql, sqlParams=params):
        print(row[0], row[1])

    print("SDM query:")
    sql = "select su.sdm_seq_unit_id, suf.file_name, su.library_name, suf.dt_modified from sdm_seq_unit_file suf inner join sdm_seq_unit su on suf.sdm_seq_unit_id = su.sdm_seq_unit_id where suf.file_type = 'FASTQ' order by su.sdm_seq_unit_id desc limit 0,10"
    for row in queryDb(dbName="sdm", sql=sql):
        print(row[0], row[1], row[2], row[3])


if __name__ == "__main__":

    try:
        testSql()
    except KeyboardInterrupt:
        sys.exit(0)
    except Exception as err:
        sys.stderr.write("ERROR: %s\n" % str(err))
        sys.exit(1)
