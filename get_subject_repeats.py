#!/usr/bin/env python

SQL_SETTINGS_FN = "mega_sql_settings"
TOTAL_SCANS_TABLE_NAME = "Total Scans"

if __name__ == "__main__":
    import os
    import sys
    import argparse
    from scai_utils import info_log, check_file

    ap = argparse.ArgumentParser(description="Find subjects who participated in multiple studies")

    #--- Read sql settings ---#
    sqlSettingsFN = os.path.join(os.getenv("HOME"), SQL_SETTINGS_FN)
    check_file(sqlSettingsFN)
    sf = open(sqlSettingsFN, "rt")
    settings = sf.read().replace("\r", "\n").split("\n")
    sf.close()

    sqlSettings = {"SQL_SERVER": settings[0],
                   "DATABASE_NAME": settings[1], 
                   "SQL_USER": settings[2], 
                   "pw": settings[3]}
    
    #--- Get data from server ---#
    import MySQLdb

    db = MySQLdb.connect(host=sqlSettings["SQL_SERVER"],
                         db=sqlSettings["DATABASE_NAME"],
                         user=sqlSettings["SQL_USER"], 
                         passwd=sqlSettings["pw"])
    cursor = db.cursor()
    cursor.execute("""SELECT `Study Code`, `Lab Code` FROM `Master Code`""")
    qRes = cursor.fetchall()

    #print(qRes) # DEBUG

    subjs = {}
    for (i0, entry) in enumerate(qRes):
        sID = entry[0]
        labCode = entry[1]

        if not (labCode in subjs):
            subjs[labCode] = []

        subjs[labCode].append(sID)

    for (i0, t_labCode) in enumerate(subjs):
        if len(subjs[t_labCode]) > 1:
            sIDs = "%s\t" % t_labCode
            for t_sID in subjs[t_labCode]:
                sIDs += t_sID + " | "
            sIDs = sIDs[:-3]
            print(sIDs)
            

    
    
