#!/usr/bin/env python

SQL_SETTINGS_FN = "mega_sql_settings"

TOTAL_SCANS_TABLE_NAME = "Total Scans"
MASTER_CODE_PREFIX = "SL"
STUDY_CODE_FIELD_NAME = "Study Code"
MASTER_CODE_FIELD_NAME = "Lab Code"

#=== A list of synonymous field names ===#
# First column is the field name in the SQL database
# The first column contains all valid field names
measNames = [["age"], 
             ["Sex", "gender"], 
             ["B-Day", "DOB", "DoB"], 
             ["Race"],
             ["Control or Patient ID", "ctrlPtID", "grpID"], 
             ["Site"], 
             ["Date Enrolled", "dateEnrolled"], 
             ["Dates Scanned", "dateScanned"]]

def print_help():
    print("Usage")
    print(__file__ + " dataFld --id projID1,projID2,... subjID1,subjID2,...")
    print("\nList of valid values for dataFld")
    
    for (i0, t_row) in enumerate(measNames):
        t_str = "\t"
        for (i1, t_item) in enumerate(t_row):
            if t_item.count(" ") > 0:
                t_str += "`%s`" % t_item
            else:
                t_str += t_item 

            if i1 < len(t_row) - 1:
                t_str += " | "

        print(t_str)

if __name__ == "__main__":
    import os
    import sys
    from scai_utils import *

    bVerbose = (sys.argv.count("-v") > 0)

    #=== Print help message ===#
    if len(sys.argv) == 1 or \
       sys.argv.count("-h") > 0 or sys.argv.count("--help") > 0:
        print_help()
        sys.exit(0)

    #=== Determine field name ===#
    ifld = sys.argv[1]

    bFound = False
    for (i0, t_row) in enumerate(measNames):
        for (i1, t_item) in enumerate(t_row):
            if ifld.lower() == t_item.lower():
                bFound = True
                fld = t_row[0]
                break;
    
    if not bFound:
        raise Exception, "Unrecognized data field name: %s" % ifld


    #=== Determine subject info mode (--mc or --id) ===#
    if sys.argv[2] == "--mc":
        sMode = "mc" # Master code (e.g., --mc 1,165)
        mcs = sys.argv[3].split(",")
    elif sys.argv[2] == "--id":
        sMode = "id" # Study and subject IDs (e.g., --id CAT,RHY 12,AWS_M01)
        studyIDs = sys.argv[3].split(",")
        subjIDs = sys.argv[4].split(",")

        if len(studyIDs) != len(subjIDs):
            raise Exception, "Unequal number of entries in studyIDs and subjIDs"
    else:
        raise Exception, "Unrecognized subject info mode: %s" % sys.argv[2]

    
    #=== Establish SQL server connection ===#
    import MySQLdb

    sqlSettingsFN = os.path.join(os.getenv("HOME"), SQL_SETTINGS_FN)
    check_file(sqlSettingsFN)
    sf = open(sqlSettingsFN, "rt")
    settings = sf.read().replace("\r", "\n").split("\n")
    sf.close()

    SQL_SERVER = settings[0]
    DATABASE_NAME = settings[1]
    SQL_USER = settings[2]
    pw = settings[3]

    if bVerbose:
        info_log("SQL settings:")
        info_log("\tSQL_SERVER = %s" % SQL_SERVER)
        info_log("\tDATABASE_NAME = %s" % DATABASE_NAME)
        info_log("\tSQL_USER = %s" % SQL_USER)

    db = MySQLdb.connect(host=SQL_SERVER, db=DATABASE_NAME,
                         user=SQL_USER, passwd=pw)
    if bVerbose:
        info_log("Connectoin to database %s at %s under user name %s has been established successfully." % (DATABASE_NAME, SQL_SERVER, SQL_USER))

    #=== Query SQL ===#
    cursor = db.cursor()

    if sMode == "id":
        sqlCmd = """SELECT `%s`, `%s` FROM `%s`""" % \
                 (STUDY_CODE_FIELD_NAME, fld, \
                  TOTAL_SCANS_TABLE_NAME)
        #print(sqlCmd)
        cursor.execute(sqlCmd)

    qRes = cursor.fetchall()
    #print(qRes)
    
    #=== SQL clean up ===#
    db.close()
    if bVerbose:
        info_log("Database connect closed.")

    #=== Extract data from query result ===#
    from get_alt_ids import get_alt_ids

    res = []
    if sMode == "id":
        for (i0, t_studyID) in enumerate(studyIDs):
            t_subjID = subjIDs[i0]

            altIDs = get_alt_ids(t_studyID, t_subjID)

            bFound = False
            for (i1, t_row) in enumerate(qRes):
                if t_row[0] == None or t_row[1] == None:
                    continue

                for (i2, t_altID) in enumerate(altIDs):
                    if t_row[0] == t_altID:
                        bFound = True
                        foundRow = t_row
                        break
                
                if bFound:
                    break
            
            if not bFound:
                res.append("")
                continue

            res.append(foundRow[1])
            
    for (i0, t_res) in enumerate(res):
        print(t_res)
