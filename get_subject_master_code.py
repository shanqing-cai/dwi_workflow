#!/usr/bin/env python

SQL_SETTINGS_FN = "mega_sql_settings"

def get_subject_master_code(sqlSettings, studyIDs, subjIDs, bVerbose):
    import MySQLdb

    SQL_SERVER = sqlSettings["SQL_SERVER"]
    DATABASE_NAME = sqlSettings["DATABASE_NAME"]
    SQL_USER = sqlSettings["SQL_USER"]
    pw = sqlSettings["pw"]

    if bVerbose:
        info_log("SQL settings:")
        info_log("\tSQL_SERVER = %s" % SQL_SERVER)
        info_log("\tDATABASE_NAME = %s" % DATABASE_NAME)
        info_log("\tSQL_USER = %s" % SQL_USER)

    db = MySQLdb.connect(host=SQL_SERVER, db=DATABASE_NAME,
                         user=SQL_USER, passwd=pw)
    if bVerbose:
        info_log("Connectoin to database %s at %s under user name %s has been established successfully." % (DATABASE_NAME, SQL_SERVER, SQL_USER))

    #=== SQL Query ===#
    cursor = db.cursor()
    cursor.execute("""SELECT `Study Code`, `Lab Code` FROM `Master Code`""")
    qRes = cursor.fetchall()

    #=== (Conditional) load the "Total Scans" table for SDAP ===#
    if studyIDs.count("SDAP") > 0:
        cursor.execute("""SELECT `Subject Data Directory Name`, `Lab Code` FROM `Total Scans`""")
        tsRes = cursor.fetchall()
        #print(tsRes)
        
    #=== SQL clean up ===#
    db.close()
    if bVerbose:
        info_log("Database connect closed.")
    
    #=== Get master codes ===#
    from get_alt_ids import get_alt_ids

    masterCodes = []
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

        if (not bFound) and t_studyID == "SDAP": # Try "Total Scans" entries
            for (i1, t_row) in enumerate(tsRes):
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
            masterCodes.append(-1)
            continue

        

        masterCodes.append(int(foundRow[1].replace("SL", "")))

    return masterCodes

if __name__ == "__main__":
    import os
    import sys
    import argparse
    from scai_utils import info_log, check_file

    ap = argparse.ArgumentParser()
    ap.add_argument("studyIDs", help="Study IDs (e.g., SEQPDS,innerspeech)")
    ap.add_argument("subjIDs", help="Subject IDs (e.g., SEQ03P13,S10)")
    ap.add_argument("-v", dest="bVerbose", action="store_true", 
                    help="Verbose mode")

    if len(sys.argv) == 1:
        ap.print_help()
        sys.exit(0)

    args = ap.parse_args()

    studyIDs = args.studyIDs.split(",")
    subjIDs = args.subjIDs.split(",")

    if len(studyIDs) != len(subjIDs):
        raise Exception, "Unequal number of entries in studyIDs and subjIDs"

    if args.bVerbose:
        info_log("# of subjects queried = %d" % len(subjIDs))

    #=== Establish SQL server connection ===#
    sqlSettingsFN = os.path.join(os.getenv("HOME"), SQL_SETTINGS_FN)
    check_file(sqlSettingsFN)
    sf = open(sqlSettingsFN, "rt")
    settings = sf.read().replace("\r", "\n").split("\n")
    sf.close()

    
    sqlSettings = {"SQL_SERVER": settings[0],
                   "DATABASE_NAME": settings[1], 
                   "SQL_USER": settings[2], 
                   "pw": settings[3]}
        
    masterCodes = get_subject_master_code(sqlSettings, studyIDs, subjIDs, 
                                          args.bVerbose)


    for (i0, t_mc) in enumerate(masterCodes):
        print(t_mc)

        
    

    
