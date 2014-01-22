#!/usr/bin/env python

SQL_SERVER = "two.bu.edu"
DATABASE_NAME = "MEGA"
SQL_USER = "scai"
PASSWD_FN = "passwd"



if __name__ == "__main__":
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
        raise Exception, "Unequal number of entries in subjIDs and subjIDs"

    if args.bVerbose:
        info_log("# of subjects queried = %d" % len(subjIDs))

    #=== Establish SQL server connection ===#
    import MySQLdb

    check_file(PASSWD_FN)
    pwf = open(PASSWD_FN, "rt")
    pw = pwf.read()
    pwf.close()

    db = MySQLdb.connect(host=SQL_SERVER, db=DATABASE_NAME,
                         user=SQL_USER, passwd=pw)
    if args.bVerbose:
        info_log("Connectoin to database %s at %s under user name %s has been established successfully." % (DATABASE_NAME, SQL_SERVER, SQL_USER))

    #=== SQL Query ===#
    cursor = db.cursor()
    cursor.execute("""SELECT `Study Code`, `Lab Code` FROM `Master Code`""")
    qRes = cursor.fetchall()

    #print(qRes)

    #=== SQL clean up ===#
    db.close()
    if args.bVerbose:
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
        
        if not bFound:
            masterCodes.append(-1)
            continue

        masterCodes.append(int(foundRow[1].replace("SL", "")))

    for (i0, t_mc) in enumerate(masterCodes):
        print(t_mc)

        
    

    
