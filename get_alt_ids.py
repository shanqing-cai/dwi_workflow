
def get_alt_study_ids(pnTable, t_studyID):
    altStudyIDs = []

    for (i0, t_row) in enumerate(pnTable):
        if t_row[0] == t_studyID:
            for (j0, item) in enumerate(t_row[1:]):
                altStudyIDs.append(item)

    return altStudyIDs



def get_alt_ids(t_studyID, t_subjID):
    pnTable = [["innerspeech", "InSp"]]

    altIDs = []

    altIDs.append(t_studyID + t_subjID)
    altIDs.append(t_studyID + "_" + t_subjID)

    if t_subjID.startswith("S") and t_subjID[1:].isdigit():
        bs = True
        altIDs.append(t_studyID + t_subjID[1:])
        altIDs.append(t_studyID + "_" + t_subjID[1:])
    else:
        bs = False

    if len(t_subjID) > 6 and t_subjID.startswith("SEQ") and \
       (t_subjID[5] == "P" or t_subjID[5] == "C"):
        altIDs.append(t_subjID.replace("SEQ", "SEQPDS"))

    if len(t_subjID) > 4 and t_subjID.startswith("FRS") and \
       t_subjID[3:].isdigit():
        altIDs.append(t_subjID)
        
    if len(t_subjID) > 5 and t_subjID.startswith("CCRS") and \
       t_subjID[4:].isdigit():
        altIDs.append(t_subjID)

    if len(t_subjID) > 4 and t_subjID.startswith("SEQ") and \
       t_subjID[3:].isdigit():
        altIDs.append(t_subjID)

    #=== TODO: SDAP subjects ===#
    

    #=== Alternative study IDs ===#
    altStudyIDs = get_alt_study_ids(pnTable, t_studyID)
    if len(altStudyIDs) > 0:
        for (i2, t_altStudyID) in enumerate(altStudyIDs):
            altIDs.append(t_altStudyID + t_subjID)
            altIDs.append(t_altStudyID + "_" + t_subjID)
        
        if bs:
            altIDs.append(t_altStudyID + t_subjID[1:])
            altIDs.append(t_altStudyID + "_" + t_subjID[1:])
            
    return altIDs
