def modify_tracula_cfg_file(fsHome, subjectsDir, \
                            baseCfgFN, dtRoot, subjID, dcmRoot, dcmFN, \
                            bvalsFN, bvecsFN, nb0, doEddy, doRotVecs, thrBet, \
                            outCfgFN, logFN=None):
    from scai_utils import info_log
    import os

    #=== Read and modify the base tracula cfg file ===#
    cfgf = open(baseCfgFN, "r")
    cfgt = cfgf.read().replace('\r', '\n').split('\n')
    cfgf.close()
        
    newcfgt = ""
    for (i0, t_line) in enumerate(cfgt):
        if len(t_line) == 0:
            newcfgt += t_line + '\n'
            continue
        else:
            t_items = t_line.split(' ')
            if len(t_items) >= 3:
                if t_items[1] == "dtroot":
                    t_line = t_line.replace(t_items[-1], dtRoot)
                elif t_items[1] == "subjlist":
                    t_line = t_line.replace(t_items[-1], \
                                            "(%s)" % subjID)
                elif t_items[1] == "dcmroot":
                    t_line = t_line.replace(t_items[-1], dcmRoot)
                elif t_items[1] == "dcmlist":
                    dcmlst = dcmFN
                    t_line = t_line.replace(t_items[-1], "(%s)" % dcmlst)
                elif t_items[1] == "bvalfile":
                    t_line = t_line.replace(t_items[-1], bvalsFN)
                elif t_items[1] == "bvecfile":
                    t_line = t_line.replace(t_items[-1], bvecsFN)
                elif t_items[1] == "nb0":
                    t_line = t_line.replace(t_items[-1], "%d" % nb0)
                elif t_items[1] == "doeddy":
                    t_line = t_line.replace(t_items[-1], \
                                           "%d" % doEddy)
                elif t_items[1] == "dorotvecs":
                    t_line = t_line.replace(t_items[-1], \
                                            "%d" % doRotVecs)
                elif t_items[1] == "thrbet":
                    t_line = t_line.replace(t_items[-1], \
                                           "%.2f" % thrBet)
                elif t_items[1] == "trainfile":
                    t_line = t_line.replace(t_items[-1], \
                                            os.path.join(fsHome, \
                                                         "trctrain", \
                                                         "trainlist.txt"))
                elif t_items[1] == "SUBJECTS_DIR":
                    t_line = t_line.replace(t_items[-1], subjectsDir)

        newcfgt += t_line + '\n'

    #=== Write to new cfg file ===#
    cfg_f = open(outCfgFN, "w")
    cfg_f.write(newcfgt)
    cfg_f.close()

expectFiles = ["dyads1_dispersion.nii.gz", "dyads1.nii.gz",
               "dyads2_dispersion.nii.gz", "dyads2.nii.gz",
               "mean_dsamples.nii.gz", "mean_f1samples.nii.gz",
               "mean_f2samples.nii.gz", 
               "mean_ph1samples.nii.gz", "mean_ph2samples.nii.gz", 
               "mean_th1samples.nii.gz", 
               "mean_th2samples.nii.gz", 
               "merged_f1samples.nii.gz", "merged_f2samples.nii.gz", 
               "merged_ph1samples.nii.gz", "merged_ph2samples.nii.gz", 
               "merged_th1samples.nii.gz", "merged_th2samples.nii.gz", 
               "nodif_brain_mask.nii.gz"]

def check_bedp_complete(bedpDir):
    import os
    from scai_utils import info_log
    
    #print("Checking bedp completion in directory: %s" % bedpDir)

    # Removed: mean_f1samples, mean_S0samples
    
    r = True
    r = r and os.path.isdir(bedpDir)

    for (i0, efn) in enumerate(expectFiles):
        r = r and os.path.isfile(os.path.join(bedpDir, efn))

        if not r:
            info_log("INFO: Cannot find expected file: %s" % efn)
            return r

    return r

def get_tracula_settings(cfgFN):
    tracCfg = {}
    
    from scai_utils import check_file, read_text_file

    check_file(cfgFN)
    txt = read_text_file(cfgFN)

    bFoundPathList = False

    for (i0, t_line) in enumerate(txt):
        if t_line.strip().startswith("set pathlist ="):
            bFoundPathList = True
            break

    paths = []
    if bFoundPathList:
        i1 = i0
        ended = False
        while not ended:
            items = txt[i1].split(" ")

            if len(items) == 0:
                i1 += 1
                continue

            if items[-1] == "\\":
                items = items[:-1]

            ended = txt[i1].strip().endswith(")")
            i1 += 1

            for t_item in items:
                t_item = t_item.replace("(", "").replace(")", "")

                if t_item.startswith("lh.") or t_item.startswith("rh") or \
                   t_item.startswith("fmajor") or t_item.startswith("fminor") \
                   and len(t_item) > 6:
                    paths.append(t_item)
                    
    tracCfg["paths"] = paths
    #print(txt)
    
    return tracCfg

# Subroutine for getting the data from a given track of a given subject
def get_dpath_track_data(trackDir):
    # List of data fields to extract
    dataFlds = ["Count", "Volume", 
                "Len_Min", "Len_Max", "Len_Avg", "Len_Center", 
                "AD_Avg", "AD_Avg_Weight", "AD_Avg_Center", 
                "RD_Avg", "RD_Avg_Weight", "RD_Avg_Center", 
                "MD_Avg", "MD_Avg_Weight", "MD_Avg_Center", 
                "FA_Avg", "FA_Avg_Weight", "FA_Avg_Center"]

    dat = {}

    import os
    from scai_utils import check_dir, check_file, read_text_file, info_log

    check_dir(trackDir)
    
    # Read from pathstats.overall.txt
    pso = os.path.join(trackDir, "pathstats.overall.txt")
    check_file(pso)

    txt = read_text_file(pso)

    for (i0, t_line) in enumerate(txt):
        t_line = t_line.strip()
        
        if t_line.startswith("#"):
            continue

        t_items = t_line.split(" ")
        if len(t_items) == 2:
            if dataFlds.count(t_items[0]) == 1:
                dat[t_items[0]] = float(t_items[1])
                
    # Check for completeness
    for (i0, t_fld) in enumerate(dataFlds):
        if not (t_fld in dat):
            info_log("WARNING: could not find data field %s in file %s" % \
                     (t_fld, pso), bWarn=True)

    # Read from pathstats.byvoxel.txt
    vso = os.path.join(trackDir, "pathstats.byvoxel.txt")
    check_file(vso)

    txt = read_text_file(vso)

    import numpy as np
    bvDat = {}
    cols = {}
    
    # Determine the number of data points
    npts = 0
    lidx0 = -1
    bFoundHeader = False
    for (i0, t_line) in enumerate(txt):
        t_line = t_line.strip()
        if len(t_line) == 0:
            continue
        
        if not (t_line.startswith("#") or t_line.startswith("x")):
            npts += 1
            if lidx0 == -1:
                lidx0 = i0
            
        if t_line.startswith("x"):
            bFoundHeader = True
            t_items = t_line.split(" ")
            cols["x"] = t_items.index("x")
            cols["y"] = t_items.index("y")
            cols["z"] = t_items.index("z")
            cols["AD"] = t_items.index("AD")
            cols["RD"] = t_items.index("RD")
            cols["MD"] = t_items.index("MD")
            cols["FA"] = t_items.index("FA")
            cols["AD_Avg"] = t_items.index("AD_Avg")
            cols["RD_Avg"] = t_items.index("RD_Avg")
            cols["MD_Avg"] = t_items.index("MD_Avg")
            cols["FA_Avg"] = t_items.index("FA_Avg")
            
    if not bFoundHeader:
        raise Exception, "Cannot find header column in file %s" % vso


    txt = txt[lidx0 : lidx0 + npts]
    #print(txt)

    # Allocate space
    bvDat["x"] = np.zeros(npts)
    bvDat["y"] = np.zeros(npts)
    bvDat["z"] = np.zeros(npts)
    bvDat["AD"] = np.zeros(npts)
    bvDat["RD"] = np.zeros(npts)
    bvDat["MD"] = np.zeros(npts)
    bvDat["FA"] = np.zeros(npts)
    bvDat["AD_Avg"] = np.zeros(npts)
    bvDat["RD_Avg"] = np.zeros(npts)
    bvDat["MD_Avg"] = np.zeros(npts)
    bvDat["FA_Avg"] = np.zeros(npts)

    keys = bvDat.keys()

    for (i0, t_line) in enumerate(txt):
        t_items = t_line.split(" ")
        
        for t_fld in keys:
            bvDat[t_fld][i0] = float(t_items[cols[t_fld]])

    dat["byVoxel"] = bvDat

    return dat
