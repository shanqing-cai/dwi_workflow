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
    
    # Removed: mean_f1samples, mean_S0samples
    
    r = True
    r = r and os.path.isdir(bedpDir)

    for (i0, efn) in enumerate(expectFiles):
        r = r and os.path.isfile(os.path.join(bedpDir, efn))

        if not r:
            info_log("INFO: Cannot find expected file: %s" % efn, logFN=logFN)
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
                if t_item.startswith("lh.") or t_item.startswith("rh") and \
                   len(t_item) > 6:
                    paths.append(t_item.replace("(", "").replace(")", ""))
                    
    tracCfg["paths"] = paths
    #print(txt)
    
    return tracCfg
