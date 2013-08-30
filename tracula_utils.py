def modify_tracula_cfg_file(fsHome, subjectsDir, \
                            baseCfgFN, dtRoot, subjID, dcmRoot, dcmFN, \
                            bvalsFN, bvecsFN, nb0, doEddy, doRotVecs, thrBet, \
                            outCfgFN):    
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
