DEBUG = False

def get_cort_rois(aROIs, lobe="all", bSpeech=False):
#if __name__ == "__main__":
#    lobe = "all"
#    bSpeech = "speech_2g_lh_0.01"
    import numpy as np

    aROIs = np.array(aROIs)
    nROIs = len(aROIs)

    roiNames = aROIs[:, 0]

    lobeNames = aROIs[:, 1]
    uLobeNames = np.unique(lobeNames)

    isSpeech = aROIs[:, 2]
    nSpeechROIs = len(np.nonzero(isSpeech == 'S')[0])
    

    if lobe == "all":
        idxLobe = np.array(range(nROIs))
    else:
        if (len(np.nonzero(uLobeNames == lobe)[0])) == 0:
            idxLobe = np.array([])
            raise Exception, "Unrecognized lobe name: %s" % lobe
        else:
            idxLobe = np.nonzero(lobeNames == lobe)[0]
            # ROIs = roiNames[idx]

    if bSpeech == True:
        idxSpeech = np.nonzero(isSpeech == "S")[0]
    else:
        idxSpeech = np.array(range(nROIs))
       
    idx = np.intersect1d(idxLobe, idxSpeech)
    
    if DEBUG:
        print(idxLobe)
        print(idxSpeech)
        print(idx)
    ROIs = list(roiNames[idx])

    return ROIs


def get_roi_num(ctab_fn, roiHemi, roiName):
    from mri_utils import read_ctab
    import os

    # Read color table
    if not os.path.isfile(ctab_fn):
        raise Exception, "Cannot open color table (ctab): %s"%(ctab_fn)
    (ct_nums, ct_names) = read_ctab(ctab_fn)

    if not ct_names.count(roiName) == 1:
        return -1
#        raise Exception, "Not exactly one entry for ROI %s is found in%s"\
#                         %(roiName, ctab_fn)

    roi_num = ct_nums[ct_names.index(roiName)]

    if roiHemi == "lh":
        roi_num += 1000
    elif roiHemi == "rh":
        roi_num += 2000
    else:
        raise Exception, "Unrecognized hemisphere name: %s"%roiHemi

    return roi_num


def get_list_roi_nums(roiList, HEMIS, ctab):
#from aparc_utils import get_cort_rois, get_roi_num

    #= Get the number of all GM ROIs =#
    rois = []
    roiNums = []
    for (i0, hemi) in enumerate(HEMIS):
        for (i1, troi) in enumerate(roiList):
            assert(len(troi) == 4)
            tname = troi[0]
                
            rois.append("%s_%s" % (hemi, tname))
            
            tnum = get_roi_num(ctab, hemi, tname)
            if tnum < 0:
                altNames = troi[3]
                    
                cnt = 0
                while tnum < 0 and cnt < len(altNames):
                    tnum = get_roi_num(ctab, hemi, altNames[cnt])

                    if tnum > 0:
                        print("INFO: %s: Using alternative ROI name of %s: %s" \
                                  % (get_list_roi_nums.__name__, \
                                     tname, altNames[cnt]))
                        break
                        
                    cnt += 1

            if tnum < 0:
                raise Exception, "Cannot determine the ROI number for ROI %s in color table file %s" % (tname, ctab)

            roiNums.append(tnum)

    return (rois, roiNums)

def gen_parc_masks(rois, roiNums, parcVol, outDir, doVolStats=True, redo=False):
    import os
    import numpy as np
    from scipy.io import savemat as savemat
    from scai_utils import check_bin_path, check_file, \
                           cmd_stdout, saydo

    check_bin_path("fslstats")

    volStats = {"roiName": [], "nVoxels": [], "mm3": []}
    for (i0, roi) in enumerate(rois):
        roiNum = roiNums[i0]
            
        maskFN = os.path.join(outDir, "%s.diff.nii.gz" % roi)
        binCmd = "mri_binarize --i %s --match %d --o %s" % \
                 (parcVol, roiNum, maskFN)
            
        if not os.path.isfile(maskFN) or redo:
            saydo(binCmd)    
        check_file(maskFN)

        #= Volume stats =#
        (so, se) = cmd_stdout("fslstats %s -V" % maskFN)
        assert(len(se) == 0)
            
        so = so.split(" ")
        assert(len(so) >= 2)
        
        volStats["nVoxels"].append(int(so[0]))
        volStats["mm3"].append(float(so[0]))
        volStats["roiName"].append(roi)
            
    if doVolStats:
        volStats["roiName"] = np.array(volStats["roiName"])
        volStats["nVoxels"] = np.array(volStats["nVoxels"])
        volStats["mm3"] = np.array(volStats["mm3"])
    
        volStatsFN = os.path.join(outDir, "vol_stats.mat")
        savemat(volStatsFN, volStats)
        check_file(volStatsFN)

        print("INFO: %s: Saved volume stats of the mask files at\n\t%s" \
              % (gen_parc_masks.__name__, volStatsFN))
