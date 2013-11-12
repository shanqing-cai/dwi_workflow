#!/usr/bin/env python

import os
import sys
import argparse

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from scai_utils import *

from dwi_analysis_settings import DWI_ANALYSIS_DIR
from dwi_project_info import projInfo

COLORS = [[0, 0, 1], [0, 1, 0], [1, 0, 0], [0, 0.5, 1],
          [0.5, 1, 0], [0.5, 0, 1], [0.5, 0.5, 0], [0, 1, 1],
          [1, 0.25, 0.75], [0.25, 1, 0.75], [0.75, 0.25, 1], [1, 0.75, 0.25],
          [0.5, 0.5, 0.25], [0.5, 0.25, 0.75], [0.75, 0.25, 0.5], [0.75, 0.5, 0.5]]



# All analysis types available 
ANALYSIS_TYPES = ["FA", "MD"]

# All grouping schemes
GROUPING_SCHEMES = ["byStudy", "byGroup"]

def makeOptList(lst):
    str = "{"

    for (i0, t_item) in enumerate(lst):
        str += "%s, " % t_item

    str = str[:-2] + "}"

    return str

def get_tensor_measure(tmFile, meas, roi, wmDepth):
    #=== Load tensor-measure data from mat file ===#
    from scipy.io import loadmat
    tmdat = loadmat(tmfile)

    #== Locate the depth ==#
    bDepthFound = 0
    for jd in range(len(tmdat["wmDepths"])):
        if wmDepth == int(tmdat["wmDepths"][jd]):
            bDepthFound = 1
            break

    if bDepthFound == 0:
        return np.nan

    #== Locate roi ==#
    bROIFound = 0
    roiNames = list(tmdat["roiNames"])
    for (jr, t_roi) in enumerate(roiNames):
        if t_roi.strip() == roi:
            bROIFound = 1
            break

    if bROIFound == 0:
        return np.nan

    t_dat = tmdat["tensMeas"][meas][0][0][jr][jd]
    
    return t_dat

if __name__ == "__main__":
    ap = argparse.ArgumentParser("Group-level DWI analysis")
    ap.add_argument("anaType",
                    help="Analysis type: " + makeOptList(ANALYSIS_TYPES))
    ap.add_argument("grpScheme",
                    help="Grouping scheme: " + makeOptList(GROUPING_SCHEMES))
    ap.add_argument("--roi",
                    help="Region of interest (ROI) name (e.g., lh_vPMC)",
                    dest="roi")
    ap.add_argument("--wmDepth", type=int, 
                    help="White-matter depth (e.g., -1 for gm, 1, 2, 3)")
    ap.add_argument("--groups", type=str, dest="groups",
                    help="Groups to analyze " + \
                         "(compulsary for the grouping scheme byGroup; " + \
                         "optional for other groupig schcmes )")
    ap.add_argument("--exclude-subjects", type=str, dest="exclS", 
                    help="Exclude subjects, " + \
                         "e.g., --exclude-subjects CAT_12,SEQPDS_SEQ02P10")
    ap.add_argument("-v", help="Verbose mode",
                    dest="bVerbose", action="store_true")

    if len(sys.argv) == 1:
        ap.print_help()
        sys.exit(0)

    args = ap.parse_args()

    if ANALYSIS_TYPES.count(args.anaType) != 1:
        error_log("Unrecognized anlaysis type: %s" % args.anaType)

    if GROUPING_SCHEMES.count(args.grpScheme) != 1:
        error_log("Unrecognized grouping scheme: %s" % args.grpScheme)

    #=== Build list of subjects to be excluded ===#
    sExcludeList = []
    if args.exclS != None:
        sExcludeList = args.exclS.split(",")

    anaGroups = [] # Groups to be included in the analysis 
    if args.groups != None:
        anaGroups = args.groups.split(",")

    #=== Check conditionally required arguments ===#
    if args.anaType == "FA" or args.anaType == "MD":
        if args.roi == None:
            error_log("The option --roi is required for the anlaysis type %s" %
                      args.anaType)
        if args.wmDepth == None:
            error_log("The option --wmDepth is required for the analysis type %s" % args.anaType)

    if args.grpScheme == "byGroup" \
       and (args.groups == None or len(args.groups) == 0):
        error_log("--groups not specified under grouping scheme: %s" \
                  % args.grpScheme)

    #=== Load data ===#    
    dat = [] # Data
    grp = [] # Group info: it can hold things like study ID and group ID
    
    for (i0, t_proj) in enumerate(projInfo["name"]):
        ns = len(projInfo["subjIDs"][i0]) # Number of subjects
        nsa = 0 # Number of subjects with data available
        
        if args.bVerbose:
            info_log("Loading data from project %s" % t_proj)

        for (i1, t_sid) in enumerate(projInfo["subjIDs"][i0]):
            t_psid = "%s_%s" % (t_proj, t_sid)

            #== Exclude subject per input argument ==#
            if sExcludeList.count(t_psid) > 0:
                if args.bVerbose:
                    info_log("Excluding subject %s (request in options)" \
                             % t_psid)
                continue

            #== Check the group identity of the subject ==#
            t_grp_id = projInfo["groupIDs"][i0][i1]
            if len(anaGroups) > 0 and anaGroups.count(t_grp_id) == 0:
                if args.bVerbose:
                    info_log("Excluding subject %s (outside groups to be analyzed)" % t_psid)
                    continue

            sDir = os.path.join(DWI_ANALYSIS_DIR, t_psid)

            #=== Load data ===#
            if args.anaType == "FA" or args.anaType == "MD":
                # Tensor measure file
                tmfile = os.path.join(sDir, "tensor", "tensor_measures.mat")

                if not os.path.isfile(tmfile):
                    continue

                t_dat = get_tensor_measure(tmfile, args.anaType, args.roi,
                                           args.wmDepth)
                if np.isnan(t_dat):
                    continue
                
                dat.append(t_dat)
               
            else:
                error_log("Analaysis type %s has not been implemented yet" \
                          % args.anaType)

            #=== Add group identity info ===#
            if args.grpScheme == "byStudy":
                grp.append(i0)
            elif args.grpScheme == "byGroup":
                grp.append(anaGroups.index(t_grp_id))
            else:
                error_log("Unrecognized grouping scheme: %s" \
                          % args.grpScheme)
            
            nsa += 1

        if args.bVerbose:
            info_log("\tData loaded from %s of %s subjects" % (nsa, ns))

    #=== Clean-up data ===#
    dat = np.array(dat)
    grp = np.array(grp)


    #=== Statistical analysis and visualizaton ===#
    import scipy.stats as stats
    if args.grpScheme == "byStudy":
#        plt.figure()

        ugrps = np.unique(grp)

        plt.figure(facecolor=[1.0, 1.0, 1.0])
        datArray = []
        for (i0, t_grp) in enumerate(ugrps):
            idx = [i for (i, x) in enumerate(grp) if x == t_grp]
            
            if len(idx) > 0:
                datArray.append(dat[idx])
                plt.plot(idx, dat[idx], "o-", color=COLORS[i0])
                plt.text(idx[0] + 0.5, np.min(dat[idx]) * 0.95,
                         projInfo["name"][t_grp],
                         color=COLORS[i0])

        plt.xlabel("Subject #")
        plt.ylabel(args.anaType + ": " + args.roi)
        
        (ANOVA_F, ANOVA_p) = stats.f_oneway(*datArray)

        xs = plt.xlim()
        ys = plt.ylim()

        plt.text(xs[0] + 0.05 * (xs[1] - xs[0]),
                 ys[1] - 0.075 * (ys[1] - ys[0]),
                 "One-way ANOVA: F = %f; p = %e" % (ANOVA_F, ANOVA_p))

        plt.show()

    elif args.grpScheme == "byGroup":
        plt.figure(facecolor=[1.0, 1.0, 1.0])

        ugrps = np.unique(grp)

        for (i0, t_grp) in enumerate(anaGroups):
            idx = [i for (i, x) in enumerate(grp) if x == i0]

            if len(idx) > 0:
                plt.plot(i0 + 1, np.mean(dat[idx]), "o", color=COLORS[i0])
                plt.plot([i0 + 1, i0 + 1], 
                         np.mean(dat[idx]) + \
                            np.array([-1, 1]) * np.std(dat[idx]) \
                               / np.sqrt(len(idx)),
                         "-", color=COLORS[i0])
        plt.xlim([0, len(anaGroups) + 1])
        plt.xticks(np.linspace(1, len(anaGroups), len(anaGroups)))
        plt.xlabel("Groups ")
        # plt.set_xticklabels(anaGroups)
#                plt.text(idx[0], dat[idx[0]] * 1.2, anaGroups[i0],
#                         color=COLORS[i0])

        plt.show()
        
    else:
        error_log("Grouping scheme %s has not been implemented yet" \
                  % args.grpScheme)
            

        
            
