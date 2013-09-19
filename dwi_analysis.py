#!/speechlab/software/EPD7/epd-7.1-1-rh5-x86_64/bin/python

import os
import sys
import glob

import argparse
import tempfile

from scai_utils import *
from mri_utils import *

from dwi_project_info import projInfo
from dwi_analysis_settings import DWI_ANALYSIS_DIR, \
    BASE_TRACULA_CFG, FIX_BVECS_SCRIPT, \
    TRACULA_DOEDDY, TRACULA_DOROTVECS, TRACULA_THR_BET, \
    INIT_FLIRT_MATS

from scan_for_dwi import extract_dwi_from_dicom

ALL_STEPS = ["convert", "dtiprep", "postqc", \
             "tracula_prep", "tracula_bedp", \
             "fix_coreg", \
             "inspect_tensor", "inspect_coreg", \
             "parcellate", "inspect_parc", \
             "probtrackx", "cort_conn_mat"]

HEMIS = ["lh", "rh"]

if __name__ == "__main__":
    stepsHelp = ""
    for str in ALL_STEPS:
        stepsHelp += "%s | " % str
    stepsHelp = stepsHelp[:-3]
    
    ap = argparse.ArgumentParser(description="DWI workflow")
    ap.add_argument("projName", help="Project name (e.g., CAT, SEQPDS, STUT)")
    ap.add_argument("subjID", help="Subject ID")
    #ap.add_argument("--fsSubjID", dest="fsSubjID", default="", \
    #                required=False, \
    #                help="FreeSurfer subject ID (required by the steps: tracula_prep and inspect_coreg)")
    ap.add_argument("step", help="Steps: {%s}" % stepsHelp)
    ap.add_argument("--parc", dest="parcName", default="", 
                    help="Cortical parcellation name (required by step parcellate; e.g., aparc12)")
    ap.add_argument("--redo", dest="bRedo", action="store_true", \
                    help="Force redoing time-consuming steps (default=False; not fully implemented yet)")
    ap.add_argument("--seed", dest="seed", nargs=2, 
                    help="Seed ROI for probtrackx: seed_roi, seed_type (e.g., lh_vPMC gm)")
    ap.add_argument("--targ", dest="targ", nargs=2, 
                    help="Target ROI for probtrackx: seed_roi, seed_type (e.g., lh_vSC gm. Use lh_all or rh_all for all ROIs in one hemisphere of the parcellation paradigm)")
    ap.add_argument("--hemi", dest="hemi", 
                    help="Hemisphere (required by the cort_conn_mat step)")
    ap.add_argument("--maskType", dest="maskType", 
                    help="Type of mask for tractography (e.g, gm, wm1mm), required by the cort_conn_mat step")
    ap.add_argument("--speech", dest="bSpeech", action="store_true", 
                    help="Use speech subnetwork (option for step cort_conn_mat")
    ap.add_argument("--skip-mris-ca-label", dest="bSkipMRISCALabel",
                    action="store_true",
                    help="Skip the mris_ca_label sub-step in step parcellate")

    #                help="Seed ROI for probtrackx")
    #ap.add_argument("--targ", dest="targ", type=str, \
    #                help="Target ROI for probtrackx (required by step probtrackx, e.g., lh_vSC")

    if len(sys.argv) == 1:
        ap.print_help()
        sys.exit(0)

    args = ap.parse_args()

    #=== Locate entry in projInfo ===#
    pidx = projInfo["name"].index(args.projName)
    sidx = projInfo["subjIDs"][pidx].index(args.subjID)

    rawFormat = projInfo["rawFormat"][pidx].upper()
    rawFN = projInfo["rawPath"][pidx]
    
    FS_SUBJECTS_DIR = projInfo["fsSubjectsDir"][pidx]
    check_dir(FS_SUBJECTS_DIR)

    #== Locate the input 4d data ==#
    rawFN = rawFN.replace("{subjID}", "%s")
    rawFN %= args.subjID
    rawFN = glob.glob(rawFN)

    assert(len(rawFN) == 1)
    rawFN = rawFN[0]

    if rawFormat == "DICOM":
        check_dir(rawFN)
    elif rawFormat == "NRRD":
        check_file(rawFN)
    elif rawFormat == "NGZ":
        check_file(rawFN)
    else:
        error_log("Unrecognized raw input format: %s" % rawFormat, \
                  logFN=logFileName)

    #== Create directory ==#
    check_dir(DWI_ANALYSIS_DIR)
    sID = "%s_%s" % (args.projName, args.subjID)
    sDir = os.path.join(DWI_ANALYSIS_DIR, sID)
    check_dir(sDir, bCreate=True)

    #=== Prepare log file name ===#
    if ALL_STEPS.count(args.step) == 0:
        error_log("Unrecognized step name: %s" % args.step, logFN=logFileName)
    
    logDir = os.path.join(sDir, "logs")
    check_dir(logDir, bCreate=True)
    
    logFileName = os.path.join(logDir, args.step + ".log")
    info_log("logFileName = %s" % logFileName)

    #== Prepare directories and some file names ==#
    dicomDir = os.path.join(sDir, "dicom")
    check_dir(dicomDir, bCreate=True, logFN=logFileName)
    
    dtiprepDir = os.path.join(sDir, "dtiprep")
    check_dir(dtiprepDir, bCreate=True, logFN=logFileName)
        
    outputNrrd = os.path.join(dtiprepDir, "raw_dwi.nrrd")

    qcedNrrd = os.path.join(dtiprepDir, "raw_dwi_QCed.nrrd")
    qcXML = os.path.join(dtiprepDir, "raw_dwi_XMLQCResult.xml")
    qcReport = os.path.join(dtiprepDir, "raw_dwi_QCReport.txt")

    qcedNGZ = os.path.join(dtiprepDir, "dwi_qced.nii.gz")
    qcedBVals = os.path.join(dtiprepDir, "dwi_qced.bvals")
    qcedBVecs = os.path.join(dtiprepDir, "dwi_qced.bvecs")

    dmriDir = os.path.join(sDir, "dmri")
    faFN = os.path.join(dmriDir, "dtifit_FA.nii.gz")
    v1FN = os.path.join(dmriDir, "dtifit_V1.nii.gz")
    lowbBrainFN = os.path.join(dmriDir, "lowb_brain.nii.gz")
    brainAnatOrigFN = os.path.join(dmriDir, "brain_anat_orig.nii.gz")

    xfmsDir = os.path.join(dmriDir, "xfms")
    diff2anatorig_mat = os.path.join(xfmsDir, "diff2anatorig.bbr.mat")
    anatorig2diff_mat = os.path.join(xfmsDir, "anatorig2diff.bbr.mat")

    traculaCfgFN = os.path.join(sDir, "tracula.cfg")

    annotDir = os.path.join(sDir, "annot")

    tracksDir = os.path.join(sDir, "tracks")

    connDir = os.path.join(sDir, "conn")

    #== Determine the FreeSurfer subject ID ==#
    assert(len(projInfo["subjIDs"][pidx]) == len(projInfo["fsSubjIDs"][pidx]))
    fsSubjID = projInfo["fsSubjIDs"][pidx][sidx]

    traculaDir = os.path.join(DWI_ANALYSIS_DIR, fsSubjID)
    tracula_dlabelDir = os.path.join(traculaDir, "dlabel")
    tracula_dmriDir = os.path.join(traculaDir, "dmri")
    tracula_scriptsDir = os.path.join(traculaDir, "scripts")

    tracula_bedpDir = os.path.join(traculaDir, "dmri.bedpostX")
    bedpDir = os.path.join(sDir, "dmri.bedpostX")

    #== Locate the bvals and bvecs files ==#
    if len(projInfo["bvalsPath"][pidx]) > 0:
        bvalsFN = \
            projInfo["bvalsPath"][pidx].replace("{subjID}", "%s") % args.subjID
        bvalsFN = glob.glob(bvalsFN)
        assert(len(bvalsFN) == 1)
        bvalsFN = bvalsFN[0]
        check_file(bvalsFN, logFN=logFileName)

        bvecsFN = \
            projInfo["bvecsPath"][pidx].replace("{subjID}", "%s") % args.subjID
        bvecsFN = glob.glob(bvecsFN)
        assert(len(bvecsFN) == 1)
        bvecsFN = bvecsFN[0]
        check_file(bvecsFN, logFN=logFileName)
    else:
        #== DICOM files without accompanying bvals and bvecs files ==#
        if args.step == "convert":
            out = extract_dwi_from_dicom(rawFN, dicomDir)

            assert(len(out["dwi4d"]) == len(out["bvals"]))
            assert(len(out["dwi4d"]) == len(out["bvecs"]))

            if len(out["dwi4d"]) == 0:
                error_log("Cannot find any DWI series in input DICOM directory: %s" % rawFN,
                          logFN=logFileName)
            elif len(out["dwi4d"]) > 1:
                print("WARNING: more than one diffusion series found in DICOM directory: %s. Using the last one.")
                
            rawFormat = "NGZ"
            rawFN = out["dwi4d"][-1]
            bvalsFN = out["bvals"][-1]
            bvecsFN = out["bvecs"][-1]
        else:
            rawFormat = "NGZ"
            rawFN = glob.glob(os.path.join(dicomDir, "*.nii.gz"))[0]
            bvalsFN = glob.glob(os.path.join(dicomDir, "*.bval*"))[0]
            bvecsFN = glob.glob(os.path.join(dicomDir, "*.bvec*"))[0]

        check_file(rawFN, logFN=logFileName)
        check_file(bvalsFN, logFN=logFileName)
        check_file(bvecsFN, logFN=logFileName)


    #=== Determine the rotMat (for postqc) ===#
    assert(len(projInfo["subjIDs"][pidx]) == len(projInfo["rotMat"][pidx]))
    rotMat = projInfo["rotMat"][pidx][sidx]

    #=== Check SUBJECTS_DIR ===#
    fsSubjectsDir = os.getenv("SUBJECTS_DIR")
    check_dir(fsSubjectsDir, logFN=logFileName)

    if fsSubjectsDir != FS_SUBJECTS_DIR:
        error_log("Your environmental SUBJECTS_DIR (%s) does not match the variable FS_SUBJECTS_DIR (%s) in dwi_analysis_settings.py" % (fsSubjectsDir, FS_SUBJECTS_DIR),
                  logFN=logFileName)

    #==== Main branches ====#
    info_log("Step %s started" % args.step, logFN=logFileName)
    
    if args.step == "convert":
        #=== Run DWIConvert ===#
        #== Check the path to DWIConvert ==#
        check_bin_path("DWIConvert", logFN=logFileName)
        
        from dtiprep_utils import format_bvals_bvecs
        format_bvals_bvecs(bvalsFN, bvecsFN)

        if rawFormat == "DICOM":
            #== Copy the dicoms ==#    
            saydo("cp -v %s %s" % (os.path.join(rawFN, "*.dcm"), dicomDir),
                  logFN=logFileName)

            #== Copy the bvecs and bvals ==#
            saydo("cp -v %s %s" % (bvalsFN, dicomDir), logFN=logFileName)
            bvalsFN = os.path.join(dicomDir, os.path.split(bvalsFN)[-1])

            saydo("cp -v %s %s" % (bvecsFN, dicomDir), logFN=logFileName)
            bvecsFN = os.path.join(dicomDir, os.path.split(bvecsFN)[-1])

            dwiCvtCmd = "DWIConvert --inputDicomDirectory %s " % dicomDir + \
                        "--inputBValues %s " % bvalsFN + \
                        "--inputBVectors %s " % bvecsFN + \
                        "--conversionMode DicomToNrrd " + \
                        "--outputVolume %s " % outputNrrd

        elif rawFormat == "NGZ":
            dwiCvtCmd = "DWIConvert --inputVolume %s " % rawFN + \
                        "--inputBValues %s " % bvalsFN + \
                        "--inputBVectors %s " % bvecsFN + \
                        "--conversionMode FSLToNrrd " + \
                        "--outputVolume %s " % outputNrrd
            
        else:
            error_log("Step %s: Unsupported input format %s" % \
                      (args.step, rawFormat),
                      logFN=logFileName)

        saydo(dwiCvtCmd, logFN=logFileName)
        check_file(outputNrrd, logFN=logFileName)

        info_log("outputNrrd = %s" % outputNrrd, logFN=logFileName)

        # 
            
    elif args.step == "dtiprep":
        #=== Run dtiprep ===#
        check_file(outputNrrd, logFN=logFileName)

        #== Check the path to DTIPrep ==#
        (so, se) = cmd_stdout("which DTIPrep")
        if len(se) > 0 or len(so) == 0:
            error_log("Cannot find the path to executable: DTIPrep",
                      logFN=logFileName)

        logFN = os.path.join(dtiprepDir, "dtiprep_notes.txt")
        dtiprepCmd = "DTIPrep -w %s -p default -d -c -n %s" % \
                     (outputNrrd, logFN)

        saydo(dtiprepCmd, logFN=logFileName)
        check_file(qcedNrrd, logFN=logFileName)
        check_file(qcXML, logFN=logFileName)
        check_file(qcReport, logFN=logFileName)

    elif args.step == "postqc":
        #=== Convert the qc'ed nrrd to nifti ===#
        check_bin_path("DWIConvert", logFN=logFileName)

        qcedBVecs_unc = os.path.join(dtiprepDir, "dwi_qced.bvecs.uncorrected")

        check_file(qcedNrrd, logFN=logFileName)
        cvtCmd = "DWIConvert --inputVolume %s " % qcedNrrd + \
                 "--conversionMode NrrdToFSL " + \
                 "--outputBValues %s " % qcedBVals + \
                 "--outputBVectors %s " % qcedBVecs_unc + \
                 "--outputVolume %s " % qcedNGZ

        saydo(cvtCmd, logFN=logFileName)
        check_file(qcedNGZ, logFN=logFileName)
        check_file(qcedBVals, logFN=logFileName)
        check_file(qcedBVecs_unc, logFN=logFileName)

        #=== Correct the b-vectors for FSL and FreeSurfer ===#
        from dtiprep_utils import correctbvec4fsl
        correctbvec4fsl(qcedNGZ, qcedBVecs_unc, qcedBVecs, rotMat)
        check_file(qcedBVecs, logFN=logFileName)
        

    elif args.step == "tracula_prep":
        #=== Tracula: step prep ===#
        import numpy as np
        from tracula_utils import modify_tracula_cfg_file

        #=== Check the version of FreeSurfer ===#
        from freesurfer_utils import check_fs_ver
        check_fs_ver(5.1, mode="eqgt")
        
        check_bin_path("trac-all", logFN=logFileName)

        #=== Check input FreeSurfer subject ID ===#
        if len(fsSubjID) == 0:
            error_log("Input argument --fsSubjID must be set for step %s" % args.step,
                      logFN=logFileName)

        check_file(qcedNGZ, logFN=logFileName)
        check_file(qcedBVals, logFN=logFileName)
        check_file(qcedBVecs, logFN=logFileName)

        check_file(BASE_TRACULA_CFG, logFN=logFileName)

        fsSDir = os.path.join(fsSubjectsDir, fsSubjID)
        check_dir(fsSDir, logFN=logFileName)

        #=== Determine nb0 ===#
        bv = np.genfromtxt(qcedBVals)
        nb0 = len(np.nonzero(bv == 0.0))

        fsHome = os.getenv("FREESURFER_HOME")
        modify_tracula_cfg_file(fsHome, fsSubjectsDir, \
                                BASE_TRACULA_CFG, DWI_ANALYSIS_DIR, \
                                fsSubjID, \
                                dtiprepDir, os.path.split(qcedNGZ)[-1], \
                                qcedBVals, qcedBVecs, nb0, \
                                TRACULA_DOEDDY, TRACULA_DOROTVECS, \
                                TRACULA_THR_BET, traculaCfgFN)
        check_file(traculaCfgFN, logFN=logFileName)
        print("INFO: TRACULA config file for subject %s written to: %s" % \
              (sID, traculaCfgFN))

        #=== Execute tracula ===#
        tracall_cmd = "trac-all -c %s -prep" % traculaCfgFN

        saydo(tracall_cmd, logFN=logFileName)
        
        check_dir(traculaDir, logFN=logFileName)
        check_dir(tracula_dlabelDir, logFN=logFileName)
        check_dir(tracula_dmriDir, logFN=logFileName)
        check_dir(tracula_scriptsDir, logFN=logFileName)

        #=== Move the tracula files back to the dwi directory, if necessary ===#
        if sID != fsSubjID:
            saydo("mv %s %s/" % (tracula_dlabelDir, sDir), logFN=logFileName)
            saydo("mv %s %s/" % (tracula_dmriDir, sDir), logFN=logFileName)
            saydo("mv %s %s/" % (tracula_scriptsDir, sDir), logFN=logFileName)
            saydo("rmdir %s" % traculaDir, logFN=logFileName)


    elif args.step == "tracula_bedp":
        from freesurfer_utils import check_fs_ver
        check_fs_ver(5.1, mode="eqgt")

        check_bin_path("trac-all", logFN=logFileName)
        
        #== Check some prerequisite files ==#
        check_file(traculaCfgFN, logFN=logFileName)
        check_file(lowbBrainFN, logFN=logFileName)
        check_file(diff2anatorig_mat, logFN=logFileName)

        s_dlabelDir = os.path.join(sDir, "dlabel")
        s_dmriDir = os.path.join(sDir, "dmri")
        s_scriptsDir = os.path.join(sDir, "scripts")

        check_dir(s_dlabelDir, logFN=logFileName)
        check_dir(s_dmriDir, logFN=logFileName)
        check_dir(s_scriptsDir, logFN=logFileName)

        #== Copy files into tracula directory if the file subject ID #
        #   differs between DWI and FreeSurfer ==#
        if sID != fsSubjID:
            if os.path.isdir(traculaDir):
                error_log("Tracula directory %s already exists. Remove it before proceeding" % traculaDir,
                          logFN=logFileName)
            
            saydo("mkdir %s" % traculaDir, logFN=logFileName)
            saydo("cp -r %s %s" % (s_dlabelDir, tracula_dlabelDir),
                  logFN=logFileName)
            saydo("cp -r %s %s" % (s_dmriDir, tracula_dmriDir),
                  logFN=logFileName)
            saydo("cp -r %s %s" % (s_scriptsDir, tracula_scriptsDir),
                  logFN=logFileName)

        tracall_cmd = "trac-all -c %s -bedp" % traculaCfgFN
        saydo(tracall_cmd, logFN=logFileName)

        check_dir(tracula_bedpDir, logFN=logFileName)

        from tracula_utils import check_bedp_complete
        if not check_bedp_complete(tracula_bedpDir):
            error_log("It appears that bedpostx failed.", logFN=logFileName)

        if sID != fsSubjID:
            saydo("cp -r %s %s/" % (tracula_dlabelDir, sDir), logFN=logFileName)
            saydo("cp -r %s %s/" % (tracula_dmriDir, sDir), logFN=logFileName)
            saydo("cp -r %s %s/" % (tracula_scriptsDir, sDir),
                  logFN=logFileName)
            saydo("cp -r %s %s/" % (tracula_bedpDir, sDir), logFN=logFileName)
            saydo("rm -rf %s" % traculaDir, logFN=logFileName)

        check_dir(bedpDir, logFN=logFileName)
        
        if check_bedp_complete(bedpDir):
            info_log("bedpostx completed successfully", logFN=logFileName)
        else:
            error_log("It appears that bedpostx failed.", logFN=logFileName)

    elif args.step == "inspect_tensor":
        check_dir(dmriDir, logFN=logFileName)
        check_file(faFN, logFN=logFileName)
        check_file(v1FN, logFN=logFileName)

        #=== Check fslview ===#
        check_bin_path("fslview", logFN=logFileName)
        
        viewCmd = "fslview %s %s" % (faFN, v1FN)
        saydo(viewCmd, logFN=logFileName)

    elif args.step == "inspect_coreg":
        check_dir(xfmsDir, logFN=logFileName)

        #=== tkregister2 path check ===#
        check_bin_path("tkregister2", logFN=logFileName)

        #=== Check input FreeSurfer subject ID ===#
        if len(fsSubjID) == 0:
            error_log("Input argument --fsSubjID must be set for step %s" % args.step,
                      logFN=logFileName)
        
        #=== Check the diffusion to anatomical coregistration ==#
        check_file(lowbBrainFN, logFN=logFileName)

        #=== Locate the FreeSurfer T1 file ===#
        fsT1FN = os.path.join(FS_SUBJECTS_DIR, fsSubjID, "mri", "T1.mgz")
        check_file(fsT1FN, logFN=logFileName)

        fsT1NGZ = os.path.join(FS_SUBJECTS_DIR, fsSubjID, \
                               "mri", "T1.nii.gz")
        cvtCmd = "mri_convert %s %s" % (fsT1FN, fsT1NGZ)
        saydo(cvtCmd, logFN=logFileName)
        check_file(fsT1NGZ, logFN=logFileName)

        d2ao = os.path.join(xfmsDir, "diff2anatorig.bbr.mat")
        faAnatOrig = os.path.join(dmriDir, "dtifit_FA_anatorig.nii.gz")
        flirt_apply_xfm(faFN, fsT1NGZ, d2ao, faAnatOrig)
        check_file(faAnatOrig, logFN=logFileName)

        tmpIdentity = tempfile.mktemp() + ".dat"
        tkrCmd = "tkregister2 --targ %s --mov %s --identity --reg %s " \
                 % (fsT1NGZ, faAnatOrig, tmpIdentity) + \
                 "--s %s --surfs " % (fsSubjID)
        saydo(tkrCmd, logFN=logFileName)
        saydo("rm -f %s" % tmpIdentity, logFN=logFileName)

    elif args.step == "fix_coreg":
        #=== Files and directories check ===#
        check_dir(xfmsDir, logFN=logFileName)
        check_file(diff2anatorig_mat, logFN=logFileName)

        check_file(lowbBrainFN, logFN=logFileName)
        check_file(brainAnatOrigFN, logFN=logFileName)

        assert(len(INIT_FLIRT_MATS) > 0)
        initMat = INIT_FLIRT_MATS[0]
        check_file(initMat, logFN=logFileName)

        #=== Programs: path checks ===#
        check_bin_path("tkregister2", logFN=logFileName)
        check_bin_path("flirt", logFN=logFileName)

        #=== flirt for initialization ===#
        flirtInitOut = os.path.join(dmriDir, "lowb_brain_flirt2anatorig.nii.gz")
        flirtInitMat = os.path.join(xfmsDir, "flirt_diff2anatorig.mat")
        flirtCmd = "flirt -in %s -ref %s -out %s -omat %s -init %s" \
                   % (lowbBrainFN, brainAnatOrigFN, flirtInitOut, \
                      flirtInitMat, initMat)

        saydo(flirtCmd, logFN=logFileName)
        check_file(flirtInitOut, logFN=logFileName)
        check_file(flirtInitMat, logFN=logFileName)

        #=== Use tkregister2 to convert the init mat to dat ===#
        flirtInitDat = os.path.join(xfmsDir, "flirt_diff2anatorig.dat")
        cvtCmd = "tkregister2 --mov %s --fsl %s --reg %s --s %s --noedit" \
            % (lowbBrainFN, flirtInitMat, flirtInitDat, fsSubjID)
        saydo(cvtCmd, logFN=logFileName)
        check_file(flirtInitDat, logFN=logFileName)

        #=== Backup the old xfm ===#
        diff2anatorig_mat_backup = os.path.join(xfmsDir, \
                                                "diff2anatorig.bbr.mat.old")
        saydo("mv %s %s" % (diff2anatorig_mat, diff2anatorig_mat_backup), logFN=logFileName)
        check_file(diff2anatorig_mat_backup, logFN=logFileName)

        #=== Run bbregister2 ===#
        diff2anatorig_dat = os.path.join(xfmsDir, "diff2anatorig.bbr.dat")
        bbrCmd = "bbregister --s %s --init-reg %s --dti --mov %s " \
                 % (fsSubjID, flirtInitDat, lowbBrainFN) + \
                 "--reg %s --fslmat %s " \
                 % (diff2anatorig_dat, diff2anatorig_mat)

        saydo(bbrCmd, logFN=logFileName)
        check_file(diff2anatorig_dat, logFN=logFileName)
        check_file(diff2anatorig_mat, logFN=logFileName)

    elif args.step == "parcellate":
        #=== Cortical parcellation using FreeSurfer ===#
        if len(args.parcName) == 0:
            error_log("parcellation name (--parc) must be supplied for step %s" % args.step,
                      logFN=logFileName)

        from dwi_analysis_settings \
            import SURF_CLASSIFIERS, PARC_FS_VER

        #== Check FreeSurfer version ==#
        from freesurfer_utils import check_fs_ver
        check_fs_ver(PARC_FS_VER, mode="eq")

        parcIdx = SURF_CLASSIFIERS["name"].index(args.parcName)
        list_py = SURF_CLASSIFIERS["list_py"][parcIdx]
        gcsw = SURF_CLASSIFIERS["gcs"][parcIdx]
        ctab = SURF_CLASSIFIERS["ctab"][parcIdx]

        if not args.bSkipMRISCALabel:
            check_bin_path("mris_ca_label", logFN=logFileName)
        
        check_file(ctab, logFN=logFileName)

        check_dir(annotDir, bCreate=True, logFN=logFileName)

        labelDir = os.path.join(FS_SUBJECTS_DIR, fsSubjID, "label")

        gcss = {}
        annotFNs = {}
        #== Step 1: mris_ca_label ==#
        for hemi in HEMIS:
            annotFNs[hemi] = os.path.join(labelDir, \
                                          "%s.%s.annot" \
                                          % (hemi, args.parcName))
            
            if not args.bSkipMRISCALabel:
                gcss[hemi] = gcsw.replace("{hemi}", hemi)
                check_file(gcss[hemi], logFN=logFileName)
                
                labelCmd = "mris_ca_label -t %s %s %s " % (ctab, fsSubjID, hemi) + \
                           "sphere.reg %s %s" % (gcss[hemi], annotFNs[hemi])

                #if not os.path.isfile(annotFNs[hemi]) or args.bRedo:
                saydo(labelCmd, logFN=logFileName)

            check_file(annotFNs[hemi], logFN=logFileName)

        #== Step 2: mri_aparc2aseg ==#
        from dwi_analysis_settings import WM_DEPTHS
        wmDepths = WM_DEPTHS

        check_bin_path("mri_aparc2aseg", logFN=logFileName)

        from mri_utils import invert_fsl_xfm_mat
        from mri_utils import flirt_apply_xfm
        invert_fsl_xfm_mat(diff2anatorig_mat, anatorig2diff_mat)
        
        for i0 in range(len(wmDepths) + 1):
            if i0 == len(wmDepths):                
                parcVol = os.path.join(annotDir, "%s.nii.gz" % args.parcName)
            else:
                depth = wmDepths[i0]
                assert(type(depth) == int)
                assert(depth > 0)
                parcVol = os.path.join(annotDir, \
                                       "%s.%smm.nii.gz" \
                                       % (args.parcName, depth))

            segCmd = "mri_aparc2aseg --s %s --o %s " % (fsSubjID, parcVol) + \
                     "--annot %s --labelwm " % (args.parcName) + \
                     "--annot-table %s " % (ctab)
            
            if i0 != len(wmDepths):
                segCmd += "--wmparc-dmax %d " % depth
        
            if not os.path.isfile(parcVol) or args.bRedo:
                saydo(segCmd, logFN=logFileName)
            
            check_file(parcVol, logFN=logFileName)
        
            #== Step 3: generate diffusion-space version of the parc volume ==#
            check_file(diff2anatorig_mat, logFN=logFileName)
            check_file(lowbBrainFN, logFN=logFileName)

            if i0 == len(wmDepths):
                parcVolDiff = os.path.join(annotDir, \
                                           "%s.diff.nii.gz" % args.parcName)
            else:
                parcVolDiff = os.path.join(annotDir, \
                                           "%s.%dmm.diff.nii.gz" \
                                           % (args.parcName, depth))

            if not os.path.isfile(parcVolDiff) or args.bRedo:
                flirt_apply_xfm(parcVol, lowbBrainFN, \
                                anatorig2diff_mat, parcVolDiff, \
                                interpMeth="nearestneighbour")
            check_file(parcVolDiff, logFN=logFileName)

        #== Step 4: Generate the individual ROI masks in the diffusion space ==#
        #==         This includes GM and WM ==#
        roiList = __import__(list_py)
        roiList = roiList.aROIs
        
        from aparc_utils import get_list_roi_nums
        (rois, roiNums) = get_list_roi_nums(roiList, HEMIS, ctab)

        import numpy as np
        roiNums = np.array(roiNums)

        parcMaskDir = os.path.join(annotDir, args.parcName)
        check_dir(parcMaskDir, bCreate=True, logFN=logFileName)

        #== Step 4.1: Gray matter (GM) ==#
        gmDir = os.path.join(parcMaskDir, "gm")
        check_dir(gmDir, bCreate=True, logFN=logFileName)

        from aparc_utils import gen_parc_masks
        parcVolDiff = os.path.join(annotDir, \
                                   "%s.diff.nii.gz" % args.parcName)
        gen_parc_masks(rois, roiNums, parcVolDiff, gmDir, \
                       doVolStats=True, redo=args.bRedo, logFN=logFileName)
        
        #== Step 4.2: White matter (WM) ==#
        wm_roiNums = roiNums + 2000;
        for (i0, depth) in enumerate(wmDepths):
            assert(type(depth) == int)
            wmDir = os.path.join(parcMaskDir, "wm%dmm" % depth)
            check_dir(wmDir, bCreate=True, logFN=logFileName)
            
            parcVolDiff = os.path.join(annotDir, \
                                       "%s.%dmm.diff.nii.gz" \
                                       % (args.parcName, depth))            
            gen_parc_masks(rois, wm_roiNums, parcVolDiff, wmDir, \
                           doVolStats=True, redo=args.bRedo, logFN=logFileName)
            
    elif args.step == "inspect_parc":
        if len(args.parcName) == 0:
            error_log("parcellation name (--parc) must be supplied for step %s" % args.step,
                      logFN=logFileName)

        parcVolDiff = os.path.join(annotDir, \
                                  "%s.diff.nii.gz" \
                                  % (args.parcName))
        check_file(parcVolDiff, logFN=logFileName)
        check_file(faFN, logFN=logFileName)

        tmpReg = tempfile.mktemp() + ".dat"

        tkrCmd = "tkregister2 --targ %s --mov %s --identity --reg %s " \
                 % (faFN, parcVolDiff, tmpReg)
        saydo(tkrCmd, logFN=logFileName)

        saydo("rm -rf %s " % (tmpReg), logFN=logFileName)
        
    elif args.step == "probtrackx":
        #=== Probabilistic tractography ===#
        from dwi_analysis_settings import SURF_CLASSIFIERS
        parcIdx = SURF_CLASSIFIERS["name"].index(args.parcName)
        list_py = SURF_CLASSIFIERS["list_py"][parcIdx]
        roiList = __import__(list_py)
        roiList = roiList.aROIs

        if args.parcName == "" or args.parcName == None:
            error_log("Compulsory --parc option not supplied for step %s" \
                      % args.step,
                      logFN=logFileName)
        
        if SURF_CLASSIFIERS["name"].count(args.parcName) == 0:
            error_log("Parcellation '%s' is not found in dwi_analysis_settings.py" \
                      % args.parcName,
                      logFN=logFileName)

        if args.seed == None:
            error_log("Compulsory --seed option not supplied for step %s" \
                      % args.step,
                      logFN=logFileName)
        
        seedList = []
        if args.seed[0].endswith("_all"):
            # All ROIs in one hemisphere (of the parcellation paradigm)
            t_hemi = args.seed[0].replace("_all", "")
            if HEMIS.count(t_hemi) == 0:
                error_log("Unrecognized hemisphere name: %s" % t_hemi,
                          logFN=logFileName)
            for (i0, troi) in enumerate(roiList):
                seedList.append("%s_%s" % (t_hemi, troi[0]))
        else:
            # Single seed ROI
            seedList.append(args.seed[0])
        seedList.sort()
        
        #== Target files ==#
        parcDir = os.path.join(annotDir, args.parcName)
        if args.targ != None:
            typeDir = os.path.join(parcDir, args.targ[1])
            check_dir(typeDir, logFN=logFileName)
            volStatsFN = os.path.join(typeDir, "vol_stats.mat")
            check_file(volStatsFN, logFN=logFileName)

            targMask = os.path.join(typeDir, "%s.diff.nii.gz" % args.targ[0])
            check_file(targMask, logFN=logFileName)
        else:
            targMask = None
        
        #== Check directories and files ==#
        check_dir(annotDir, logFN=logFileName)
        
        check_dir(parcDir, logFN=logFileName)

        typeDir = os.path.join(parcDir, args.seed[1])
        check_dir(typeDir, logFN=logFileName)
        volStatsFN = os.path.join(typeDir, "vol_stats.mat")
        check_file(volStatsFN, logFN=logFileName)
        
        check_dir(tracksDir, bCreate=True, logFN=logFileName)

        parcTracksDir = os.path.join(tracksDir, args.parcName)
        check_dir(parcTracksDir, bCreate=True, logFN=logFileName)

        #== bedp merged files and brain mask ==#
        check_dir(bedpDir, logFN=logFileName)

        from tracula_utils import check_bedp_complete
        check_bedp_complete(bedpDir)

        bedpBase = os.path.join(bedpDir, "merged")

        brainMask = os.path.join(bedpDir, "nodif_brain_mask.nii.gz")
        check_file(brainMask, logFN=logFileName)

        #== Iterate through all seed ROIs ==#
        from tractography import run_probtrackx

        for (i0, seedROI) in enumerate(seedList):
            #= Seed files =#
            seedMask = os.path.join(typeDir, "%s.diff.nii.gz" % seedROI)
            check_file(seedMask, logFN=logFileName)
            
            #== Prepare output directory ==#
            if args.targ != None:
                outDir = os.path.join(parcTracksDir, 
                              "%s_%s_to_%s_%s" % (seedROI, args.seed[1], 
                                                  args.targ[0], args.targ[1]))
            else:
                outDir = os.path.join(parcTracksDir, 
                                      "%s_%s" % (seedROI, args.seed[1]))

            check_dir(outDir, bCreate=True, logFN=logFileName)

            #== Do the work ==#
        
            run_probtrackx(seedMask, targMask, bedpBase, brainMask, outDir, 
                           doSeedNorm=True, doSize=True, 
                           doTargMaskedFDT=True, 
                           ccStop=False, 
                           bRedo=args.bRedo,
                           logFN=logFileName)
    elif args.step == "cort_conn_mat":
        #=== Calculate the cortical connectivity matrix ===#
        #===     from the probtrackx results ===#
        from dwi_analysis_settings import SURF_CLASSIFIERS
        
        if args.parcName == "" or args.parcName == None:
            error_log("Compulsory --parc option not supplied for step %s" \
                      % args.step,
                      logFN=logFileName)
        
        if SURF_CLASSIFIERS["name"].count(args.parcName) == 0:
            error_log("Parcellation '%s' is not found in dwi_analysis_settings.py" % args.parcName,
                      logFN=logFileName)

        if args.hemi == None or args.hemi == "":
            error_log("Required option hemi is not supplied during step %s" % args.step,
                      logFN=logFileName)
        assert(HEMIS.count(args.hemi) == 1)

        if args.maskType == None or args.maskType == "":
            error_log("Required option maskType is not supplied during step %s" % args.step,
                      logFN=logFileName)

        parcIdx = SURF_CLASSIFIERS["name"].index(args.parcName)
        list_py = SURF_CLASSIFIERS["list_py"][parcIdx]

        roiList = __import__(list_py)
        roiList = roiList.aROIs

        #=== Check the existence of all masks (diffusion-space) ===#
        parcDir = os.path.join(annotDir, args.parcName)
        typeDir = os.path.join(parcDir, args.maskType)

        check_dir(parcDir, logFN=logFileName)
        check_dir(typeDir, logFN=logFileName)

        check_dir(tracksDir, logFN=logFileName)

        parcTracksDir = os.path.join(tracksDir, args.parcName)
        check_dir(parcTracksDir, logFN=logFileName)

        #=== File name of the output ===#
        check_dir(connDir, bCreate=True, logFN=logFileName)

        connFN = "%s_%s_%s" \
                 % (args.parcName, args.maskType, args.hemi)
        if args.bSpeech:
            connFN += ".speech"
        connFN += ".mat"
        connFN = os.path.join(connDir, connFN)

        from tractography import generate_cort_conn_mat
        generate_cort_conn_mat(roiList, typeDir, parcTracksDir, args.hemi, 
                               args.bSpeech,  args.maskType, connFN,
                               logFN=logFileName)

        

    else:
        error_log("Unrecognized step: %s" % args.step, logFN=logFileName)
