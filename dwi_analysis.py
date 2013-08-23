#!/speechlab/software/EPD7/epd-7.1-1-rh5-x86_64/bin/python

import os
import sys
import glob
import argparse
from scai_utils import *

from dwi_project_info import projInfo
from dwi_analysis_settings import DWI_ANALYSIS_DIR, \
    BASE_TRACULA_CFG, FS_SUBJECTS_DIR, \
    FIX_BVECS_SCRIPT, \
    TRACULA_DOEDDY, TRACULA_DOROTVECS, TRACULA_THR_BET

from scan_for_dwi import extract_dwi_from_dicom

ALL_STEPS = {"convert", "dtiprep", "postqc", \
             "tracula_prep", "tracula_bedp"}

if __name__ == "__main__":
    stepsHelp = ""
    for str in ALL_STEPS:
        stepsHelp += "%s | " % str
    stepsHelp = stepsHelp[:-3]
    
    ap = argparse.ArgumentParser(description="DWI workflow")
    ap.add_argument("projName", help="Project name (e.g., CAT, SEQPDS, STUT)")
    ap.add_argument("subjID", help="Subject ID")
    ap.add_argument("--fsSubjID", dest="fsSubjID", default="", \
                    required=False, \
                    help="FreeSurfer subject ID (required by the tracula_prep step)")
    ap.add_argument("step", help="Steps: {%s}" % stepsHelp)

    if len(sys.argv) == 1:
        ap.print_help()
        sys.exit(0)

    args = ap.parse_args()

    #=== Locate entry in projInfo ===#
    pidx = projInfo["name"].index(args.projName)
    sidx = projInfo["subjIDs"][pidx].index(args.subjID)

    rawFormat = projInfo["rawFormat"][pidx].upper()
    rawFN = projInfo["rawPath"][pidx]

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
        raise Exception, "Unrecognized raw input format: %s" % rawFormat

    #== Create directory ==#
    check_dir(DWI_ANALYSIS_DIR)
    sID = "%s_%s" % (args.projName, args.subjID)
    sDir = os.path.join(DWI_ANALYSIS_DIR, sID)
    check_dir(sDir, bCreate=True)

    #== Prepare directories and some file names ==#
    dicomDir = os.path.join(sDir, "dicom")
    check_dir(dicomDir, bCreate=True)    
    
    dtiprepDir = os.path.join(sDir, "dtiprep")
    check_dir(dtiprepDir, bCreate=True)
        
    outputNrrd = os.path.join(dtiprepDir, "raw_dwi.nrrd")

    qcedNrrd = os.path.join(dtiprepDir, "raw_dwi_QCed.nrrd")
    qcXML = os.path.join(dtiprepDir, "raw_dwi_XMLQCResult.xml")
    qcReport = os.path.join(dtiprepDir, "raw_dwi_QCReport.txt")

    qcedNGZ = os.path.join(dtiprepDir, "dwi_qced.nii.gz")
    qcedBVals = os.path.join(dtiprepDir, "dwi_qced.bvals")
    qcedBVecs = os.path.join(dtiprepDir, "dwi_qced.bvecs")

    #== Locate the bvals and bvecs files ==#
    if len(projInfo["bvalsPath"][pidx]) > 0:
        bvalsFN = \
            projInfo["bvalsPath"][pidx].replace("{subjID}", "%s") % args.subjID
        bvalsFN = glob.glob(bvalsFN)
        assert(len(bvalsFN) == 1)
        bvalsFN = bvalsFN[0]
        check_file(bvalsFN)

        bvecsFN = \
            projInfo["bvecsPath"][pidx].replace("{subjID}", "%s") % args.subjID
        bvecsFN = glob.glob(bvecsFN)
        assert(len(bvecsFN) == 1)
        bvecsFN = bvecsFN[0]
        check_file(bvecsFN)
    else:
        #== DICOM files without accompanying bvals and bvecs files ==#
        if args.step == "convert":
            out = extract_dwi_from_dicom(rawFN, dicomDir)

            assert(len(out["dwi4d"]) == len(out["bvals"]))
            assert(len(out["dwi4d"]) == len(out["bvecs"]))

            if len(out["dwi4d"]) == 0:
                raise Exception, "Cannot find any DWI series in input DICOM directory: %s" % rawFN
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

        check_file(rawFN)
        check_file(bvalsFN)
        check_file(bvecsFN)

    #==== Main branches ====#
    if args.step == "convert":
        #=== Run DWIConvert ===#
        #== Check the path to DWIConvert ==#
        check_bin_path("DWIConvert")
        #(so, se) = cmd_stdout("which DWIConvert")
        #if len(se) > 0 or len(so) == 0:
        #    raise Exception, "Cannot find the path to executable: DWIConvert"

        if rawFormat == "DICOM":
            #== Copy the dicoms ==#    
            saydo("cp -v %s %s" % (os.path.join(rawFN, "*.dcm"), dicomDir))

            #== Copy the bvecs and bvals ==#
            saydo("cp -v %s %s" % (bvalsFN, dicomDir))
            bvalsFN = os.path.join(dicomDir, os.path.split(bvalsFN)[-1])

            saydo("cp -v %s %s" % (bvecsFN, dicomDir))
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
            raise Exception, "Step %s: Unsupported input format %s" % \
                             (args.step, rawFormat)

        saydo(dwiCvtCmd)
        check_file(outputNrrd)

        print("outputNrrd = %s" % outputNrrd)
            
    elif args.step == "dtiprep":
        #=== Run dtiprep ===#
        check_file(outputNrrd)

        #== Check the path to DTIPrep ==#
        (so, se) = cmd_stdout("which DTIPrep")
        if len(se) > 0 or len(so) == 0:
            raise Exception, "Cannot find the path to executable: DTIPrep"

        logFN = os.path.join(dtiprepDir, "dtiprep_notes.txt")
        dtiprepCmd = "DTIPrep -w %s -p default -d -c -n %s" % \
                     (outputNrrd, logFN)

        saydo(dtiprepCmd)
        check_file(qcedNrrd)
        check_file(qcXML)
        check_file(qcReport)

    elif args.step == "postqc":
        #=== Convert the qc'ed nrrd to nifti ===#
        check_bin_path("DWIConvert")

        qcedBVecs_unc = os.path.join(dtiprepDir, "dwi_qced.bvecs.uncorrected")

        check_file(qcedNrrd)
        cvtCmd = "DWIConvert --inputVolume %s " % qcedNrrd + \
                 "--conversionMode NrrdToFSL " + \
                 "--outputBValues %s " % qcedBVals + \
                 "--outputBVectors %s " % qcedBVecs_unc + \
                 "--outputVolume %s " % qcedNGZ

        saydo(cvtCmd)
        check_file(qcedNGZ)
        check_file(qcedBVals)
        check_file(qcedBVecs_unc)

        #=== Correct the b-vectors for FSL and FreeSurfer ===#
        #from dtiprep_utils import correctbvec4fsl
        #correctbvec4fsl(qcedNGZ, qcedBVecs_unc, qcedBVecs)
        #check_file(qcedBVecs)
    
        check_file(FIX_BVECS_SCRIPT)
        fixBVecsCmd = "%s %s %s %s" % \
                      (FIX_BVECS_SCRIPT, qcedNGZ, \
                       qcedBVecs_unc, qcedBVecs)
        check_file(qcedBVecs)

    elif args.step == "tracula_prep":
        #=== Tracula: step prep ===#
        import numpy as np
        from tracula_utils import modify_tracula_cfg_file

        #=== Check the version of FreeSurfer ===#
        fsHome = os.getenv("FREESURFER_HOME")
        (so, se) = cmd_stdout("cat %s" \
                              % os.path.join(fsHome, "build-stamp.txt"))
        assert(len(se) == 0)
        assert(len(so) > 0)
        vItem = so.replace("\n", "").split("-v")[-1]

        assert(vItem.count(".") == 2)
        vItem = float("%s.%s" % (vItem.split(".")[0], vItem.split(".")[1]))

        if vItem < 5.1:
            raise Exception, "It appears that a FreeSurfer version older than 5.1 is being used. To switch version, do . fss ver (e.g., . fss 5.3.0)"
        
        check_bin_path("trac-all")

        #=== Check input FreeSurfer subject ID ===#
        if len(args.fsSubjID) == 0:
            raise Exception, \
                "Input argument --fsSubjID must be set for step %s" % args.step

        check_file(qcedNGZ)
        check_file(qcedBVals)
        check_file(qcedBVecs)

        check_file(BASE_TRACULA_CFG)

        #=== Check SUBJECTS_DIR ===#
        fsSubjectsDir = os.getenv("SUBJECTS_DIR")
        check_dir(fsSubjectsDir)

        fsSDir = os.path.join(fsSubjectsDir, args.fsSubjID)
        check_dir(fsSDir)

        #=== Determine nb0 ===#
        bv = np.genfromtxt(qcedBVals)
        nb0 = len(np.nonzero(bv == 0.0))

        traculaCfgFN = os.path.join(sDir, "tracula.cfg")

        modify_tracula_cfg_file(fsHome, fsSubjectsDir, \
                                BASE_TRACULA_CFG, DWI_ANALYSIS_DIR, \
                                args.fsSubjID, \
                                dtiprepDir, os.path.split(qcedNGZ)[-1], \
                                qcedBVals, qcedBVecs, nb0, \
                                TRACULA_DOEDDY, TRACULA_DOROTVECS, \
                                TRACULA_THR_BET, traculaCfgFN)
        check_file(traculaCfgFN)
        print("INFO: TRACULA config file for subject %s written to: %s" % \
              (sID, traculaCfgFN))

        #=== Execute tracula ===#
        tracall_cmd = "trac-all -c %s -prep" % traculaCfgFN

        saydo(tracall_cmd)

        traculaDir = os.path.join(DWI_ANALYSIS_DIR, args.fsSubjID)
        tracula_dlabelDir = os.path.join(traculaDir, "dlabel")
        tracula_dmriDir = os.path.join(traculaDir, "dmri")
        tracula_scriptsDir = os.path.join(traculaDir, "scripts")
        
        check_dir(traculaDir)
        check_dir(tracula_dlabelDir)
        check_dir(tracula_dmriDir)
        check_dir(tracula_scriptsDir)

        #=== Move the tracula files back to the dwi directory, if necessary ===#
        if sID != args.fsSubjID:
            saydo("mv %s %s/" % (tracula_dlabelDir, sDir))
            saydo("mv %s %s/" % (tracula_dmriDir, sDir))
            saydo("mv %s %s/" % (tracula_scriptsDir, sDir))
            saydo("rmdir %s" % traculaDir)
                
    

    else:
        raise Exception, "Unrecognized step: %s" % args.step
