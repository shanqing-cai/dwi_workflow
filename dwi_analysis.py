#!/speechlab/software/EPD7/epd-7.1-1-rh5-x86_64/bin/python

import os
import sys
import glob
import argparse
from scai_utils import *

from dwi_project_info import projInfo
from dwi_analysis_settings import DWI_ANALYSIS_DIR

from scan_for_dwi import extract_dwi_from_dicom

ALL_STEPS = {"convert", "dtiprep"}

if __name__ == "__main__":
    stepsHelp = ""
    for str in ALL_STEPS:
        stepsHelp += "%s | " % str
    stepsHelp = stepsHelp[:-3]
    
    ap = argparse.ArgumentParser(description="DWI workflow")
    ap.add_argument("projName", help="Project name (e.g., CAT, SEQPDS, STUT)")
    ap.add_argument("subjID", help="Subject ID")
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
    
    #==== Main branches ====#
    if args.step == "convert":
        #=== Run DWIConvert ===#
        #== Check the path to DWIConvert ==#
        (so, se) = cmd_stdout("which DWIConvert")
        if len(se) > 0 or len(so) == 0:
            raise Exception, "Cannot find the path to executable: DWIConvert"

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

    else:
        raise Exception, "Unrecognized step: %s" % args.step
