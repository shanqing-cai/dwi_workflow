#!/speechlab/software/EPD7/epd-7.1-1-rh5-x86_64/bin/python

import os
import sys
import glob
import argparse
import tempfile
from scai_utils import *

def extract_dwi_from_dicom(inDicomDir, outDir):
    #=== Work ===#
    #== Check the availability of dcm2nii and dicomhdr ==#
    check_bin_path("dcm2nii")
    check_bin_path("dicom_hdr")
    
    check_dir(inDicomDir)

    #== Run dcm2nii ==#
    tmpDir = tempfile.mktemp()
    saydo("mkdir %s" % tmpDir)
    check_dir(tmpDir)

    dcm2nii_cmd = "dcm2nii -a n -o %s %s" \
                  % (os.path.abspath(tmpDir), \
                     os.path.join(os.path.abspath(inDicomDir), "*"))
    saydo(dcm2nii_cmd)

    #== Scan for bval and bvec files ==#
    bvals = glob.glob(os.path.join(tmpDir, "*.bval"))
    bvals.sort()

    if len(bvals) == 0:
        print("Cannot find any diffusion-weighted series in input dicom directory: %s" % inDicomDir)
        saydo("rm -rf %s" % tmpDir)

    out = {"dwi4d": [], "bvals": [], "bvecs": []}

    outDir = os.path.abspath(outDir)
    check_dir(outDir, bCreate=True)

    for (i0, bval_fn) in enumerate(bvals):
        fn_base = os.path.split(bval_fn)[1]
        fn_path = os.path.split(bval_fn)[0]
        dwi4d = os.path.join(fn_path, fn_base.replace(".bval", ".nii.gz"))
        check_file(dwi4d)
        
        bvec_fn = os.path.join(fn_path, fn_base.replace(".bval", ".bvec"))
        check_file(bvec_fn)

        saydo("cp %s %s" % (bval_fn, outDir))
        check_file(os.path.join(outDir, fn_base))
        out["bvals"].append(os.path.join(outDir, fn_base))

        saydo("cp %s %s" % (bvec_fn, outDir))
        bvec_fn_1 = os.path.join(outDir, fn_base.replace(".bval", ".bvec"))
        check_file(bvec_fn_1)
        out["bvecs"].append(bvec_fn_1)
        
        saydo("cp %s %s" % (dwi4d, outDir))
        dwi4d_fn_1 = os.path.join(outDir, \
                                  fn_base.replace(".bval", ".nii.gz"))
        check_file(dwi4d_fn_1)
        out["dwi4d"].append(dwi4d_fn_1)

    saydo("rm -rf %s" % tmpDir)

    return(out)
    

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Unpack a dicom directory and extract files that belong to the diffusion series (Requires: dcm2nii (mricron), dicom_hdr (AFNI)")
    ap.add_argument("inDicomDir", help="Input DICOM dir (contains .dcm files)")
    ap.add_argument("outDir", help="Output directory")

    if len(sys.argv) == 1:
        ap.print_help()
        sys.exit(0)

    args = ap.parse_args()
    
    out = extract_dwi_from_dicom(args.inDicomDir, args.outDir)
    
    print(out)
    
