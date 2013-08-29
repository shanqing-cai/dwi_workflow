#!/speechlab/software/EPD7/epd-7.1-1-rh5-x86_64/bin/python

import os
import sys
import argparse

import nibabel as nib
import numpy as np

from scai_utils import *

def get_trackvis_mat(aff):
    """based on input from ruopeng @ tracvkis

    """
    for i in range(10):
#aff = aff.dot(np.linalg.inv(np.eye(3) + 3*aff.T.dot(aff)).dot(3*np.eye(3) + aff.T.dot(aff)))
        aff = 0.5 * (aff + np.linalg.inv(aff.T))

    return np.dot(aff, np.array([[1,0,0],[0,-1,0],[0,0,1]]))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Fix the erroneous bvecs (gradient vectors) using the solution provided by the TrackVis folks.")
    ap.add_argument("imgfile", help="Input raw 4D DWI image")
    ap.add_argument("bvecs_in", help="Input bvecs file")
    ap.add_argument("bvecs_out", help="Output bvecs file")

    if len(sys.argv) == 1:
        ap.print_help()
        sys.exit(0)

    args = ap.parse_args()
    imgfile = args.imgfile
    bvecs_in = args.bvecs_in
    bvecs_out = args.bvecs_out
    
    check_file(imgfile)
    check_file(bvecs_in)

    img = nib.load(imgfile)
    bvecs = np.genfromtxt(bvecs_in)
        
    aff = img.get_affine().copy()[:3, :3]
    mat = get_trackvis_mat(aff)

    bvecs1 = mat.dot(bvecs.T).T
    
    np.savetxt(bvecs_out, bvecs1, fmt="%.6f")
    check_file(bvecs_out)
    
    print("Saved fixed bvecs to file: %s" % bvecs_out)
