def correctbvec4fsl(dwifile, bvec, bvec1, rotMat, logFN=None):
# Based on the script (by Satrajit Ghosh?): qc2rotatefsl.py

    import nibabel as nib
    import numpy as np
    from nipype.utils.filemanip import split_filename
    from scai_utils import info_log, error_log
 
    aff = nib.load(dwifile).get_affine()[:3, :3]
    for i in range(10):
        #aff = aff.dot(np.linalg.inv(np.eye(3) + 3*aff.T.dot(aff)).dot(3*np.eye(3) + aff.T.dot(aff)))
        aff = 0.5 * (aff + np.linalg.inv(aff.T))

    # DTIPrep output in nifti
    # mat = np.dot(aff, np.array([[1,0,0],[0,1,0],[0,0,-1]])) # STUT, RHY
    # mat = np.dot(aff, np.array([[1,0,0],[0,-1,0],[0,0,1]])) # CAT
    mat = np.dot(aff, np.array(rotMat)) # CAT

    bvecs = np.genfromtxt(bvec)
    if bvecs.shape[1] != 3:
        bvecs = bvecs.T
    bvecs = mat.dot(bvecs.T).T
    # outfile = '%s_forfsl.bvec' % split_filename(bvec)[1]
    outfile = bvec1
    np.savetxt(outfile, bvecs, '%.17g %.17g %.17g')
    return outfile

def format_bvals_bvecs(bvalsFN, bvecsFN, logFN=None):
    from scai_utils import cmd_stdout, info_log, error_log
    import numpy as np

    #=== Check the format of bvals and bvecs ===#

    import numpy as np

    (so, se) = cmd_stdout("wc -l %s" % bvecsFN)
    if len(se) > 0 or len(so) == 0:
        error_log("Cannot perform wc on bvecs file: %s" % bvecsFN, logFN=logFN)
    
    ln = int(so.split(" ")[0])
    info_log("ln = %d" % ln, logFN=logFN)
        
    if ln < 3:
        error_log("Unrecognized format in bvecs file: %s" % bvecsFN, \
                  logFN=logFN)
    elif ln == 3:
        #== Convert bvecs file ==#
        bvecs = np.genfromtxt(bvecsFN)
            
        assert(len(bvecs) == ln)
        bvecs = bvecs.T
            
        np.savetxt(bvecsFN, bvecs, fmt="%.15f")
        
        bvecs = np.genfromtxt(bvecsFN)
        lbv = len(bvecs)
        assert(lbv > 3)
        
        info_log("INFO: Swapped rows and columns in bvecs file: %s\n" \
                 % bvecsFN,
                 logFN=logFN)

        #== Convert bvals file ==#
        bvals = np.genfromtxt(bvalsFN).T
        np.savetxt(bvalsFN, bvals, fmt="%.15f")
        info_log("INFO: Swapped rows and columns in bvecs file: %s\nc" \
                 % bvecsFN,
                 logFN=logFN)
