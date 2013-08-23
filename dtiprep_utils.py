def correctbvec4fsl(dwifile, bvec, bvec1):
# Based on the script (by Satrajit Ghosh?): qc2rotatefsl.py

    import nibabel as nib
    import numpy as np
    from nipype.utils.filemanip import split_filename
 
    aff = nib.load(dwifile).get_affine()[:3, :3]
    for i in range(10):
        #aff = aff.dot(np.linalg.inv(np.eye(3) + 3*aff.T.dot(aff)).dot(3*np.eye(3) + aff.T.dot(aff)))
        aff = 0.5 * (aff + np.linalg.inv(aff.T))
    mat = np.dot(aff, np.array([[1,0,0],[0,1,0],[0,0,-1]])) # DTIPrep output in nifti
    bvecs = np.genfromtxt(bvec)
    if bvecs.shape[1] != 3:
        bvecs = bvecs.T
    bvecs = mat.dot(bvecs.T).T
    # outfile = '%s_forfsl.bvec' % split_filename(bvec)[1]
    outfile = bvec1
    np.savetxt(outfile, bvecs, '%.17g %.17g %.17g')
    return outfile
