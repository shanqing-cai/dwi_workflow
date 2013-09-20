def calculate_roi_tensor_measures(parcName, roiList, wmDepths,
                                  TENSOR_MEASURES, HEMIS,
                                  dmriDir, annotDir, 
                                  tensMeasMatFN, logFileName):
    #=== Load masks (gm and wm of different depths) ===#
    import os
    import nibabel as nb
    import numpy as np
    from scai_utils import check_file, check_dir, info_log, error_log
    
    mask_shapes = []
        
    nDepths = len(wmDepths)
    roiNames = []
    nzIdx = []
    for i0 in range(nDepths):
        nzIdx.append([])
            
    bSpeech = []

    parcDir = os.path.join(annotDir, parcName)
    for (i0, wmDepth) in enumerate(wmDepths):
        if wmDepth == -1:
            parcTypeDir = os.path.join(parcDir, "gm")
            info_log("Loading gray-matter masks")
        else:
            parcTypeDir = os.path.join(parcDir, "wm%dmm" % wmDepth)
            info_log("Loading white-matter masks of %d-mm depth" \
                     % wmDepth)
            
        for (i1, troi) in enumerate(roiList):
            for (i2, hemi) in enumerate(HEMIS):
                maskFN = os.path.join(parcTypeDir, \
                                      "%s_%s.diff.nii.gz" % (hemi, troi[0]))
                check_file(maskFN, logFN=logFileName)
        
                t_img = nb.load(maskFN)
                t_img_dat = t_img.get_data()

                mask_shapes.append(np.shape(t_img_dat))

                t_img_dat = np.ndarray.flatten(t_img_dat)
                nzIdx[i0].append(np.nonzero(t_img_dat)[0])

                if wmDepth == -1:
                    if troi[2] == 'N':
                        bSpeech.append(0)
                    else:
                        bSpeech.append(1)
                        
                    roiNames.append("%s_%s" % (hemi, troi[0]))

    #=== Check that the dimensions of all mask images match ===#
    if len(np.unique(mask_shapes)) != 1:
        error_log("Non-unique matrix size among the mask files",
                  logFN=logFN)

    #=== Load the dtifit_* files and extract the measures ===#
    nROIs = len(roiNames)
    assert(len(bSpeech) == nROIs)
        
    tensMeas = {}

    mask_shape = np.unique(mask_shapes)

    check_dir(dmriDir, logFN=logFileName)
    for (i0, measName) in enumerate(TENSOR_MEASURES):
        tensMeas[measName] = np.zeros([nROIs, nDepths])

        for (i1, t_depth) in enumerate(wmDepths):
            if t_depth == -1:
                info_log("Extracting tensor measure %s from gray matter" \
                         % measName)
            else:
                info_log("Extracting tensor measure %s from %d-mm deep white matter" % (measName, t_depth))
                
            assert(len(nzIdx[i1]) == nROIs)

            for (i2, troi) in enumerate(roiNames):
                measImg = os.path.join(dmriDir, \
                                       "dtifit_%s.nii.gz" % measName)
                check_file(measImg, logFN=logFileName)

                t_img = nb.load(measImg)
                if not list(mask_shape[0]) == list(np.shape(t_img)):
                    error_log("The diffusion tensor measure volume %s (%s) does not have a dimension that matches those of the mask files" \
                              % (measImg, measName))
                    
                t_img_dat = np.ndarray.flatten(t_img.get_data())

                tensMeas[measName][i2, i1] = \
                                       np.mean(t_img_dat[nzIdx[i1][i2]])

    #=== Write data to file ===#
    from scipy.io import savemat
    res = {"roiNames": roiNames,
           "bSpeech": bSpeech,
           "parcName": parcName,
           "maskShape": mask_shape,
           "wmDepths": wmDepths, 
           "tensMeas": tensMeas}

    savemat(tensMeasMatFN, res)
    check_file(tensMeasMatFN, logFN=logFileName)

    info_log("Tensor measures (%d types) and associated data were saved at: %s"
             % (len(TENSOR_MEASURES), tensMeasMatFN),
             logFN=logFileName)
