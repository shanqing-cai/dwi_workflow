def run_probtrackx(seedMask, targMask, bedpBase, brainMask, outDir, 
                   doSeedNorm=True, doSize=True, 
                   doTargMaskedFDT=True, 
                   ccStop=False, 
                   bRedo=False):
#=========================================================#
# Mode 1: from seed to targ
#         Specify both seedMask and targMask
#
# Mode 2: from seed to all
#         Specify only seedMask; set targMask=None
#
# Options:
#         ccStop: Use corpus callosum stop mask
#
#=========================================================#
    import os
    from scai_utils import check_file, check_dir, check_bin_path, \
                           saydo, cmd_stdout
    from mri_utils import nz_voxels
    
    #== Get seed and targ nvox ==#
    check_file(seedMask)

    (seed_nVoxels, seed_mm3) = nz_voxels(seedMask)
    seed_nVoxels = float(seed_nVoxels)
    assert(seed_nVoxels > 0)
    
    if targMask != None:
        check_file(targMask)
        
        (targ_nVoxels, targ_mm3) = nz_voxels(targMask)
        targ_nVoxels = float(targ_nVoxels)
        assert(targ_nVoxels > 0)
        
    check_bin_path("probtrackx")
        
    check_dir(outDir)

    if targMask != None:
        #= Prepare waypoint file =#
        wpfn = os.path.join(outDir, "waypoints.txt")
        wptext = os.path.abspath(targMask) + "\n"

        wpf = open(wpfn, "w")
        wpf.write(wptext)
        wpf.close()
        check_file(wpfn)
        
    cmd = 'probtrackx --mode=seedmask -x %s ' % seedMask + \
          '-s %s ' % bedpBase + \
          '-m %s ' % brainMask + \
          '-l -c 0.2 -S 2000 --steplength=0.5 ' + \
          '-P 5000 ' + \
          '--forcedir --opd --pd --dir=%s ' % outDir

    if targMask != None:
        cmd += "--stop=%s --waypoints=%s " % (targMask, wpfn)
        
    fdt_paths_fn = os.path.join(outDir, "fdt_paths.nii.gz")

    if not os.path.isfile(fdt_paths_fn) or bRedo:
        saydo(cmd)

    #== Check for probtrackx completion ==#
    check_file(fdt_paths_fn)

    #== Save probtrackx command ==#
    cmd_fn = os.path.join(outDir, "command.txt")
    cmd_f = open(cmd_fn, "wt")
    cmd_f.write("%s\n" % cmd)
    cmd_f.close()
    check_file(cmd_fn)

    #== Generate seed size-normalized fdt_paths ==#
    fdt_paths_norm_fn = os.path.join(outDir, "fdt_paths_norm.nii.gz")
    check_bin_path("fslmaths")

    norm_cmd = "fslmaths -dt float %s -div %d %s -odt float" % \
               (fdt_paths_fn, seed_nVoxels, fdt_paths_norm_fn)
    if not os.path.isfile(fdt_paths_norm_fn) or bRedo:
        saydo(norm_cmd)
        
    check_file(fdt_paths_norm_fn)

    if doSize:
        #== Write to seed size file ==#
        seed_size_fn = os.path.join(outDir, 'seed_size.txt')
        seed_size_f = open(seed_size_fn, 'w')
        seed_size_f.write("%d %f" % (int(seed_nVoxels), seed_mm3))
        seed_size_f.close()
        check_file(seed_size_fn)

        print("INFO: Saved seed size data to file: %s" % seed_size_fn)

        if targMask != None:
            #== Write to targ size file ==#
            targ_size_fn = os.path.join(outDir, 'targ_size.txt')
            targ_size_f = open(targ_size_fn, 'w')
            targ_size_f.write("%d %f" % (int(targ_nVoxels), targ_mm3))
            targ_size_f.close()
            check_file(targ_size_fn)

            print("INFO: Saved targ size data to file: %s" % targ_size_fn)

    if (targMask != None) and doTargMaskedFDT:
        #== Get target masked tract density ==#
        check_bin_path("fslstats")
        (so, se) = cmd_stdout("fslstats %s -k %s -m" \
                              % (fdt_paths_norm_fn, targMask))
        assert(len(se) == 0)
        so = so.split()
        assert(len(so) >= 1)
        targ_masked_norm_fdt = float(so[0])

        targ_masked_norm_fdt_fn = \
            os.path.join(outDir, "targ_masked_norm_fdt.txt")
        tmnff = open(targ_masked_norm_fdt_fn, "wt")
        tmnff.write("%f" % targ_masked_norm_fdt)
        tmnff.close()

        check_file(targ_masked_norm_fdt_fn)

        print("INFO: Saved target-masked normalized FDT value tofile: %s" \
              % targ_masked_norm_fdt_fn)

def check_probtrackx_complete(trackResDir, mode, doSeedNorm=True, doSize=True):
# Mode: seedOnly or seedTarg
    import os
    from scai_utils import check_file
    
    ALL_MODES = ["seedOnly", "seedTarg"]
    assert(ALL_MODES.count(mode) == 1)

    fdt = os.path.join(trackResDir, "fdt_paths.nii.gz")
    check_file(fdt)
    
    if doSeedNorm:
        fdt_norm = os.path.join(trackResDir, "fdt_paths_norm.nii.gz")
        check_file(fdt_norm)

    if doSize:
        seed_size_fn = os.path.join(trackResDir, "seed_size.txt")
        check_file(seed_size_fn)

        if mode == "seedTarg":
            targ_size_fn = os.path.join(targResDir, "targ_size.txt")
            check_file(targ_size_fn)

    

def generate_cort_conn_mat(roiList, parcTypeDir, parcTracksDir, hemi, 
                           arg_bSpeech, maskType, connFN):
    import os
    import numpy as np
    import nibabel as nb
    from scai_utils import check_file, check_dir
    

    mask_shapes = []
    roiNames = []
    nzIdx = []
    bSpeech = []
    for (i0, troi) in enumerate(roiList):
        targROI = troi[0]
        maskFN = os.path.join(parcTypeDir, \
                              "%s_%s.diff.nii.gz" % (hemi, targROI))
        check_file(maskFN)
        
        t_img = nb.load(maskFN)
        t_img_dat = t_img.get_data()

        mask_shapes.append(np.shape(t_img_dat))
            
        t_img_dat = np.ndarray.flatten(t_img_dat)
            
        nzIdx.append(np.nonzero(t_img_dat)[0])
        roiNames.append(troi[0])
        if troi[2] == 'N':
            bSpeech.append(0)
        else:
            bSpeech.append(1)

    roiNames = np.array(roiNames)
    bSpeech = np.array(bSpeech)
    nzIdx = np.array(nzIdx)
    if arg_bSpeech:
        roiNames = roiNames[np.nonzero(bSpeech)[0]]
        nzIdx = nzIdx[np.nonzero(bSpeech)[0]]


    nROIs = len(roiNames)
    assert(len(nzIdx) == nROIs)
    if len(np.unique(mask_shapes)) != 1:
        raise Exception, "Non-unique matrix size among the mask files"
    imgShape = np.unique(mask_shapes)[0]

    #=== Check the completion of seed-only probtrackx ===#
    #===     and calculate the conn matrix ===#
    connMat = np.zeros([nROIs, nROIs])

    for (i0, troi) in enumerate(roiNames):
        seedROI = troi
        trackResDir = os.path.join(parcTracksDir, 
                                   "%s_%s_%s" % \
                                   (hemi, seedROI, maskType))
                                   
        check_probtrackx_complete(trackResDir, "seedOnly", 
                              doSeedNorm=True, doSize=True)
        
        fdt_norm = os.path.join(trackResDir, "fdt_paths_norm.nii.gz")
        t_img = nb.load(fdt_norm)
        t_img_dat = t_img.get_data()
            
        assert(list(np.shape(t_img_dat)) == list(imgShape))
        t_img_dat = np.ndarray.flatten(t_img_dat)

        for (i1, troi1) in enumerate(roiNames):
            connMat[i0, i1] = np.mean(t_img_dat[nzIdx[i1]])
        
    #=== Make symmetric ===#
    connMat = 0.5 * (connMat + connMat.T)

    #=== Write result .mat file ===#
    from scipy.io import savemat
    res = {"roiNames": roiNames, "connMat": connMat}
    savemat(connFN, res)
    check_file(connFN)
        
    print("INFO: Connectivity matrix and associated data were saved at: %s" % (connFN))
