#==========================================================#
# rawFormat: supported: 
#    DICOM
#       DICOM file with or without accompanying bvals & bvecs files
#       If the bvals & bvecs files are unavailable, leave bvalsPath and bvecsPath empty
#
#    ngz
#       .nii.gz files with accompanying bvals and bvecs files
#
#==========================================================#

projInfo = {"name": ["CAT", \
                     "STUT", \
                     "innerspeech", \
                     "SEQPDS"],
            "subjIDs": [["12", "13", "14", "15", "20", "21", "22", "24", \
                         "28", "29", "38", "45", "46", "47", "52", "56", \
                         "58"], \
                        ["S01", "S02", "S03", "S04", "S05", \
                         "S06", "S07", "S08", "S09", "S10", \
                         "S11", "S12", "S13", "S14", "S15", \
                         "S16", "S17", "S18", "S19", "S20", \
                         "S21", "S22", "S23", "S24", "S25", \
                         "S26", "S27", "S28", "S29", "S30", \
                         "S31", "S32", "S33", "S34", "S35", \
                         "S36", "S37", "S38", "S39"], \
                        ["S10", "S11", "S12", "S13", \
                         "S15", "S16", "S18", "S19", "S20", "S21", "S24"]],
            "rawFormat": ["DICOM", \
                          "ngz", \
                          "DICOM"],
            "rawPath": ["/speechlab/5/carrie/cat/dti/{subjID}", \
                        "/speechlab/5/scai/BACKUP/DATA/{subjID}/diffusion/*.nii.gz", \
                        "/speechlab/3/overduin/data/2???????/{subjID}/mri"], \
            "bvalsPath": ["/speechlab/5/carrie/cat/dti/{subjID}/bvals", \
                          "/speechlab/5/scai/BACKUP/DATA/{subjID}/diffusion/*.bvals", \
                          ""], \
            "bvecsPath": ["/speechlab/5/carrie/cat/dti/{subjID}/bvecs", \
                          "/speechlab/5/scai/BACKUP/DATA/{subjID}/diffusion/*.bvecs", \
                          ""]}
