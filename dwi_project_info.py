#==========================================================#
# rawFormat: supported: 
#    DICOM
#       DICOM file with or without accompanying bvals & bvecs files
#       If the bvals & bvecs files are unavailable, leave bvalsPath and bvecsPath empty.
#       The program will extract the bvals and bvecs.
#
#    ngz
#       .nii.gz files with accompanying bvals and bvecs files
#
#==========================================================#
# fsSubjIDs: FreeSurfer subject IDs
#
#==========================================================#
# rotMat is the matrix used in the postqc step
#    It is a 3x3 matrix with 1's or -1's along the diagonal.
#    In most circumstances, it should be identical across all subjects in a given project.
#    At this point it can only be figured out by trial and error.
#    Shanqing will work on a more satisfactory solution later.
#
#==========================================================#

import numpy as np

projInfo = {"name": ["CAT",
                     "STUT",
                     "innerspeech",
                     "SEQPDS",
                     "RHY",
					 "FRS",
					 "CCRS",
					 "SEQ"], \
            "subjIDs": [["12", "13", "14", "15", "20", "21", "22", "24",
                         "28", "29", "38", "45", "46", "47", "52", "56",
                         "58"], \
                        ["S01", "S02", "S03", "S04", "S05",
                         "S06", "S07", "S08", "S09", "S10",
                         "S11", "S12", "S13", "S14", "S15",
                         "S16", "S17", "S18", "S19", "S20",
                         "S21", "S22", "S23", "S24", "S25",
                         "S26", "S27", "S28", "S29", "S30",
                         "S31", "S32", "S33", "S34", "S35",
                         "S36", "S37", "S38", "S39"], \
                        ["S10", "S11", "S12", "S13",
                         "S15", "S16", "S18", "S19",
                         "S20", "S21", "S24"], \
                        ["SEQ01P01", "SEQ03P03", "SEQ04P04", "SEQ01P05", "SEQ02P07",
                           "SEQ03P08", "SEQ04P09", "SEQ02P10", "SEQ01P11", "SEQ02P12",
                           "SEQ03P13", "SEQ04P14", "SEQ01P15", "SEQ02P16", "SEQ03P17",
                           "SEQ04P18",
                           "SEQ03C01", "SEQ04C02", "SEQ01C03", "SEQ03C04", "SEQ04C05",
                           "SEQ01C06", "SEQ03C07", "SEQ02C08", "SEQ01C10", "SEQ01C11",
                           "SEQ01C12", "SEQ03C13", "SEQ02C14", "SEQ04C15", "SEQ04C16",
                           "SEQ02C17"], \
                        ["ANS_M01"], \
						["FRS03001", "FRS03003", "FRS03005", "FRS03006", "FRS03007", 
						"FRS03008", "FRS03011", "FRS03016"], \
						["CCRS03015", "CCRS03017", "CCRS03009", "CCRS03010", "CCRS03013"], \
						["SEQ03001", "SEQ03003", "SEQ03005", "SEQ03006",
						"SEQ03009", "SEQ03010", "SEQ03011", "SEQ03012", "SEQ03013",
						"SEQ03015", "SEQ03018"]], \
            "fsSubjIDs": [["CAT_12", "CAT_13", "CAT_14", "CAT_15", "CAT_20", "CAT_21", "CAT_22", "CAT_24",
                           "CAT_28", "CAT_29", "CAT_38", "CAT_45", "CAT_46", "CAT_47", "CAT_52", "CAT_56",
                           "CAT_58"], \
                          ["STUT_S01", "STUT_S02", "STUT_S03", "STUT_S04", "STUT_S05",
                           "STUT_S06", "STUT_S07", "STUT_S08", "STUT_S09", "STUT_S10",
                           "STUT_S11", "STUT_S12", "STUT_S13", "STUT_S14", "STUT_S15",
                           "STUT_S16", "STUT_S17", "STUT_S18", "STUT_S19", "STUT_S20",
                           "STUT_S21", "STUT_S22", "STUT_S23", "STUT_S24", "STUT_S25",
                           "STUT_S26", "STUT_S27", "STUT_S28", "STUT_S29", "STUT_S30",
                           "STUT_S31", "STUT_S32", "STUT_S33", "STUT_S34", "STUT_S35",
                           "STUT_S36", "STUT_S37", "STUT_S38", "STUT_S39"], \
                          ["innerspeech.03", "innerspeech.04", "innerspeech.05", "innerspeech.06",
                           "innerspeech.08", "innerspeech.09", "innerspeech.11", "innerspeech.12",
                           "innerspeech.13", "innerspeech.14", "innerspeech.16"], \
                          ["1", "2", "3", "4", "5",
						  "6", "7", "8", "9", "10",
                         "11", "12", "13", "14", "15",
                         "16", "17", "18", "19", "20",
                         "21", "22", "23", "24", "25",
                         "26", "27", "28", "29", "30",
                         "31", "32"], \
                          ["RHY_ANS_M01"], \
						  ["1", "3", "5", "6", "7",
						  "8", "11", "16"], \
						  ["1", "4", "12", "13", "15"], \
						  ["1", "3", "4", "5",
						  "6", "7", "8", "9", "10",
                         "11", "12"]], \
            "rawFormat": ["DICOM",
                          "ngz",
                          "DICOM",
                          "DICOM",
						  "DICOM",
						  "DICOM",
						  "DICOM",
                          "DICOM"],
            "rawPath": ["/speechlab/5/carrie/cat/dti/{subjID}",
                        "/speechlab/5/scai/BACKUP/DATA/{subjID}/diffusion/*.nii.gz",
                        "/speechlab/3/overduin/data/2???????/{subjID}/mri",
                        "/speechlab/home/dsbeal/SEQPDS_orig/{subjID}/TrioTim*",
                        "/speechlab/5/scai/RHY/RAWDATA/{subjID}",
						"/speechlab/2/jsegawa/FRS03/rawdata/{subjID}/TrioTim*",
						"/speechlab/4/disk2/jsegawa/CCRS03/rawdata/{subjID}/TrioTim*",
						"/speechlab/2/jsegawa/SEQ03_pass1/RAWDATA/{subjID}/TrioTim*"],
            "bvalsPath": ["/speechlab/5/carrie/cat/dti/{subjID}/bvals",
                          "/speechlab/5/scai/BACKUP/DATA/{subjID}/diffusion/*.bvals",
                          "", 
						  "",
						  "",
						  "",
                          "",
                          ""], \
            "bvecsPath": ["/speechlab/5/carrie/cat/dti/{subjID}/bvecs",
                          "/speechlab/5/scai/BACKUP/DATA/{subjID}/diffusion/*.bvecs",
                          "",
						  "",
						  "",
						  "",
                          "", 
                          ""], \
            "rotMat": [[np.array([[1,0,0],[0,-1,0],[0,0,1]])] * 17, 
                       [np.array([[1,0,0],[0,1,0],[0,0,-1]])] * 39, 
                       [np.array([[1,0,0],[0,-1,0],[0,0,-1]])] * 11, 
                       [np.array([[1,0,0],[0,1,0],[0,0,-1]])] * 32, 
                       [np.array([[-1,0,0],[0,1,0],[0,0,1]])] * 1,
					   [np.array([[1,0,0],[0,1,0],[0,0,-1]])] * 32,
					   [np.array([[1,0,0],[0,1,0],[0,0,-1]])] * 32,
					   [np.array([[1,0,0],[0,1,0],[0,0,-1]])] * 32]}
