DWI_ANALYSIS_DIR = "/speechlab/subjects/pool_dwi/"
BASE_TRACULA_CFG = "/speechlab/software/pool_dwi/base.cfg"
# FS_SUBJECTS_DIR = "/speechlab/subjects/pool"

FIX_BVECS_SCRIPT = "/speechlab/software/pool_dwi/fix_bvecs.py"

# INIT_FLIRT_MATS must be an array
INIT_FLIRT_MATS = ["/speechlab/software/pool_dwi/_defaults/d2a.mat"]

#=== FS auto parcellation paradigms ===#
# name: name of the parcellation 
# gcs: classifier gcs file name, use {hemi} for hemisphere name
# ctab: color table name
SURF_CLASSIFIERS = {"name": ["aparc12"], \
                    "gcs": ["/speechlab/software/pool_dwi/_defaults/{hemi}.slFRSatlas17.gcs"], \
                    "ctab": ["/speechlab/software/pool_dwi/_defaults/slFRS17.ctab"]}
# In order to ensure a homogenous parcellation across all subjects, 
# the script will stipulate that the specified version of FreeSurfer
# is used for mris_ca_label
PARC_FS_VER = "5.0.0" 

TRACULA_DOEDDY = 0
TRACULA_DOROTVECS = 0
TRACULA_THR_BET = 0.3
