def check_fs_ver(reqVer, mode="eq"):
#=== Check the version of FreeSurfer ===#
# Two modes: 
#    "eq": exact match
#    "eqgt": equal to or greater than 
    import os
    from scai_utils import cmd_stdout

    fsHome = os.getenv("FREESURFER_HOME")
    (so, se) = cmd_stdout("cat %s" \
                          % os.path.join(fsHome, "build-stamp.txt"))
    assert(len(se) == 0)
    assert(len(so) > 0)
    vItem = so.replace("\n", "").split("-v")[-1]

    assert(vItem.count(".") == 2)
    # vItem = float("%s.%s" % (vItem.split(".")[0], vItem.split(".")[1]))

    switchInstr = "To switch version, do . fss ver (e.g., . fss 5.3.0)";

    if mode == "eqgt":
        if type(reqVer) != float:
            raise Exception, "%s: reqVer must be a float value under mode %s" \
                             % (check_fs_ver.__name__, mode)

        vItem = float("%s.%s" % (vItem.split(".")[0], vItem.split(".")[1]))

        if vItem < reqVer:
            raise Exception, "It appears that a FreeSurfer version older than %.1f is being used. %s" % (reqVer, switchInstr)
    elif mode == "eq":
        if type(reqVer) != str:
            raise Exception, "%s: reqVer must be a string value under mode %s" \
                             % (check_fs_ver.__name__, mode)
        if vItem != reqVer:
            raise Exception, "It appears that a FreeSurfer version (%s) other than the required version (%s) is being used. %s" % (vItem, reqVer, switchInstr)
    else:
        raise Exception, "%s: Unrecognized mode: %s" \
                         % (check_fs_ver.__name__, mode)
            
            
