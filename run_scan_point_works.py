import os
from TS_find_mirror import optTS

wfpath=os.path.join(os.getcwd(),"scan_opt_sn2Cl")
k=0
while 1:
    wpath=os.path.join(wfpath,f"work{k}")
    print(wpath)
    if not os.path.exists(wpath):
        break
    os.chdir(wpath)
    if(k>0):
        optTS(xyz_path=os.path.join("to_opt.xyz"), threshold_rel=8, threshold_force=0.0001, mirror_coef=0.4, print_output=True, maxstep=10**4, programm=dict(name="xtb", force_constant= 6, acc=0.0001),do_preopt=True,step_along=0)
    k+=1