
import os,shutil
import TS_find_mirror

BtS_base="0\nvacuum\n"

folder_work="scan_opt_Sn2Cl_many_dirs"
wd=os.getcwd()

wpath=os.path.join(wd,folder_work)

listworks=os.listdir(wpath)

for work in listworks:
    listfiles=os.listdir(os.path.join(wpath,work))
    for file in listfiles:
        if file != "bonds_to_search" and file != "to_opt.xyz":
            os.remove(os.path.join(wpath,work,file))
    
for work in listworks:
    try:
        os.chdir(wd)
        TS_find_mirror.optTS(os.path.join(folder_work, work, "to_opt.xyz"), threshold_rel=8, threshold_force=0.0001, print_output=True, maxstep=10**3, programm=dict(name="xtb", acc=0.0001, force_constant= 6))
    except:
        None

print ("Benchmark completed!")