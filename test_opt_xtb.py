import os
from pyxyz import Confpool
from TS_find import optTS
def opt_and_compare(rpath):
    optTS(rpath, "to_opt.xyz",optimized_cap=0.00001, ratio=60, maxstep=1000,mode="autostrict")
    p = Confpool()
    p.include_from_file(os.path.join(rpath,"TS.xyz"))
    p.include_from_file(os.path.join(rpath,"xtbopt.xyz"))
    p.generate_connectivity(0, mult=1.3, ignore_elements=[])
    p.generate_isomorphisms()
    rmsd_value= p[0].rmsd(p[len(p)-1])[0]
    print(rmsd_value)
    return rmsd_value

def print_result(result):
    color_by_value=lambda val :"\033"+ ("[92mgood"if val<0.001 else "[93mnot bad" if val<0.01 else "[91mbad") + "\033[00m"
    for key in result.keys():
        print(f'{key} \t{"{:9.7f}".format(result[key])} {color_by_value(result[key])}')
        

rpaths=["da_test", "ep_test","sn2_test","apw_test"]
result={}
for rpath in rpaths:
    print("   ")
    print(rpath)
    result[rpath]=opt_and_compare(os.path.join(os.getcwd(),"tests",rpath))
print_result(result)