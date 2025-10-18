import os,shutil
import json
import copy

BtS_base="0\nvacuum\n"

folder_from="configurations"
folder_work="tests_neb"
wd=os.getcwd()

fpath=os.path.join(wd,folder_from)
wpath=os.path.join(wd,folder_work)
if os.path.exists(wpath):
    shutil.rmtree(wpath)
os.mkdir(wpath)
with open("active_bonds.json","r") as file:
    active_bonds=json.load(file)
#print(active_bonds)

"""_listfiles=os.listdir(fpath)
listfiles=[]
for name in _listfiles:
    if "sp" in name:
       listfiles.append(name)
"""

for name in active_bonds.keys():
    dir_name=os.path.join(wpath,name)
    os.mkdir(dir_name)
    shutil.copy(os.path.join(fpath,name+"-sp.xyz"), os.path.join(dir_name,"to_opt.xyz"))
    
    BtS_content=copy.deepcopy(BtS_base)
    for arr in active_bonds[name]["associations"]:
        BtS_content+=f"b {arr[0]+1} {arr[1]+1} 1\n"
    for arr in active_bonds[name]["dissociations"]:
        BtS_content+=f"b {arr[0]+1} {arr[1]+1} -1\n"
    with open(os.path.join(dir_name,"bonds_to_search"),"w+") as BtS:
        BtS.writelines(BtS_content)
    