import os
rpath=os.path.join(os.getcwd(),"sn2Cl_scan")
structsname="scan_structs_points.xyz"
bts_content="1\nwater\nb 1 8 1\nb 1 5 -1"
with open(os.path.join(rpath,structsname), "r") as xyzlog:
    lines_xyzs=xyzlog.readlines()
num=int(lines_xyzs[0])

wpath=os.path.join(os.getcwd(),"scan_opt_sn2")
if not os.path.exists(wpath):
    os.mkdir(wpath)

n=0#string number in lines_xyzs
k=0#xyz number
while n+num+2<=len(lines_xyzs) and k<500:#k<500 just for not folder overflow
    xyz_lines=lines_xyzs[n:n+num+2]
    path_wf=os.path.join(wpath,f"work{k}")
    if not os.path.exists(path_wf):
        os.mkdir(path_wf)
    path_w=os.path.join(path_wf,"to_opt.xyz")
    with open(path_w,"w+") as work:
        work.writelines(xyz_lines)
    path_bts=os.path.join(path_wf,"bonds_to_search")
    with open(path_bts,"w+") as bts:
        bts.writelines(bts_content)
    k+=1
    n+=num+2
    
    