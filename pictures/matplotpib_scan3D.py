import numpy as np
import matplotlib.pyplot as plt
import os

def read_scan(scanpath):
    with open(scanpath, "r") as scanfile:
        file_strs=scanfile.readlines()
    x_split=file_strs[0].split()
    x=np.linspace(float(x_split[0]),float(x_split[1]), 1+int(x_split[2]))

    y_split=file_strs[1].split()
    y=np.linspace(float(y_split[0]),float(y_split[1]), 1+int(y_split[2]))

    z=[]
    for line in file_strs[2:]:
        line=line[:-1]
        linesplit=line.split(" ")
        linesplit=[x for x in linesplit if (x and x!=" " and x!="\n")]
        if linesplit!=[]:
            z.append([])
            for E_val in linesplit:
                z[len(z)-1].append(float(E_val))
    return x,y,z

def read_way(waypath):
    wx,wy,wz=[],[],[]
    with open(waypath, "r") as file:
        way_strs=file.readlines()
    for w_str in way_strs:
        strsplit=w_str.split()
        wx.append(float(strsplit[0]))
        wy.append(float(strsplit[1]))
        wz.append(float(strsplit[2]))
    return wx,wy,wz 

NAME="sn2Cl"
initial_cwd = os.getcwd()

ways=[]
wfpath=os.path.join(os.getcwd(),f"scan_opt_{NAME}")
k=0
while 1:
    wpath=os.path.join(wfpath,f"work{k}")
    print(wpath)
    if not os.path.exists(wpath):
        break
    if 1:
        waypath=os.path.join(wpath,"way_log.txt")
        ways.append(read_way(waypath))
    k+=1

fig = plt.figure()

ax = fig.add_subplot(projection='3d')
for i,way in enumerate(ways):
    if(1):
        ax.plot(way[1],way[0],way[2],color=(i/24,0,0),linewidth=1,zorder=2*i+2)
        ax.scatter(way[1][0],way[0][0],way[2][0],color=(0,0,0),linewidth=1,zorder=2*i+3)
        ax.scatter(way[1][-1],way[0][-1],way[2][-1],color=(0,0,1),marker="+",sizes=[100],linewidth=3,zorder=2*i+3)

plt.xlim(1.5, 3.5) 
plt.ylim(1.5, 3.5)
plt.xlabel("     d(C, O)",loc="right")
plt.ylabel("               d(C, Cl)",loc="top")
plt.xticks(np.arange(1.5, 3.5, 0.5))
plt.yticks(np.arange(1.5, 3.5, 0.5))

ax.set_zlim(3.14/4,3.14)
ax.set_zlabel("a(Cl, C, O)")
fig.savefig(f"pictures/fig_3d_{NAME}", dpi=600)
