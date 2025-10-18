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
    
    zmin=100000
    for i in range(len(z)):
        for j in range(len(z[0])):
            if(x[i]>1.6 and y[i]>1.6 and x[i]):
                zmin=min(zmin,z[i][j])
            if z[i][j] < -12.870:
                z[i][j] = -12.870
    for i in range(len(z)):
        for j in range(len(z[0])):
            z[i][j] = (z[i][j]-zmin)*627.509-4

    return x,y,z

def read_way(waypath):
    wx,wy,wz=[],[],[]
    with open(waypath, "r") as file:
        way_strs=file.readlines()
    for w_str in way_strs:
        strsplit=w_str.split()
        wx.append(float(strsplit[0]))
        wy.append(float(strsplit[1]))
        #wz.append(float(strsplit[2]))
    return wx,wy,wz 

NAME="Sn2Cl_many_dirs"
initial_cwd = os.getcwd()
rpath=os.path.join(initial_cwd,f"{NAME}_scan")

ix,iy,z=read_scan(os.path.join(rpath,"scan_global.txt"))
x,y=np.meshgrid(ix,iy)

ways=[]
wfpath=os.path.join(os.getcwd(),f"scan_opt_{NAME}")
k=0
while 1:
    wpath=os.path.join(wfpath,f"work{k}")
    print(wpath)
    if not os.path.exists(wpath):
        break
    waypath=os.path.join(wpath,"way_log.txt")
    ways.append(read_way(waypath))
    print(len(ways[-1][0]))
    k+=1

plt.axes().set_aspect(1)  

#B_SKIP=5
NLAYERS=23
A_SKIP=13

colors_contour=[]

#for i in range(B_SKIP):
#    colors_contour.append((0,0,0))

REL_WHITE=0.7
for i in range(NLAYERS):
    colors_contour.append(((i/NLAYERS)**0.7,0.5*(i/NLAYERS)**0.7,0.5-0.5*i/NLAYERS))
    #colors_contour[-1]=(colors_contour[-1][0]*REL_WHITE+1-REL_WHITE, colors_contour[-1][1]*REL_WHITE+1-REL_WHITE, colors_contour[-1][2]*REL_WHITE+1-REL_WHITE)
    

import matplotlib.colors as mcolors
cmap_custom = mcolors.ListedColormap(colors_contour)

contour_plot=plt.contourf(x+100, y, z,NLAYERS,zorder=2,cmap=cmap_custom)#this plot only used for colorbar 
plt.colorbar(contour_plot,label="kcal/mol")

plt.contour(x,y,z,NLAYERS+A_SKIP,zorder=1,linewidths=0.5,colors=colors_contour)

for i,way in enumerate(ways):
    if(1):
        plt.plot(way[1],way[0],color=(1,1,1),linewidth=3,zorder=2*i+2)
        plt.plot(way[1],way[0],color=(0.7*np.sin(np.pi*(0.333+i/len(ways)))**2, 0.7*np.sin(np.pi*i/len(ways))**2,0.7*np.sin(np.pi*(0.666+i/len(ways)))**2 ),linewidth=1,zorder=2*i+2)
        plt.scatter(way[1][0],way[0][0],color=(0,0,0),linewidth=1,zorder=2*i+3)
        if(i!=7):
            plt.scatter(way[1][-1],way[0][-1],color=(1,0,0),marker="o",sizes=[25],linewidth=3,zorder=2*i+3)


plt.xlim(x.min(), 2.2) 
plt.ylim(1.8, 2.5) 
plt.xlabel("forming bond")
plt.ylabel("breaking bond")

plt.savefig(f"pictures/fig_scan_{NAME}", dpi=300)
