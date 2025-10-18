import numpy as np
import plotly.graph_objects as go
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

initial_cwd = os.getcwd()
rpath=os.path.join(initial_cwd,"sn2Cl_scan")

ix,iy,z=read_scan(os.path.join(rpath,"scan_points.txt"))
x,y=np.meshgrid(ix,iy)
fig = go.Figure()
fig.add_surface(x=x, y=y, z=z, colorscale='Reds_r', colorbar_thickness=25, colorbar_len=0.75, opacity=0.7)
#wx,wy,wz=read_way(os.path.join(rpath,"way_log.txt"))

line_marker = dict(color='blue', width=3)

#fig.add_scatter3d(x=wx, y=wy, z=wz, mode='lines', line=line_marker, name='')



#Define the second family of coordinate lines

line_marker = dict(color='#101010', width=3)
x,y=np.meshgrid(iy,ix)
for xx, yy, zz in zip(x, y, z):
    fig.add_scatter3d(x=xx, y=yy, z=zz, mode='lines', line=line_marker, name='')

for yy, xx, zz in zip(x, y, np.array(z).T):
    fig.add_scatter3d(x=xx, y=yy, z=zz, mode='lines', line=line_marker, name='')    
    

fig.update_layout(width=700, height=700, showlegend=False)
fig.write_html("fig.html")
