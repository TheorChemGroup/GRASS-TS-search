import numpy as np
import pandas as pd
import plotly.graph_objects as go
xmin=-3.5
xmax=1.5
ymin=-4
ymax=1
xi = np.linspace(xmin, xmax, 101)
yi = np.linspace(ymin, ymax, 101)
X ,Y  = np.meshgrid(xi,yi)

xi = np.linspace(xmin, xmax, 10)
yi = np.linspace(ymin, ymax, 10)
XL,YL = np.meshgrid(xi,yi)

e=2.718281828



fig = go.Figure(data=[
    
   # go.Surface(x=X,y=Y,z=0.5*(+0.25*X**2 + e**X + 0.2*Y**2 + e**Y-1.205), opacity=0.5,colorscale=[[0, "rgb(32,120,220)"], [1, "rgb(32,120,220)"]]),
    go.Surface(x=X,y=Y,z=0.5*(-0.25*X**2 - e**X + 0.2*Y**2 + e**Y      ), opacity=0.99, colorscale=[[0, "rgb(220,220,220)"], [1, "rgb(220,220,220)"]])
])
line_marker = dict(color='#202020', width=3)
#top grid
i=0
for xx, yy in zip(X, Y):
    if i%2==0 or i==len(X)-1:
        fig.add_scatter3d(x=xx, y=yy, z=0.5*(-0.25*xx**2 - e**xx + 0.2*yy**2 + e**yy+0.01)        , mode='lines', opacity=0.2 if i!=60 else 0.99, line=line_marker, name='')
    i+=1
i=0
for xx, yy in zip(X, Y):
    if i%2==0 or i==len(X)-1:
        fig.add_scatter3d(x=xx, y=yy, z=0.5*(-0.25*xx**2 - e**xx + 0.2*yy**2 + e**yy-0.01)        , mode='lines', opacity=0.2 if i!=60 else 0.99, line=line_marker, name='')
    i+=1
i=0
for yy, xx in zip(Y.T, X.T):
    if i%2==0 or i==len(Y)-1:
        fig.add_scatter3d(x=xx, y=yy, z=0.5*(-0.25*xx**2 - e**xx + 0.2*yy**2 + e**yy+0.01)        , mode='lines', opacity=0.2, line=line_marker, name='')    
    i+=1
i=0
for yy, xx in zip(Y.T, X.T):
    if i%2==0 or i==len(Y)-1:
        fig.add_scatter3d(x=xx, y=yy, z=0.5*(-0.25*xx**2 - e**xx + 0.2*yy**2 + e**yy-0.01)        , mode='lines', opacity=0.2, line=line_marker, name='')    
    i+=1
#bottom grid
line_marker = dict(color='#2020E1', width=3)
i=0
for xx, yy in zip(X, Y):
    if i%2==0 or i==len(X)-1:
        fig.add_scatter3d(x=xx, y=yy, z=0.5*(+0.25*xx**2 + e**xx + 0.2*yy**2 + e**yy -1.205 +0.01), mode='lines', opacity=0.2 if i!=60 else 0.99, line=line_marker, name='')
    i+=1
i=0
for xx, yy in zip(X, Y):
    if i%2==0 or i==len(X)-1:
        fig.add_scatter3d(x=xx, y=yy, z=0.5*(+0.25*xx**2 + e**xx + 0.2*yy**2 + e**yy -1.205 -0.01), mode='lines', opacity=0.2 if i!=60 else 0.99, line=line_marker, name='')
    i+=1
i=0
for yy, xx in zip(Y.T, X.T):
    if i%2==0 or i==len(Y)-1:
        fig.add_scatter3d(x=xx, y=yy, z=0.5*(+0.25*xx**2 + e**xx + 0.2*yy**2 + e**yy -1.205 +0.01), mode='lines', opacity=0.2, line=line_marker, name='')    
    i+=1
i=0
for yy, xx in zip(Y.T, X.T):
    if i%2==0 or i==len(Y)-1:
        fig.add_scatter3d(x=xx, y=yy, z=0.5*(+0.25*xx**2 + e**xx + 0.2*yy**2 + e**yy -1.205 -0.01), mode='lines', opacity=0.2, line=line_marker, name='')    
    i+=1
fig.update_layout(scene={
                "camera": {
                    "projection": {
                        "type": "orthographic"
                    }
                  }})
fig.write_html("fig.html")
