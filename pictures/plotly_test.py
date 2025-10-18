import numpy as np
import pandas as pd
import plotly.graph_objects as go

fig = go.Figure(data=[
    go.Mesh3d(
        x=[0, 0, 1, 1, 0, 0, 1, 1],
        y=[0, 1, 1, 0, 0, 1, 1, 0],
        z=[0, 0, 0, 0, 1, 1, 1, 1],
        color='gray',
        i=[7, 0, 0, 0, 4, 4, 6, 6, 4, 0, 3, 2],
        j=[3, 4, 1, 2, 5, 6, 5, 2, 0, 1, 6, 3],
        k=[0, 7, 2, 3, 6, 7, 1, 1, 5, 5, 7, 6],
        opacity=0.5,
        hoverinfo='skip'
    ),
    go.Scatter3d(
        x=[0.2, 0.3, 0.4],
        y=[0.2, 0.3, 0.4],
        z=[0.2, 0.3, 0.4],
        mode="markers",
        hoverinfo='text',
        text=['1', '2', '3'],
        marker=dict(size=5, symbol="circle",
                    )
    ),
])
fig.update_layout()
fig.write_html("fig.html")