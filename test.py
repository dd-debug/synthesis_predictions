'''
Created on 2020.12.3

@author: jiadongc
'''
# import pandas as pd
# import holoviews as hv
# from holoviews import opts, dim
# from bokeh.sampledata.les_mis import data
# from bokeh.plotting import show, output_file
# hv.extension('bokeh')
# hv.output(size=200)
# links = pd.DataFrame(data['links'])
# print(links)
# print(links.head(3))
# hv.Chord(links)
# nodes = hv.Dataset(pd.DataFrame(data['nodes']), 'index')
# 
# print(nodes)
# nodes.data.head()
#   
# chord = hv.Chord((links, nodes)).select(value=(5, None))
# chord.opts(
#     cmap='Category20', edge_cmap='Category20', edge_color=dim('source').str(), 
#                labels='name', node_color=dim('index').str(), width=500, height=500)
# # chord.opts(
# #     node_color='index', edge_color='source', label_index='index', 
# #     cmap='Category10', edge_cmap='Category10', width=500, height=500)
# # output_file('chordtest.html')
# # show(hv.render(chord))


# from chord import Chord
# matrix = [
#     [0, 5, 6, 4, 7, 4],
#     [5, 0, 5, 4, 6, 5],
#     [6, 5, 0, 4, 5, 5],
#     [4, 4, 4, 0, 5, 5],
#     [7, 6, 5, 5, 0, 4],
#     [4, 5, 5, 5, 4, 0],
# ]
#  
# names = ["Action", "Adventure", "Comedy", "Drama", "Fantasy", "Thriller"]
# Chord(matrix, names).show()
# Chord(matrix, names).to_html()




import plotly.graph_objects as go

fig = go.Figure()

# Create scatter trace of text labels
fig.add_trace(go.Scatter(
    x=[2, 3.5, 6],
    y=[1, 1.5, 1],
    text=["Vertical Line",
          "Horizontal Dashed Line",
          "Diagonal dotted Line"],
    mode="text",
))

# Set axes ranges
fig.update_xaxes(range=[0, 7])
fig.update_yaxes(range=[0, 2.5])

# Add shapes
fig.add_shape(type="line",
    x0=1, y0=0, x1=1, y1=2,
    line=dict(color="RoyalBlue",width=3)
)
# fig.add_shape(type="line",
#     x0=2, y0=2, x1=5, y1=2,
#     line=dict(
#         color="LightSeaGreen",
#         width=4,
#         dash="dashdot",
#     )
# )
# fig.add_shape(type="line",
#     x0=4, y0=0, x1=6, y1=2,
#     line=dict(
#         color="MediumPurple",
#         width=4,
#         dash="dot",
#     )
# )
fig.update_shapes(dict(xref='x', yref='y'))
fig.show()

import plotly.graph_objects as go
import numpy as np
t = np.linspace(0, 4*np.pi, 50)
t2 = np.pi * np.arange(5)
fig = go.Figure(go.Scatter(x=t, y=np.sin(t), mode='lines'))
fig.add_trace(go.Scatter(x=t2, y=np.sin(t2), mode='markers'))
print(np.sin(t2[0]),t2[0])
point = t2[0]
fig.update_layout(annotations=[
            go.layout.Annotation(x=point,
            y=np.sin(point),
            xref="x",
            yref="y",
            text="dict Text",
            align='center',
            showarrow=False,
            yanchor='bottom',
            textangle=90,
            font_size = 10)])
fig.show()










