import sqlite3
from math import e

import plotly.graph_objs as go
import plotly.express as px
from plotly.subplots import make_subplots

import numpy as np
import pandas as pd
from scipy import signal

def frange(start, stop, step):
  i = start
  while i < stop:
    yield i
    i += step

p = []
h = []

conn = sqlite3.connect('main.db')
c = conn.cursor()

fig = go.Figure()

for x in range(200, 600, 5):
  p_l = []
  h_l = []

  c.execute('SELECT P,h FROM amm_water WHERE T=' + str(x) + ' AND P<20 ORDER BY P ASC')

  res = c.fetchall()

  for r in res:
    p_l.append(r[0] / 10)
    h_l.append(r[1])
  
  fig.add_trace(go.Scatter(x=h_l, y=p_l, name='T={0}'.format(x-273.15), mode='lines+markers', line=dict(color="blue", dash="dash"), marker=dict(size=2)))

for x in range(200, 600, 5):
  p_l = []
  h_l = []

  c.execute('SELECT P,h FROM amm_over WHERE T=' + str(x) + ' AND P<20 ORDER BY P ASC')

  res = c.fetchall()

  for r in res:
    p_l.append(r[0] / 10)
    h_l.append(r[1])
  
  fig.add_trace(go.Scatter(x=h_l, y=p_l, name='T={0}'.format(x-273.15),
                    mode='lines+markers', line=dict(color="blue", dash='dash')))

for x in frange(8.8, 10.4, 0.1):
  p_l = []
  h_l = []

  c.execute('SELECT P,h FROM amm_water WHERE s<' + str(x+0.009) + ' AND s>' + str(x-0.009) + ' AND P<20 ORDER BY P ASC')

  res = c.fetchall()

  for r in res:
    p_l.append(r[0] / 10)
    h_l.append(r[1])
  
  fig.add_trace(go.Scatter(x=h_l, y=p_l, name='S={0}'.format(x),
                    mode='lines', line=dict(color="red")))

for x in frange(8.8, 10.4, 0.1):
  p_l = []
  h_l = []

  c.execute('SELECT P,h FROM amm_over WHERE s<' + str(x+0.009) + ' AND s>' + str(x-0.009) + ' AND P<20 ORDER BY P ASC')

  res = c.fetchall()

  for r in res:
    p_l.append(r[0] / 10)
    h_l.append(r[1])
  
  fig.add_trace(go.Scatter(x=h_l, y=p_l, name='S={0}'.format(x),
                    mode='lines', line=dict(color="red")))
                    
for x in frange(0.1, 4, 0.1):
  p_l = []
  h_l = []

  c.execute('SELECT P,h FROM amm_over WHERE v<' + str(x+0.007) + ' AND v>' + str(x-0.007) + ' AND P<20 ORDER BY P ASC')

  res = c.fetchall()

  for r in res:
    p_l.append(r[0] / 10)
    h_l.append(r[1])

  fig.add_trace(go.Scatter(x=h_l, y=p_l, name='v={0}'.format(x),
                        mode='lines', line=dict(color="green")))
  
p_l = []
h_l = []
hh_l = []

c.execute('SELECT P,hd,hdd FROM amm_wet WHERE P<20 ORDER BY P ASC')

res = c.fetchall()

for r in res:
  p_l.append(r[0] / 10)
  h_l.append(r[1])
  hh_l.append(r[2])

  fig.add_trace(go.Scatter(x=[r[1], r[2]], y=[r[0] / 10, r[0] / 10], name='P,T=const',
                  mode='lines+markers', line=dict(color="blue", dash='dash')))

fig.add_trace(go.Scatter(x=h_l, y=p_l, name='x=0',
                  mode='lines', line=dict(color="black")))

fig.add_trace(go.Scatter(x=hh_l, y=p_l, name='x=1',
                  mode='lines', line=dict(color="black")))
  
fig.update_yaxes(type="log", ticklabelstep=2, showgrid=True, gridwidth=1, gridcolor='LightPink', zeroline=False, dtick = np.log10(1))
fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor='LightPink', zeroline=False, dtick = 50)

fig.show()
