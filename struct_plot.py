#!/usr/bin/env python3
"""
-----------------------------------------------------------
 created: wwwennie
   
  make prettified plots for data associated with characterizing structure
    as found in struct_*.py scripts

  with subplot capability. maybe. haven't decided yet. 
----------------------------------------------
"""

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import seaborn as sns
import numpy as np

#================================================================================#

#------------------------------- Helper functions -------------------------------#
def color_picker(index):
  """ For when more control over colors is preferred 
      Cycles over the possible options """

  sns.set_color_codes()
  colors= ["royalblue","mediumvioletred","orange","forestgreen","black","orangered",
           "mediumslateblue","maroon","lightseagreen"]
  index = np.mod(index,len(colors))
  return colors[index]

def line_picker(index):
  """ For when more control over line styles is preferred 
      Cycles over the possible options """

  lines = ['solid','dashed','dotted','dashdot']
  index = np.mod(index,len(lines))
  return lines[index]

def marker_picker(index):
  """ For when more control over line styles is preferred 
      Cycles over the possible options """
  marker = ['.','s','v','^','*','x','D']
  index = np.mod(index,len(lines))
  return marker[index]

#------------------------------- Histogram plots -------------------------------#

def hist_plot(distributions,title,legend,xlim=[0,10],xlab="",ylab="Counts (arb. units)",bins=20,
              savefile="/home/wwwennie/Downloads/tmp.pdf"):
   """
     Make a pretty histrogram using seaborn with labels and legends

     distributions = list of lists containing distribution to plot; all distributions
                       will be plotted on same axes
     title  = title of plot
     xlabel, ylabel = axes labels
     bins = bins in histrogram
     xlim = limits of x axis

     savefile = name of file to save

   """
   sns.set(font_scale=1.5)
   sns.set_style("ticks")
   sns.set_color_codes()

   index = 0
   for dist in distributions:
      ax = sns.distplot(dist,bins=bins, kde=False, rug=False, 
                             hist_kws={"histtype":"stepfilled",
                                       "lw":5,
                                       "alpha":0.8,
                                       "color":color_picker(index)})
      index = index + 1
 
   ax.set_title(title,fontsize=16,y=1.08)
   ax.legend(legend,loc=0)
   ax.set_xlabel(xlab)
   ax.set_ylabel(ylab)
   ax.set_xlim(xlim)
   plt.savefig(savefile)
   plt.show()
   
#------------------------------- Structure difference -------------------------------#

def quiver3d(coords,uvw):
    """" Modified form of 3D matplotlib utility quiver, with scaled and colored arrows, 
         according to arrow magnitude """
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    # extract arguments for plotting
    x = coords[:,0]
    y = coords[:,1]
    z = coords[:,2]
    u = uvw[:,0]
    v = uvw[:,1]
    w = uvw[:,2]
  
    # set axes limits
    xmin = np.min(x)
    xmax = np.max(x)
    ymin = np.min(y)
    ymax = np.max(y)
    zmin = np.min(z)
    zmax = np.max(z)
    

    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)
    ax.set_zlim(zmin,zmax)
    ax.quiver(x,y,z,u,v,w)
    plt.show()
