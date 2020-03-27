#!/usr/bin/env python3
"""
-----------------------------------------------------------
 created: wwwennie
   
 quick and easy way to plot lots of macroscopic and microscopic
 average potentials
  
 as outputted by QE v6+ pp.x 
----------------------------------------------
"""

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np


# constants
rytoev = 13.6056980659
bohrtoang = 0.529177249

#================================================================================#
#------------------------------- Helper functions -------------------------------#
def color_picker(index):
  """ For when more control over colors is preferred 
      Cycles over the possible options """
  colors= ["blue","orange","green","magenta","darkcyan","cornflowerblue","r","indigo","goldenrod","navy"]
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

#------------------------------- Microscopic V_potential -------------------------------#
def micro_pot(potfile,ax=None):
 """ plot microscopic planar average
    
     ax = optional argument of axes object; for ability to plot on same axes between
          different function calls
 """
 sns.set()
 sns.set_style("ticks")

 if ax is None: 
   fig = plt.figure()
   ax = fig.gca()

 data = np.loadtxt(potfile)

 plt.plot(np.multiply(data[:,0],bohrtoang),np.multiply(data[:,1],rytoev))
 plt.xlabel("z (Ang)", fontsize=20)
 plt.ylabel(r"V$_{micro}$ (eV)", fontsize=20)
 
 return ax

#------------------------------- Macroscopic V_potential -------------------------------#
def macro_pot(potfile,ax=None):
 """ plot macroscopic planar average
    
     ax = optional argument of axes object; for ability to plot on same axes between
          different function calls
 """
 sns.set()
 sns.set_style("ticks")
 
 data = np.loadtxt(potfile)
 if ax is None: 
   fig=plt.figure()
   ax=fig.gca()
   plt.plot(np.multiply(data[:,0],bohrtoang),np.multiply(data[:,2],rytoev), linewidth=2)
   plt.xlabel("z (Ang)")
   plt.ylabel(r"V$_{macro}$ (eV)")
 else:
   ax.plot(np.multiply(data[:,0],bohrtoang),np.multiply(data[:,2],rytoev),linewidth=2)
   ax.set_xlabel("z (Ang)")
   ax.set_ylabel(r"V$_{avg}$ (eV)")

 return ax
