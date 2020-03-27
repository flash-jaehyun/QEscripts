#!/usr/bin/env python3
"""
-----------------------------------------------------------
 created: wwwennie
   
 plot the band structure;
 requires input from multiple files, including:
    bands file
    kpoints file with labelling of high-symmetry points


# example usage:
import bs

folder="/path/to/files/"
bsfile=folder+"bands.ex"
kptfile=folder+"kpath.dat"

allkpts,kptpath, bands = bs.read_bands(bsfile) 
kpts,tics,label = bs.kpt_tics(kptfile,allkpts,kptpath)
bs.plot_bs(allkpts,kptpath,bands,xtics=tics,labels=label,yrange=[0,15],vbi=103)
----------------------------------------------
"""

from matplotlib.ticker import MultipleLocator
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import csv

#-----------------------------  Helper functions -------------------------------#
def cryst2cart(vec,trmat,iflag):

  """
    convert from crystalline to cartesian coordinates and back

    vec = vector
    trmat = transformation matrix
    if iflag = 1, crystallographic to cartesian
        trmat = real-space lattice basis for atoms
              = reciprocal-space lattice basis for k-point
    if iflag = -1: the opposite 
  """
  
  vau = np.zeros(3)

  # iflag cases got flipped for single vector (compared to original in QE)
  # likely since QE version does implicit transpose for the case of many vectors
  #   checked explicitly, but didn't think too much about why 
  if iflag == -1:
     for kpol in range(3):
        vau[kpol] = trmat[kpol][0] * vec[0] + trmat[kpol][1]*vec[1] + trmat[kpol][2]*vec[2]
  elif iflag == 1:
     for kpol in range(3):
        vau[kpol] = trmat[0][kpol] * vec[0] + trmat[1][kpol]*vec[1] + trmat[2][kpol]*vec[2]
  else:
     print("I need a valid iflag option!")

  return vau

def calc_kptpath(rawkpt):
 """
   From a list of k-poits along high-symmetry path, calculate the distance from initial point
   For use in bandstructure plotting

   rawkpt: list of kpts; should use cartesian coordinates
   dist: list of distances; same units as rawkpt
  
 """
 kpts = rawkpt
 
 nkpt,nc = np.shape(kpts)
 dist = np.zeros(nkpt)
 
 for i in range(nkpt):
     if i == 0:
         old = np.array(kpts[i]) 
         new = np.array(kpts[i])
     else:
         old = np.array(kpts[i-1]) 
         new = np.array(kpts[i])
    
     dist[i] = np.linalg.norm(new-old) + dist[i-1]    
 return dist

def find_tics(allkpt,kpts):
 """
   Find the positions for where to put tics lable of high-symmetry points

   allkpt: list of kpt coords corresponding to kpath
   kpts: array; list of high-symmetry kpts

   Basically finding index of high-symmetry k-points in kpts and returning the associated 
   point in kpath

   within tolerance
   accounting for repeats
 """
 idx = []
 allkpt_copy = np.copy(allkpt)
 nrow, ncol = np.shape(allkpt) 

 for k in kpts:
    tmp_k = np.tile(np.array(k), (nrow,1)) 
    idx.append(np.where((abs(np.subtract(allkpt_copy, k)) < 1e-6).all(axis=1))[0][0])
  
    # remove first instance, in case of repeats
    allkpt_copy = np.delete(allkpt_copy,idx[-1],0)
   
 print("Searching for kpath ticks")
 print("Found at indices: ", np.add(idx,1))

 return idx

#-----------------------------  Input -------------------------------#

def read_bands_R(filename, R):
  """ Procedure to read in bands
      Assumes file format from qe-eigenvalues.pl script

      filename (str): file containing bands
      R: lattice vector matrix

     Output:
      bands (list): list of lists containing bands
      kptpath (list): list containing kptpath distances
  """

  rawdat = np.loadtxt(filename,skiprows=1)
  bands = rawdat[:,10:-1]
  kpts = rawdat[:,1:4] # crystalline coords

  # get kptpath as distances; in units of lattice vectors
  G = 2 * np.pi * np.linalg.inv(R).T
  kpts_cart = []
  for k in kpts:
     cart  = cryst2cart(k,G,1)
     kpts_cart.append(cart)
  kptpath = calc_kptpath(kpts_cart)

  # transpose bands in order to iterate through list as columns, not rows of original data
  return kpts, kptpath, np.transpose(bands)

def read_bands(filename):
  """ Procedure to read in bands
      Assumes file format from qe-eigenvalues.pl script

      filename (str): file containing bands

     Output:
      bands (list): list of lists containing bands
      kptpath (list): list containing kptpath distances
  """

  rawdat = np.loadtxt(filename,skiprows=1)
  bands = rawdat[:,10:-1]
  kptpath= rawdat[:,7] # 8th column, crystal units
  #kptpath= rawdat[:,8] # 9th column, 1/Ang units
  kpts = rawdat[:,1:4] # crystalline coords

  # transpose bands in order to iterate through list as columns, not rows of original data
  return kpts, kptpath, np.transpose(bands)

def kpt_tics(kptfile,allkpt,kptpath):
  """ Procedure to extract the k-point path ticks for labelling the x axis
      
      kptfile (str): file containing high-symmetry k-points and labels in format
             
           0.000 0.000 0.000 ! $\Gamma$
           ...

      allkpt: list of crystal coord k-points
      kptpath (list, float): kpt dist of allkpt, plottable x-axis

      assumes label is always final element

     Output: 
      kpts (list): kpts in units of kptpath file
      tics (list): distance of each kpt from initial one, marking ticks
      label (list): and the corresponding high-symmetry labels 
  """ 
  kptx = []
  kpty = []
  kptz = []
  label = []

  with open(kptfile,'r') as infile:
     csv_reader = csv.reader(infile,delimiter=' ',skipinitialspace=True)
     for line in csv_reader:
         kptx.append(float(line[0]))
         kpty.append(float(line[1]))
         kptz.append(float(line[2]))
         label.append(line[-1])

  infile.close()
  kpts = np.stack((kptx,kpty,kptz),axis=1)
  
  idx=find_tics(allkpt,kpts)
  tics = kptpath[idx]

  return kpts,tics,label

def find_gap(bands,vb,cb):
   """ Find indices correpsonding to optical gap given vb and cb index """
 
   ivbm = np.argmax(bands[vb,:])
   icbm = np.argmin(bands[cb,:])

   return ivbm, icbm
   
#-----------------------------  Plotting -------------------------------#

def plot_bs(allkpts,kptpath,bands,xtics=[],labels=[],yrange=[-50,50],eshift=0.0,vbi=1e4):
  """ Plot the bandstructure with some formatting
    
      allkpts (list): list containing all kpt coord
      kptpath (list): list containing kptpath distances
      bands (list): list of lists containing bands
      xtics (list): containing tics and high-symmetry k-point labels
      eshift (float): amount to shift bands by, e.g., E_fermi
      yrange (list): contains [ymin ymax] to plot between
      vbi (integer): band index of topmost valence band
  """
  sns.set_style("white")

  fig = plt.figure()
  ax = fig.gca()

  # plot bands
  bndnum = 0
  for band in bands:
     if bndnum <= vbi:
       ax.plot(kptpath,band,'navy')
     else:
       ax.plot(kptpath,band,'orange')
     bndnum = bndnum + 1

  # add vertical lines
  for i in np.arange(len(xtics)):
     ax.axvline(x=xtics[i],color='black',linewidth=0.5)
 
  # report band gap
  ivbm,icbm = find_gap(bands,vbi,vbi+1)
  cbm = bands[vbi+1,icbm]
  vbm = bands[vbi,ivbm]
  print("Band gap: %.3f"% (cbm-vbm))
  print("CBM: %.3f,%d @ kpt coord: "%(cbm,icbm),allkpts[icbm])
  print("VBM: %.3f,%d @ kpt coord: "%(vbm,ivbm),allkpts[ivbm])

  # format plot
  ax.set_xticks(xtics)
  ax.set_xticklabels(labels,fontsize=20)
  ax.tick_params(labelsize=20)
  ax.set_xlim(0-0.001,max(xtics)+0.001)
  ax.set_ylim(yrange)

  ml = MultipleLocator(10)
  ax.yaxis.set_minor_locator(ml)
  ax.set_ylabel("Energy (eV)",fontsize=20)
  plt.show()
  
