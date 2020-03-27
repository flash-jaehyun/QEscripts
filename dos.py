#!/usr/bin/env python3
"""
-----------------------------------------------------------
 created: wwwennie
   
 plot DOS and pDOS, with possibility of nspin
 compatible with QE v6+ outputs of dos.x and projwfc.x


# example usage
folder="/path/to/files/"
pdosfiles=['pdos_Bi_s','pdos_V_d','pdos_O_p']
legend=['Total','Bi s','V d','O p']
pdosfiles=[folder+"/pdos/"+s for s in pdosfiles]
dosfile=folder+"/dos/"+'dos.ex'

dos.pdos_plot(pdosfiles,nspin=2,legend=legend,xlim=[0,17],ylim=[0,25],scale=1.0,dos=True,dosfile=dosfile)
----------------------------------------------
"""

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

#================================================================================#
#------------------------------- Helper functions -------------------------------#
def color_picker(index):
  """ For when more control over colors is preferred 
      Cycles over the possible options """
 # colors= ["blue","orange","green","magenta","darkcyan","cornflowerblue","r","indigo","goldenrod","navy"]
 # colors= ["darkmagenta","fuchsia","darkcyan","orange","firebrick","cornflowerblue","indigo","goldenrod","red","green"]
  colors= ["fuchsia","darkmagenta","darkcyan","orange","green","cornflowerblue","darkred","goldenrod","indigo","firebrick"]
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

# helper function for generated dos above but with more smearing ad-hoc 
# via convolution with Gaussian
def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise(ValueError,"smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise(ValueError, "Input vector needs to be bigger than window size.")


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise(ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")


    s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    return y

#------------------------------- Total DOS -------------------------------#

def dos_plot(filename,nspin=1,xmin=-10,xmax=20,scale=1.0): 
   """
     plot the DOS as outputted by QE v6+ dos.x
  
     if nspin = 1:
     #  E (eV)   dos(E)     Int dos(E) 
     if nspin = 2:
     #  E (eV)   dos(E)up    dos(E)down     Int dos(E) 
  
     xmin,xmax = in eV; range of energies to plot
     scale  = how much to scale the DOS by; e.g., by number of unit cells 
   """
   sns.set()
   sns.set_style("ticks")

   fig,ax = plt.figure()
   #ax = fig.gca()

   # import data
   if (nspin==1):
     key = ['energy','dos','intdos']
   elif (nspin==2):
     key = ['energy','dosup','dosdown','intdos']

   data = pd.read_csv(filename,delim_whitespace=True,names=key,dtype=np.float32)

   # plot
   if (nspin == 1):
      plt.plot(data['energy'],data['dos'],'k')
   
   if (nspin == 2):
      plt.plot(data['energy'],data['dosup'],'k')
      plt.plot(data['energy'],np.multiply(-1.0,data['dosdown']),'k')
      
   ax.legend(loc=0)
   ax.set_xlabel("Energy (eV)")
   ax.set_ylabel("DOS")
   plt.show()
   
def pdos_plot(pdosfiles,legend='',nspin=1,xlim=[0,20],ylim=[0,20],scale=1.0,dos=True,dosfile='',lsmooth=False,shift=0.0): 
   """
     plot the partial DOS as outputted by QE v6+ dos.x
     with also option to plot total DOS
  
     if nspin = 1:
     #  E (eV)   pdos(E)      
     if nspin = 2:
     #  E (eV)   pdos(E)up    pdos(E)down     
    
     pdosfiles = list of filenames to plot
     legend = list of strings presenting name of each curve plotted
     xlim,ylim = in eV and DOS respectively; range to plot
     dos = whether to also plot total dos
     dosfile = correponding DOS file
     scale  = how much to scale the DOS by; e.g., by number of unit cells 
     smooth = boolean, whether to run smoothing procedure
     shift  = rigid shift of energies; for aligning to zero

    TODO:
      legend for nspin = 2 (double labelling happening)
   """
   sns.set()
   sns.set_style("ticks")

   fig = plt.figure()
   ax = fig.gca()

   # plot total DOS first, if desired
   if dos:   
     # import data
     if (nspin==1):
       key = ['energy','dos','intdos']
     elif (nspin==2):
       key = ['energy','dosup','dosdown','intdos']
  
     data = pd.read_csv(dosfile,delim_whitespace=True,header=1,names=key,dtype=np.float32)
     energy = np.add(data['energy'],shift)
     win_len = 11
  
     # plot
     numdat = len(energy)
     if (nspin == 1):
        if lsmooth:
           dos = smooth(data['dos'],window_len=win_len)[:numdat]
        else:
           dos = data['dos']
        plt.plot(energy,dos,'k')
      
     if (nspin == 2):
        if lsmooth:
            dosup = smooth(data['dosup'],window_len=win_len+100)[:numdat]
            dosdown = smooth(data['dosdown'],window_len=win_len+100)[:numdat]
            print("Smoothing Total Dos",win_len+30)
        else:
            dosup = data['dosup']
            dosdown = data['dosdown']
        plt.plot(energy,dosup,'k')
        plt.plot(energy,np.multiply(-1.0,dosdown),'k')

   # import data
   if (nspin==1):
     key = ['energy','pdos']
   elif (nspin==2):
     key = ['energy','pdosup','pdosdown']

   # plot partial density of states
   index=0
   for fil in pdosfiles:
     data = pd.read_csv(fil,delim_whitespace=True,header=1,names=key,dtype=np.float32)
     energy = np.add(data['energy'], shift)
     numdat = len(energy)
     color = color_picker(index)
     if (nspin == 1):
         if lsmooth:
            pdos = smooth(data['pdos'],window_len=win_len)[:numdat]
         else:
            pdos = data['pdos']
         plt.plot(energy,pdos,color)
     
     if (nspin == 2):     
         if lsmooth:
            pdosup = smooth(data['pdosup'],window_len=win_len)[:numdat]
            pdosdown = smooth(data['pdosdown'],window_len=win_len)[:numdat]
         else:
            pdosup = data['pdosup']
            pdosdown = data['pdosdown']
           
         plt.plot(energy,pdosup,color)
         plt.plot(energy,np.multiply(-1.0,pdosdown),color)
     index=index+1
    
 
   ax.legend(legend,loc=0)
   plt.xlim(xlim[0],xlim[1])
   plt.ylim(ylim[0],ylim[1])
   plt.xticks(fontsize=15)
   plt.yticks(fontsize=15)
   plt.xlabel("Energy (eV)", fontsize=20)
   plt.ylabel("DOS",fontsize=20)
   plt.show()

   fig.savefig("/home/wwwennie/Downloads/dos.pdf",bbox_inches="tight")
