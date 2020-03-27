#!/usr/bin/env python3
"""
-----------------------------------------------------------
 created: wwwennie
   
 example file to demonstrate plotting capabilities in bs.py and dos.py
----------------------------------------------
"""

import bs
import dos

folder = "./example/"

kptfile = folder+"kpath.dat" # contains high-symm points
bandfile = folder+"bands_example.dat"

allkpts,kptpath, bands = bs.read_bands(bandfile) 
kpts, tics,label = bs.kpt_tics(kptfile,allkpts,kptpath)

bs.plot_bs(allkpts,kptpath,bands,xtics=tics,labels=label,yrange=[-5,5],vbi=103,eshift=-7.74)

#-------------------------

pdosfiles=[folder+i for i in ["pdos_Bi_s","pdos_O_p","pdos_V_d"]]
legend = ["Total","Bi s","O p", "V d"]
dosfile = "dos_example.dat"
dos.pdos_plot(pdosfiles,legend=legend,nspin=1,xlim=[-5,5],ylim=[0,30],dosfile=folder+dosfile,lsmooth=True,shift=-8.61)
