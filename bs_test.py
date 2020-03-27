#!/usr/bin/env python3
"""
-----------------------------------------------------------
 created: wwwennie
   
 file for testing contents of bs.py
----------------------------------------------
"""

import bs

kptfile = "./kpath.dat" # contains high-symm points
bandfile = "./bands.daniel"

allkpts,kptpath, bands = bs.read_bands(bandfile) 
kpts, tics,label = bs.kpt_tics(kptfile,allkpts,kptpath)

bs.plot_bs(kptpath,bands,xtics=tics,labels=label,yrange=[0,15],vbi=103)

