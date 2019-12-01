# QEscripts
Helper scripts for working with [Quantum Espresso](https://www.quantum-espresso.org/)
Tested with QE v6.+

Dependencies: [ase](https://wiki.fysik.dtu.dk/ase/index.html]), [pymatgen](https://pymatgen.org/)

Scripts for easing intermediate and plotting workflow involving QE:
 - struct: scripts for manipulating and analyzing QE-input structures
 - importing QE input and output files for structure manipulation in ase and pymatgen (dependences: [Xcrysden](http://www.xcrysden.org/)
    1) statistics on structure
        histogram of bond length, bond angle, polyhedra statistics, e.g., 
        
        ![](./example-bond-histrogram.png)
        
    2) introducing point defects and impuritie
   
