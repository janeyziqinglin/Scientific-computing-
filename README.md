# Scientific-computing 
👉 Linear Algebra   
This is an excellent tool if you want to brush up your linear algebra and learn some basic python programming.
The folder contains guided jupyter Notebooks, using libraries such as Numpy, Scipy and other scientific computing 
library.

👉 Radial distribution function (RDF)\
Radial Distribution Function (RDF) is the probability of finding a pair of atoms a
distance r apart relative to the probability for a completely uniform distribution. It is most commonly used in gasses,
liquids, and solutions, since it can be used to calculate thermodynamic properties such as the internal energy 
and pressure of the system. But is relevant at any size scale, such as packing of colloids, and is useful in complex
heterogeneous media, such as the distribution of ions around DNA. 

The following code shows an example of obtaining RDF for snapshots (128) extracted from a molecular dynamics simulation of water, reported in directory `pbe400_128`. 

The code is composed of four sections:
   - (1) INTRO
   - (2) GET LIST OF FILE NAMES
   			The atomic structure from the XYZ file in directory `pbe400_128` are loaded
   - (3) COMPUTE RDF 
   			The total average, minimum, and maximum RDF for OO, OH, and HH are collected.
   - (4) PLOT 
   			The results are plotted.
   			![rdf_parallel](https://user-images.githubusercontent.com/105125897/173940530-94b4c123-62a4-4c6a-b826-bc281a5e24da.png)

*To speed up the calculation, the code is parallelized and the sample plot is obtained by running the code on four cores.  
mpirun -np 4 python3 rdf.py

Reference: 
J. P. Hansen and I. R. McDonald, Theory of Simple Liquids,  2nd  Ed. (Academic Press, New York, 1986); D. A. McQuarrie, Statistical Mechanics. (Harper & Row, New York, 1976).

