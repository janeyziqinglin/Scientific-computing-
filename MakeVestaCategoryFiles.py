#!/usr/bin/env python
# making vesta category files
"""
"""

__copyright__ = "GPL"
__license__ = "Python"

import os
import time
import sys
import re
import numpy
from numpy import *
import operator

def periodic(arg1,arg2):
	#returns periodic difference
      dx=arg2-arg1
      while (dx<-0.5):
                 dx=dx+1.0
      while (dx>0.5):
                 dx=dx-1.0
      return dx 


nOctants=8
nFirstOctant=12
ndim=nFirstOctant*nOctants
nTypes=8
Protons=numpy.zeros((nTypes,ndim),dtype=numpy.int)

indx=0
#reading a file to populate arrays
for line in open("ProtonTypes.dat"):
        dummy=line.split()
        for i in range(0,nFirstOctant):
                Protons[indx][i]=eval(dummy[i]);
        indx=indx+1

if (indx!=nTypes):
        print "indx is != nTypes\n"
        exit()

#Fill out the rest of the array
for indx in range(0,nTypes):
  for i in range(0,nFirstOctant):
      for j in range(1,nOctants):
              Protons[indx][i+nFirstOctant*j]=Protons[indx][i]+96*j

colors=["yellow","green","blue","orange","teal","brown","pink","gray"]
print "TwiceRemovedO_OHpointsCloseO\n",colors[0],Protons[0]
print "TwiceRemovedO_OHpointsFarO\n",colors[1],Protons[1]
print "AdjacentO_OHpointsCloseO\n",colors[2],Protons[2]
print "AdjacentO_OHpointsFarO\n",colors[3],Protons[3]
print "OnceRemovedO_DopantFace_OHpointsCloseO\n",colors[4],Protons[4]
print "OnceRemovedO_DopantFace_OHpointsFarO\n",colors[5],Protons[5]
print "OnceRemovedO_NoDopantFace_OHpointsCloseO\n",colors[6],Protons[6]
print "OnceRemovedO_NoDopantFace_OHpointsFarO\n",colors[7],Protons[7]

lattice=numpy.zeros((3,3))
lengths=numpy.zeros(3)
angles=numpy.zeros(3)

energy=[]
convergence=[]
atomLabel=["Ba","Zr","O","Sc","H"]
atomDescription=["  2.2400 190 101 116  30 239  44 204  0","  1.6000   0 255   0   0 255   0 204  0","  0.7400 254   3   0 254   3   0 204  0","  1.6400 181  99 171 181  99 171 204  0","  0.4600 "]
hCategory=["255 255 128 255 204 204 204  0","128 255 128 255 204 204 204  0","  0 128 255 255 204 204 204  0","255 128  64 255 204 204 204  0"," 64 128 128 255 204 204 204  0","128  64  64 255 204 204 204  0","255 128 255 255 204 204 204  0","192 192 192 255 204 204 204  0"]

nameoffile="AllProtonsAndBackBone.vasp"

numberOfAtomsInContcar=1088
contcarCoord=numpy.zeros((numberOfAtomsInContcar+1,3))
numberOfAtomTypes=5
natoms=numpy.zeros(numberOfAtomTypes,dtype=numpy.int)
j=0
for line in open(nameoffile):
	if j==1:
		dummy=line.split()
		latFactor=eval(dummy[0])
	if j>1 and j <5:
		dummy=line.split()
		lattice[j-2][0]=latFactor*eval(dummy[0])
		lattice[j-2][1]=latFactor*eval(dummy[1])
		lattice[j-2][2]=latFactor*eval(dummy[2])
		lengths[j-2]=numpy.linalg.norm(lattice[j-2])
	if j==5:
		angles[0]=numpy.arccos(numpy.dot(lattice[0],lattice[1])/lengths[0]/lengths[1])*180/numpy.pi
		angles[1]=numpy.arccos(numpy.dot(lattice[1],lattice[2])/lengths[2]/lengths[1])*180/numpy.pi
		angles[2]=numpy.arccos(numpy.dot(lattice[0],lattice[2])/lengths[0]/lengths[2])*180/numpy.pi
		dummy=line.split()
	if j==6:
		dummy=line.split()
		for nn in range(0,numberOfAtomTypes):
		  natoms[nn]=eval(dummy[nn])	
	if (j>7 and j<=7+numberOfAtomsInContcar):  
	    dummy=line.split()
            contcarCoord[j-7][0]=eval(dummy[0]) 	
            contcarCoord[j-7][1]=eval(dummy[1]) 	
            contcarCoord[j-7][2]=eval(dummy[2]) 	
   	    #print j-7,contcarCoord[j-7][0],contcarCoord[j-7][1],contcarCoord[j-7][2]
	j=j+1
nameoffile="AllColors.vesta"
vesta=open(nameoffile,'w')
vesta.write("#VESTA_FORMAT_VERSION 3.3.0\nCRYSTAL\nTITLE\nAllProtonsColors\nGROUP\n1 1 P 1\nSYMOP\n")
vesta.write(" 0.000000  0.000000  0.000000  1  0  0   0  1  0   0  0  1   1\n")
vesta.write(" -1.0 -1.0 -1.0  0 0 0  0 0 0  0 0 0\nTRANM 0\n")
vesta.write(" 0.000000  0.000000  0.000000  1  0  0   0  1  0   0  0  1\nLTRANSL\n -1\n")
vesta.write(" 0.000000  0.000000  0.000000  0.000000  0.000000  0.000000\n")
vesta.write("LORIENT\n -1   0   0   0   0\n 1.000000  0.000000  0.000000  1.000000  0.000000  0.000000\n")
vesta.write(" 0.000000  0.000000  1.000000  0.000000  0.000000  1.000000\nLMATRIX\n")
vesta.write(" 1.000000  0.000000  0.000000  0.000000\n 0.000000  1.000000  0.000000  0.000000\n")
vesta.write(" 0.000000  0.000000  1.000000  0.000000\n 0.000000  0.000000  0.000000  1.000000\n")
vesta.write(" 0.000000  0.000000  0.000000\n")
vesta.write("CELLP\n  "+str(lengths[0])+" "+str(lengths[1])+" "+str(lengths[2])+" "+str(angles[1])+" "+str(angles[2])+" "+str(angles[0])+"\n")
vesta.write("  0.000000   0.000000   0.000000   0.000000   0.000000   0.000000\nSTRUC\n")
atomStart=1 
for itype in range(0,numberOfAtomTypes):
	atomEnd=atomStart+natoms[itype]
  	for nn in range(atomStart,atomEnd):
     		vesta.write("  "+str(nn)+"  "+atomLabel[itype]+" "+atomLabel[itype]+str(nn-atomStart+1)+"  1.0000   "+str(contcarCoord[nn][0])+" "+str(contcarCoord[nn][1])+" "+str(contcarCoord[nn][2])+"    1a       1\n")
		vesta.write("                            0.000000   0.000000   0.000000  0.00\n")
	atomStart=atomEnd
vesta.write("  0 0 0 0 0 0 0\nTHERI 0\n")
atomStart=1
for itype in range(0,numberOfAtomTypes):
	atomEnd=atomStart+natoms[itype]
        for nn in range(atomStart,atomEnd):
		vesta.write(" "+str(nn)+"  "+atomLabel[itype]+str(nn-atomStart+1)+"  1.000000\n")
	atomStart=atomEnd
vesta.write("  0 0 0\nSHAPE\n  0       0       0       0   0.000000  0   192   192   192   192\n")
vesta.write("BOUND\n    -0.1      1.1      -0.1      1.1      -0.1      1.1\n  0   0   0   0  0\n")
vesta.write("SBOND\n  1     O     H    0.00000    1.20000  0  1  0  0  1  0.250  2.000 127 127 127\n")
vesta.write("  2     Sc     O     0.00000    2.42029  0  1  1  0  1  0.250  2.000 127 127 127\n")
vesta.write("  3    Zr     O    0.00000    3.09651  0  1  1  0  1  0.250  2.000 127 127 127\n")
vesta.write("  0 0 0 0\nSITET\n")
atomStart=1
for itype in range(0,numberOfAtomTypes):
        atomEnd=atomStart+natoms[itype]
        for nn in range(atomStart,atomEnd):
		if itype==numberOfAtomTypes-1:
		  #check which atom type
		  for ktypefind in range(0,nTypes): 
			for kk in range(0,ndim):
				if nn-atomStart+1==Protons[ktypefind][kk]:
					iH=ktypefind
                  vesta.write(" "+str(nn)+"  "+atomLabel[itype]+str(nn-atomStart+1)+atomDescription[itype]+hCategory[iH]+"\n")
		else:
                  vesta.write(" "+str(nn)+"  "+atomLabel[itype]+str(nn-atomStart+1)+atomDescription[itype]+"\n")
	print atomStart,natoms[itype],atomEnd
        atomStart=atomEnd

vesta.write("  0 0 0 0 0 0\nVECTR\n 0 0 0 0 0\nVECTT\n 0 0 0 0 0\nSPLAN\n  0   0   0   0\n")
vesta.write("LBLAT\n -1\nLBLSP\n -1\n")
vesta.write("DLATM\n0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63  -1\nDLBND\n -1\nDLPLY\n -1\n")
vesta.write("PLN2D\n  0   0   0   0\nATOMT\n  1         Ba  2.2400 190 101 116  30 239  44 204\n")
vesta.write("  2         Zr  1.6000   0 255   0   0 255   0 204\n  3          O  0.7400 254   3   0 254   3   0 204\n")
vesta.write("  4          Sc 1.6400 181  99 171 181  99 171 204  \n  5          H  0.4600 255 204 204 255 204 204 204\n  0 0 0 0 0 0\n")
vesta.write("SCENE\n-0.990443 -0.074351  0.116164  0.000000\n-0.070893  0.996917  0.033633  0.000000\n")
vesta.write("-0.118306  0.025075 -0.992660  0.000000\n 0.000000  0.000000  0.000000  1.000000\n")
vesta.write("  0.000   0.000\n  0.000\n  1.425\nHBOND 0 2\n\n")
vesta.write("STYLE\nDISPF 39850946\nMODEL   0  1  0\nSURFS   0  1  1\nSECTS  96  1\n")
vesta.write("FORMS   0  1\nATOMS   0  0  1\nBONDS   1\nPOLYS   1\nVECTS 1.000000\n")
vesta.write("FORMP\n  1  1.0   0   0   0\nATOMP\n 24  24   0  50  2.0   0\nBONDP\n")
vesta.write("  1  16  0.250  2.000 127 127 127\nPOLYP\n 204 1  1.000 180 180 180\n")
vesta.write("ISURF\n  0   0   0   0\nTEX3P\n  1 0.00000E+000 1.00000E+000\nSECTP\n")
vesta.write("  1 0.00000E+000 1.00000E+000 0.00000E+000\nHKLPP\n 192 1  1.000 255   0 255\n")
vesta.write("UCOLP\n   0   1  1.000   0   0   0\nCOMPS 1\nLABEL 1    12  1.000 0\nPROJT 1  0.962\n")
vesta.write("BKGRC\n 255 255 255\nDPTHQ 1 -0.5000  3.5000\nLIGHT0 1\n")
vesta.write(" 1.000000  0.000000  0.000000  0.000000\n 0.000000  1.000000  0.000000  0.000000\n")
vesta.write(" 0.000000  0.000000  1.000000  0.000000\n 0.000000  0.000000  0.000000  1.000000\n")
vesta.write(" 0.000000  0.000000 20.000000  0.000000\n 0.000000  0.000000 -1.000000\n  26  26  26 255\n")
vesta.write(" 179 179 179 255\n 255 255 255 255\nLIGHT1\n 1.000000  0.000000  0.000000  0.000000\n")
vesta.write(" 0.000000  1.000000  0.000000  0.000000\n 0.000000  0.000000  1.000000  0.000000\n")
vesta.write(" 0.000000  0.000000  0.000000  1.000000\n 0.000000  0.000000 20.000000  0.000000\n")
vesta.write(" 0.000000  0.000000 -1.000000\n   0   0   0   0\n   0   0   0   0\n   0   0   0   0\n")
vesta.write("LIGHT2\n 1.000000  0.000000  0.000000  0.000000\n 0.000000  1.000000  0.000000  0.000000\n")
vesta.write(" 0.000000  0.000000  1.000000  0.000000\n 0.000000  0.000000  0.000000  1.000000\n")
vesta.write(" 0.000000  0.000000 20.000000  0.000000\n 0.000000  0.000000 -1.000000\n")
vesta.write("   0   0   0   0\n   0   0   0   0\n   0   0   0   0\n")
vesta.write("LIGHT3\n 1.000000  0.000000  0.000000  0.000000\n 0.000000  1.000000  0.000000  0.000000\n")
vesta.write(" 0.000000  0.000000  1.000000  0.000000\n 0.000000  0.000000  0.000000  1.000000\n")
vesta.write(" 0.000000  0.000000 20.000000  0.000000\n 0.000000  0.000000 -1.000000\n")
vesta.write("   0   0   0   0\n   0   0   0   0\n   0   0   0   0\n")
vesta.write("ATOMM\n 204 204 204 255\n  25.600\nBONDM\n 255 255 255 255\n 128.000\n")
vesta.write("POLYM\n  255 255 255 255\n 128.000\nSURFM\n   0   0   0 255\n 128.000\n")
vesta.write("FORMM\n 255 255 255 255\n 128.000\nHKLPM\n 255 255 255 255\n 128.000\n")
vesta.close()
  

