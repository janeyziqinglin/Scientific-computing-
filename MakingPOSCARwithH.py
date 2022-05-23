#!/usr/bin/env python
"""
"""

__copyright__ = "GPL"
__license__ = "Python"

import os
import time
import sys
import re
import numpy
from math import *
import operator

def periodic(arg1,arg2):
        #returns periodic difference
      dx=arg2-arg1
      while (dx<-0.5):
                 dx=dx+1.0
      while (dx>0.5):
                 dx=dx-1.0
      return dx;

lattice=numpy.zeros((3,3))
initialBackbonePositions=numpy.zeros((320,3))

counter=0
counterBackbone=0

for line in open("POSCARAdj_filled.vasp"):
	if counter==1:
		dummy=line.split()
		size=eval(dummy[0])
        if counter>1 and counter <5:
                dummy=line.split()
                lattice[counter-2][0]=eval(dummy[0])
                lattice[counter-2][1]=eval(dummy[1])
                lattice[counter-2][2]=eval(dummy[2])
        if counter>7:
                dummy=line.split()
                initialBackbonePositions[counterBackbone][0]=eval(dummy[0])
                initialBackbonePositions[counterBackbone][1]=eval(dummy[1])
                initialBackbonePositions[counterBackbone][2]=eval(dummy[2])
                counterBackbone=counterBackbone+1
        counter=counter+1

print "counter is "+str(counter)
print "counterBackbone is "+str(counterBackbone)

#Checking that we can print out what we read in
#print "Y doped BaZrO3 FN VO Glazer 11"                      
#print size 
#for k in range(0,3):
#  print lattice[k][0],lattice[k][1],lattice[k][2]
#print " Ba   Zr   O   Sc"
#print "  64    56    191   8"
#print "Direct"
#for k in range(0,counterBackbone):
#  print initialBackbonePositions[k][0],initialBackbonePositions[k][1],initialBackbonePositions[k][2]

nH=0
faceLabels=[]
HPos=numpy.zeros((769,3))
# reading in all possible protons and labels
#for line in open("HlabelXYZ.dat"):
for line in open("HlabelXYZStraighter.dat"):
	dummy=line.split()
	faceLabels.append(dummy[1])
	nH=nH+1
	# first set in positions 1-96
	# divided by 2 since we'll be in half the size
	HPos[nH][0]=eval(dummy[2])/2.0
	HPos[nH][1]=eval(dummy[3])/2.0
	HPos[nH][2]=eval(dummy[4])/2.0
	# shifted in the x direction by 0.5 (97-192)
	HPos[96+nH][0]=HPos[nH][0]+0.5
	HPos[96+nH][1]=HPos[nH][1]
	HPos[96+nH][2]=HPos[nH][2]
	#shifted in the y direction by 0.5 (193-288)
	HPos[192+nH][0]=HPos[nH][0]
	HPos[192+nH][1]=HPos[nH][1]+0.5
	HPos[192+nH][2]=HPos[nH][2]
	# shiftedin x and y direction by 0.5 in each (289-384)
	HPos[288+nH][0]=HPos[nH][0]+0.5
	HPos[288+nH][1]=HPos[nH][1]+0.5
	HPos[288+nH][2]=HPos[nH][2]
	# shifted in z direction by 0.5 (385-480) 
	HPos[384+nH][0]=HPos[nH][0]
	HPos[384+nH][1]=HPos[nH][1]
	HPos[384+nH][2]=HPos[nH][2]+0.5
	# shifted in the x and z direction by 0.5 (481-576)
	HPos[480+nH][0]=HPos[nH][0]+0.5
	HPos[480+nH][1]=HPos[nH][1]
	HPos[480+nH][2]=HPos[nH][2]+0.5
	#shifted in the y and z direction by 0.5 (577-672)
	HPos[576+nH][0]=HPos[nH][0]
	HPos[576+nH][1]=HPos[nH][1]+0.5
	HPos[576+nH][2]=HPos[nH][2]+0.5
	# shiftedin x and y and z direction by 0.5 in each (673-768)
	HPos[672+nH][0]=HPos[nH][0]+0.5
	HPos[672+nH][1]=HPos[nH][1]+0.5
	HPos[672+nH][2]=HPos[nH][2]+0.5
	#print "H "+dummy[1],nH,HPos[nH][0],HPos[nH][1],HPos[nH][2]

#Finding the protons closest to each of the oxygens
closestHtoeachO=numpy.zeros((192,4),dtype=numpy.int) 
nH_222=96*8
ODistToH=numpy.zeros(nH_222+1)
sortLabel=numpy.zeros(nH_222+1,dtype=numpy.int)
sortDisplacement=numpy.zeros(3)
whichDisplacement=numpy.zeros(3,dtype=numpy.int)
for n in range(120,312): #zero index in array! first O is at 120 and they end at 311 
	for m in range(1,nH_222+1):
		#print n-120+1,initialBackbonePositions[n][0],initialBackbonePositions[n][1],initialBackbonePositions[n][2]
		#print m,HPos[m][0],HPos[m][1],HPos[m][2]
		dx=periodic(initialBackbonePositions[n][0],HPos[m][0])
		dy=periodic(initialBackbonePositions[n][1],HPos[m][1])
		dz=periodic(initialBackbonePositions[n][2],HPos[m][2])
		dxC=lattice[0][0] * dx + lattice[1][0] * dy + lattice[2][0] * dz
		dyC=lattice[0][1] * dx + lattice[1][1] * dy + lattice[2][1] * dz
		dzC=lattice[0][2] * dz + lattice[1][2] * dy + lattice[2][2] * dz
		ODistToH[m]=size*((dxC ** (2)) + (dyC ** (2)) + (dzC ** (2))) ** (0.5)
		sortLabel[m]=m
		#print n-120+1,sortLabel[m],faceLabels[m-1],ODistToH[m]
	#sorting
	# bubble sort from book "Discrete Mathematics"
	#end of list is nH_222 and the range function goes one less than the end of the list
	for j in range(1, nH_222):
		#end of list if nH_222
		# one j less than the end of list is nH_222-j but range goes one less so end of 
		#range had to be nH_222-j+1
		for m in range(1, nH_222 - j+1):
			if ODistToH[nH_222 + 1 - m] < ODistToH[nH_222 - m]:
				#putting in earlier larger distance into a temporary variable and
				# switching order
				Temp = ODistToH[nH_222 - m]
				ODistToH[nH_222 - m] = ODistToH[nH_222 + 1 - m]
				ODistToH[nH_222 + 1 - m] = Temp
				# same switching for label
				iTemp = sortLabel[nH_222 - m]
				sortLabel[nH_222 - m] = sortLabel[nH_222 + 1 - m]
				sortLabel[nH_222 + 1 - m] = iTemp
	
	#print "Sorted O to H distance "
	#print "O index O | distance from closestH"
	#for j in range(1, nH_222+1):
	#	print str(j)+" "+str(sortLabel[j]) +" "+ faceLabels[sortLabel[j]-1]+" " + str(ODistToH[j])
	#print
	#Next we note the closest 4 Hs to each O
	for j in range(1,5):
		closestHtoeachO[n-120][j-1]=sortLabel[j]
		#print str(j)+" "+str(sortLabel[j]) +" "+ faceLabels[sortLabel[j]-1]+" " + str(ODistToH[j])
		#Next, we identify the direction of each OH bond
                sortDisplacement[0]=periodic(initialBackbonePositions[n][0],HPos[sortLabel[j]][0])
                sortDisplacement[1]=periodic(initialBackbonePositions[n][1],HPos[sortLabel[j]][1])
                sortDisplacement[2]=periodic(initialBackbonePositions[n][2],HPos[sortLabel[j]][2])
		#print HPos[sortLabel[j]][0],HPos[sortLabel[j]][1],HPos[sortLabel[j]][2]
		#print sortDisplacement[0],sortDisplacement[1],sortDisplacement[2]
		largest=0.0 
		khigh=0
		for kk in range (0,3):
			while (abs(sortDisplacement[kk])>largest):
				largest=abs(sortDisplacement[kk])
				khigh=kk
		#print "The largest displacement in magnitude is ",sortDisplacement[khigh]," and occurs for ",khigh
		# next, we change HPos[sortLabel[j]] to have a distance of 1.0 Angstroms and to be in one main direction	
		for kk in range(0,3):
			HPos[sortLabel[j]][kk]=initialBackbonePositions[n][kk]+numpy.sign(sortDisplacement[khigh])*lattice[khigh][kk]*0.5/(2*size) 
		print "O",n-120+1,"H ",sortLabel[j],(n-120)*4+j,HPos[sortLabel[j]][0],HPos[sortLabel[j]][1],HPos[sortLabel[j]][2]

	#print n-120+1,closestHtoeachO[n-120][0],closestHtoeachO[n-120][1],closestHtoeachO[n-120][2],closestHtoeachO[n-120][3]
	#print

#Setting the vacancy location
whichVO=raw_input("Which vacancy are you considering?(1-192)\n")
nVO=eval(whichVO)
print "The closest Hs to this vacancy are: "
print closestHtoeachO[nVO-1][0],closestHtoeachO[nVO-1][1],closestHtoeachO[nVO-1][2],closestHtoeachO[nVO-1][3]

#Making individual POSCARs and also a file with all protons
poscarAllProtons=open("AllProtonsAndBackBone.vasp","w")
poscarAllProtons.write("All Protons\n")
poscarAllProtons.write(str(size)+"\n")
for k in range(0,3):
     poscarAllProtons.write(str(lattice[k][0])+" "+str(lattice[k][1])+" "+str(lattice[k][2])+"\n")
poscarAllProtons.write (" Ba   Zr   O  Sc  H\n")
poscarAllProtons.write ("  64    56    192   8  768\n")
poscarAllProtons.write("Direct"+"\n")
for k in range(0,counterBackbone):
            poscarAllProtons.write(str(initialBackbonePositions[k][0])+" "+ str(initialBackbonePositions[k][1])+" "+str(initialBackbonePositions[k][2])+"\n")
for i in range(1,nH_222+1):
	    poscarAllProtons.write(str(HPos[i][0])+" "+str(HPos[i][1])+" "+str(HPos[i][2])+"\n")
poscarAllProtons.close()

#setting up directories and files for the specific vacancy making sure not to use the protons attached to the vacancy 
for i in range(1,nH_222+1): 
    if ((i!=closestHtoeachO[nVO-1][0]) and (i!=closestHtoeachO[nVO-1][1]) and (i!=closestHtoeachO[nVO-1][2]) and (i!=closestHtoeachO[nVO-1][3])):
	#writing the POSCAR file without distortion (1)
	nameoffile = "POSCAR"+str(i)+".vasp"
	#print nameoffile
	poscar = open(nameoffile,'w')
	poscar.write("Proton"+str(i)+"\n")                      
	poscar.write(str(size)+"\n")
	for k in range(0,3):
  		poscar.write(str(lattice[k][0])+" "+str(lattice[k][1])+" "+str(lattice[k][2])+"\n")
	poscar.write (" Ba   Zr   O  Sc  H\n")
	poscar.write ("  64    56    191   8  1\n")
	poscar.write("Direct"+"\n")
	for k in range(0,counterBackbone):
		if k != 64+56+nVO-1:
  			poscar.write(str(initialBackbonePositions[k][0])+" "+ str(initialBackbonePositions[k][1])+" "+str(initialBackbonePositions[k][2])+"\n")
	poscar.write(str(HPos[i][0])+" "+str(HPos[i][1])+" "+str(HPos[i][2])+"\n")
	poscar.close()
	#writing the script files
	nameofscript="script"+str(i)+".slurm"
	script = open(nameofscript, 'w')
	script.write("#!/bin/tcsh\n\n")
	script.write("#SBATCH -p stdmem\n")
	script.write("#SBATCH -n 20 #Number of cores\n")
	script.write("#SBATCH -N 1 #number of nodes\n")
	script.write("#SBATCH -t 48:00:00\n")
	script.write("#SBATCH --mem=100GB\n")
	script.write("#SBATCH --export=ALL\n")
	script.write("#SBATCH --job-name=\"V"+str(nVO)+"_P"+str(i)+"\"\n")
	script.write("#SBATCH --mail-user=lin24z@mtholyoke.edu\n")
	script.write("#SBATCH --mail-type=ALL\n")
	script.write("\nsource /usr/share/Modules/init/tcsh\nset echo\n\n")
	script.write("module load intel openmpi-1.6/intel\n")
	script.write("echo $SLURM_JOB_NODELIST\ncd $SLURM_SUBMIT_DIR\nsleep 5\n")
	if i <10 :
		nameDir="V"+str(nVO)+"P0"+str(i)
	else:
		nameDir="V"+str(nVO)+"P"+str(i)
	script.write("mkdir "+nameDir+"\n")
	script.write("cp INCARmin "+" "+nameDir+"/INCAR\n")
	script.write("cp KPOINTS"+" "+nameDir+"\n")
	script.write("cp POTCAR"+" "+nameDir+"\n")
	script.write("cp "+nameoffile+" "+nameDir+"/POSCAR\n")
	script.write("cd "+nameDir+"\nsleep 5\n")	
	script.write("mpiexec /home/gomez/bin/vasp5.3gamma > vasp.out\nsleep 5\n")
	script.write("rm -f CHG CHGCAR EIGENVAL vasprun.xml IBZKPT OSZICAR WAVECAR XDATCAR DOSCAR PCDAT\nsleep5\n")














