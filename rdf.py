import numpy as np
import os
import matplotlib.pyplot as plt
from mpi4py import MPI

############
# 1. INTRO #
############

class XYZfile:
    """This class loads an atomic structure from a XYZ file."""

    def __init__(self, filename):
        """Reads a XYZ file

        filename (str) : The name of the file
        """

        self.filename = filename
        with open(self.filename,'r') as file :
            # number of atoms
            self.np = int(file.readline())
            # comment
            self.comment = file.readline()
            # data
            self.data = np.zeros(shape=(self.np,),dtype=[("element","U9"),("x","f8"),("y","f8"),("z","f8")])
            for i in range(self.np):
                parts = file.readline().strip().split(' ')
                self.data[i] = (parts[0],parts[1],parts[2],parts[3])

class CubicCell:

    def __init__(self, comment):
        """Reads a cubic cell from a comment contained in the xyz file

        comment (str) : The comment
        """
        self.L = np.zeros(shape=(3,),dtype="f8")
        parts = comment.strip().split(' ')
        self.L[0] = parts[1]  # period along x
        self.L[1] = parts[5]  # period along y
        self.L[2] = parts[9]  # period along z

    def volume(self):
        """Computes the volume of the cell"""
        return self.L[0]*self.L[1]*self.L[2]  # volume of a cube

    def wrap(self,v):
        """Wraps the vector v into the unit cell

        v (3-dim ndarray) : vector
        """
        return np.remainder(v, self.L)  # wraps a vector into the cell

    def pbc_distance(self,v1,v2):
        """Computes the distance between two vectors with PBC

        v1 (3-dim ndarray) : vector
        v2 (3-dim ndarray) : vector
        """
        d = np.array(v1,dtype="f8") - np.array(v2,dtype="f8") # general difference
        d = np.remainder(d + self.L/2.0, self.L) - self.L/2.0  # minimum image difference
        return np.linalg.norm(d) # norm of the minimum image difference

def sphere(r):
    return(4/3*np.pi*r**3)

def volume2(r,L):
    x = r/L
    return (-np.pi/12*(3-36*x**2+32*x**3))*L**3

def volume3(r,L):
    x = r/L
    return (-np.pi/4. + 3*np.pi*x**2 + np.sqrt(4*x**2-2) + (1-12*x**2)*np.arctan(np.sqrt(4*x**2-2)) + 2/3*x**2*8*x*np.arctan((2*x*(4*x**2-3))/(np.sqrt(4*x**2-2)*(4*x**2+1))) )*L**3

def gauss1D(x,mu=0,sigma=1):
    return np.exp(-0.5*(x-mu)**2/sigma/sigma)/np.sqrt(2.*np.pi)/sigma

def rdf(f,Rmax,NRmax,element_center,element_distant) :

    # get the cell
    cell = CubicCell(f.comment)

    is_center = np.char.startswith(f.data[:]["element"],element_center)
    is_distant = np.char.startswith(f.data[:]["element"],element_distant)

    # count the number of center particles
    Np_center = np.sum(is_center)

    # count the number of distant particles
    Np_distant = np.sum(is_distant)

    bins = np.linspace(0,Rmax,NRmax+1,endpoint=True)
    radii = np.zeros(shape=(NRmax,),dtype="f8")
    hist = np.zeros(shape=(NRmax,),dtype="f8")
    volumes = np.zeros(shape=(NRmax,),dtype="f8")

    for i in range(NRmax):
        radii[i] = (bins[i]+bins[i+1])/2. # middle point between two bins
        if radii[i]<=cell.L[0]/2 :
            volumes[i]=sphere(bins[i+1])-sphere(bins[i])
        elif radii[i]<=cell.L[0]*np.sqrt(2)/2  :
            volumes[i]=volume2(bins[i+1],cell.L[0])-volume2(bins[i],cell.L[0])
        elif radii[i]<=cell.L[0]*np.sqrt(3)/2 :
            volumes[i]=(volume3(bins[i+1],cell.L[0])-volume3(bins[i],cell.L[0]))
        else :
            volumes[i]=0

    # accumulate gaussians in hist
    for i in range(f.np):
        if is_center[i]: #filter
            Ri = [f.data[i]["x"], f.data[i]["y"], f.data[i]["z"] ]
            for j in range(f.np):
                if is_distant[j] and (i!=j): #filter
                    Rj = [f.data[j]["x"], f.data[j]["y"], f.data[j]["z"] ]
                    d = cell.pbc_distance(Ri,Rj)
                    hist += gauss1D(radii,d,0.2)*Rmax/NRmax

    rho_avg = Np_distant / cell.volume()
    return np.divide(hist, volumes, out=np.zeros_like(hist), where=volumes!=0)/rho_avg/Np_center, radii

#############################
# 2. GET LIST OF FILE NAMES #
#############################

fnames = []

# We loop over all files that end with ".xyz"
for dirpath, dirnames, filenames in os.walk("pbe400_128", topdown=False):
    xyz_files = [f for f in filenames if f.endswith('.xyz')]
    for name in xyz_files:
        fnames.append(os.path.join(dirpath, name))

##################
# 3. COMPUTE RDF #
##################

Rmax = 9
NRmax = 100
no_files = len(fnames)

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()


g_OO_avg = np.zeros(shape=(NRmax,),dtype="f8")
g_OH_avg = np.zeros(shape=(NRmax,),dtype="f8")
g_HH_avg = np.zeros(shape=(NRmax,),dtype="f8")
g_OO_max = np.zeros(shape=(NRmax,),dtype="f8")
g_OH_max = np.zeros(shape=(NRmax,),dtype="f8")
g_HH_max = np.zeros(shape=(NRmax,),dtype="f8")
g_OO_min = np.ones(shape=(NRmax,),dtype="f8")*100
g_OH_min = np.ones(shape=(NRmax,),dtype="f8")*100
g_HH_min = np.ones(shape=(NRmax,),dtype="f8")*100

# LOOP that needs to be parallelized
#divide 128 snapshots to number of nodes


mystates=[]
for v in range(no_files):
    if rank == v % size:
        mystates.append(v)
maxlocstates=comm.allreduce(len(mystates),op=MPI.MAX)
print(rank,maxlocstates,mystates)

for i, fname in enumerate(fnames):
    if i in mystates :
        f = XYZfile(fname)
        g_OO, r_OO = rdf(f,Rmax,NRmax,"O","O")
        g_OH, r_OH = rdf(f,Rmax,NRmax,"O","H")
        g_HH, r_HH = rdf(f,Rmax,NRmax,"H","H")
    # mean:
        g_OO_avg += g_OO / no_files
        g_OH_avg += g_OH / no_files
        g_HH_avg += g_HH / no_files
    #max:
        g_OO_max = np.maximum(g_OO,g_OO_max)
        g_OH_max = np.maximum(g_OH,g_OH_max)
        g_HH_max = np.maximum(g_HH,g_HH_max)
    #min:
        g_OO_min = np.minimum(g_OO,g_OO_min)
        g_OH_min = np.minimum(g_OH,g_OH_min)
        g_HH_min = np.minimum(g_HH,g_HH_min)
        print(f"{i}/{no_files-1}, Processing file: {fname}...")
# mean_buff:
sendbufooav = g_OO_avg
sendbufohav = g_OH_avg
sendbufhhav = g_HH_avg
recvbufooav = np.zeros_like(sendbufooav)
recvbufohav = np.zeros_like(sendbufohav)
recvbufhhav = np.zeros_like(sendbufhhav)
comm.Reduce(sendbufooav, recvbufooav, op=MPI.SUM, root=0)
comm.Reduce(sendbufohav, recvbufohav, op=MPI.SUM, root=0)
comm.Reduce(sendbufhhav, recvbufhhav, op=MPI.SUM, root=0)
g_OO_avg=recvbufooav
g_OH_avg=recvbufohav
g_HH_avg=recvbufhhav

# max_buff:
sendbufoomax = g_OO_max
sendbufohmax = g_OH_max
sendbufhhmax = g_HH_max
recvbufoomax = np.zeros_like(sendbufoomax)
recvbufohmax = np.zeros_like(sendbufohmax)
recvbufhhmax = np.zeros_like(sendbufhhmax)
comm.Reduce(sendbufoomax, recvbufoomax, op=MPI.MAX, root=0)
comm.Reduce(sendbufohmax, recvbufohmax, op=MPI.MAX, root=0)
comm.Reduce(sendbufhhmax, recvbufhhmax, op=MPI.MAX, root=0)
g_OO_max=recvbufoomax
g_OH_max=recvbufohmax
g_HH_max=recvbufhhmax

# min_buff:
sendbufoomin = g_OO_min
sendbufohmin = g_OH_min
sendbufhhmin = g_HH_min
recvbufoomin = np.zeros_like(sendbufoomin)
recvbufohmin = np.zeros_like(sendbufohmin)
recvbufhhmin = np.zeros_like(sendbufhhmin)
comm.Reduce(sendbufoomin, recvbufoomin, op=MPI.MIN, root=0)
comm.Reduce(sendbufohmin, recvbufohmin, op=MPI.MIN, root=0)
comm.Reduce(sendbufhhmin, recvbufhhmin, op=MPI.MIN, root=0)
g_OO_min=recvbufoomin
g_OH_min=recvbufohmin
g_HH_min=recvbufhhmin

###########
# 4. PLOT #
###########
if rank ==0:
    file_name = "rdf_parallel.png"
    plt.plot(r_OO,g_OO_avg,marker="o",color="r",label="g_OO")
    plt.plot(r_OO,g_OO_min,color="r")
    plt.plot(r_OO,g_OO_max,color="r")
    plt.plot(r_OH,g_OH_avg,marker="o",color="b",label="g_OH")
    plt.plot(r_OH,g_OH_min,color="b")
    plt.plot(r_OH,g_OH_max,color="b")
    plt.plot(r_HH,g_HH_avg,marker="o",color="g",label="g_HH")
    plt.plot(r_HH,g_HH_min,color="g")
    plt.plot(r_HH,g_HH_max,color="g")
    plt.xlabel(r"r ($\AA$)")
    plt.xlim(0,Rmax)
    plt.ylim(0,4)
    plt.ylabel(r"g(r)")
    plt.legend()
    plt.savefig(file_name)
    print(f"File saved: {file_name}")
