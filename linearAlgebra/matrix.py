#represent matrices in Python 
from PIL import Image  #This imports the Python Image Library (PIL).
import math   #Our usual module for extra math functions
import numpy as np   #This is the main module that we will be using today.  
#It allows us to create and manipulate matrices

get_ipython().magic('matplotlib inline')   #this is included to help make our plots appear right under the cell
import matplotlib.pyplot as plt            #this is the main plotting module

#declaring matrices with numpy
A_1 = np.array([[1,2],[3,4]])  #this is a 2x2 matrix.
#the first brackets (i.e. [1,2]) constitute the first row of the matrix, and similarly for the 
#second brackets ([3,4] second row)
#Be sure to have brackets around the whole thing ... [[1,2],[3,4]] in this example
print("A_1 = \n",A_1)
#you can also make column vectors (i.e. an nx1 matrix)
v_1 = np.array([[10],[7]])  #a 2x1 column matrix
print("v_1 = \n", v_1)

#To get the shape of a matrix (i.e. the number of rows and columns ... (r,c)) use .shape :
print("The shape of A_1 is ", A_1.shape)

#we can also make a few simple matrices without specifying each matrix entry:
Zeros = np.zeros( (2,3) )  #this create a 2x3 matrix whose entries are all zero
print("Zeros = \n", Zeros)
Ones = np.ones(v_1.shape)  #this creates a matix of all ones with the same shape as v_1
print("Ones = \n", Ones)
ID = np.eye(4)   #this creates an identity matrix of size 4x4 (it must be square)
print("Identity Matrix = \n", ID)
