##################################################################
# project: Finite element with c++ 
# Code: input.py
# Developer: Reza Rahimi E-mail: Reza.Rahimi@Dal.Ca
# Date: 06-Jul-2018
# rev: 1.00
##################################################################
# This helps the user to buid finite element mode the example is a cantilever beam with a point load

from finite import Finite
#define variables for this input file
dimentionX = 10.0; numberOfElementX = 2;
dimentionY =  1.0; numberOfElementY =  1;
incr_x = dimentionX/numberOfElementX; # increment between each node 
incr_y = dimentionY/numberOfElementY;
# define the simulation number
cantilever = Finite(1)
# define Nodes
#  6 --- 7 --- 8
#  3 --- 4 --- 5
#  0 --- 1 --- 2
for j in range(numberOfElementY+1):
    for i in range(numberOfElementX+1):
        cantilever.Node(i*incr_x, j*incr_y)  # Node(x,y)
# define fix
for i in range(numberOfElementY+1):
    cantilever.Fix(i*(numberOfElementX+1),1,1) # fix end of the cantliver
# define Load
cantilever.pointLoad((numberOfElementX+1)*(numberOfElementY+1)-1,0.0,-10.0)
# define mesh
'''
    6 --- 7 --- 8
    |  3  |  4  |
    3 --- 4 --- 5
    |  1  |  2  |
    0 --- 1 --- 2
'''
xElement = numberOfElementX+1
for i in range(numberOfElementY):
    for j in range(numberOfElementX):
        cantilever.meshQuadrilateral(i*xElement+j,i*xElement+1+j,(i+1)*xElement+1+j,(i+1)*xElement+j);
# build the model
cantilever.build()
