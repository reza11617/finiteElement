##################################################################
# project: Finite element with c++ 
# Code: finite.py
# Developer: Reza Rahimi E-mail: Reza.Rahimi@Dal.Ca
# Date: 06-Jul-2018
# rev: 1.00
##################################################################
# This program reads the input file and prints a dat file for c program to use as input

class Finite(object):
    def __init__ (self,simNumber):
        self.outFile = self.start(simNumber)
        self.nodeNumber = 0
        self.fixNumber = 0
        self.pointLoadNumber = 0
        self.meshQuadrilateralNumber = 0
        
    def start(self,simNumber):
        fileName = '/home/reza/phdFiles/programing/finiteElementC/src/Cinput/sim'+ str(simNumber)+'.dat'
        return open(fileName,'w')

    def Node(self,x,y): #   Code number ,x,y
        self.outFile.write('100' + str(self.nodeNumber)+','+str(x)+','+str(y)+'\n')
        self.nodeNumber = self.nodeNumber + 1

    def Fix(self,nodeTag,Dof_x,Dof_y): # code number, nodeTag, dof_x, dof_y
        self.outFile.write('200' + str(self.fixNumber)+','+str(nodeTag) +','+str(Dof_x)+','+str(Dof_y)+'\n')
        self.fixNumber = self.fixNumber + 1

    def pointLoad(self,nodeTag,Dof_x,Dof_y): # code number, nodeTag, load dof_x, load dof_y
        self.outFile.write('301' + str(self.pointLoadNumber)+','+str(nodeTag) +','+str(Dof_x)+','+str(Dof_y)+'\n')
        self.pointLoadNumber = self.pointLoadNumber + 1

    def meshQuadrilateral(self, nodeTag1, nodeTag2, nodeTag3, nodeTag4):
        self.outFile.write('401' + str(self.meshQuadrilateralNumber) + ',' + str(nodeTag1)+ ',' + str(nodeTag2) + ',' + str(nodeTag3) + ',' + str(nodeTag4) + '\n')
        self.meshQuadrilateralNumber = self.meshQuadrilateralNumber + 1; 

    def build(self):
        self.outFile.close
