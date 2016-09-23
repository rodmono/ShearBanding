# Shear flow inhomogeneities (Two Components: R. Lugo-Frias, H. Reinken)
# This file generates raw data

from fipy   import *
from random import *
from math   import *
import os

L        = 1.
nx       = 200
dx       = L/nx
mesh     = Grid1D(nx=nx, dx=dx, overlap=5)

# Important parameters
D         = 1./1000.  # 1/Er number
shear     = 2.0
sigma     = 0.0
Theta     = -0.5
tumble    = 0.8

dt        = 0.9*dx**2/(2.*5.*D)
steps     = 1500

# Boundary conditions, not given (no--flux)
x         = mesh.cellCenters[0]

##################
# Define (A component)

q0        = CellVariable(name=r"$q0$", mesh=mesh)
q1        = CellVariable(name=r"$q1$", mesh=mesh)
q2        = CellVariable(name=r"$q2$", mesh=mesh)
q3        = CellVariable(name=r"$q3$", mesh=mesh)
q4        = CellVariable(name=r"$q4$", mesh=mesh)

##################
# Set initial conditions 

## Zero
q0.value  = 0.0
q1.value  = 0.0
q2.value  = 0.0
q3.value  = 0.0
q4.value  = 0.0

## Random
q0.setValue(random()-0.5)
q1.setValue(random()-0.5)
q2.setValue(random()-0.5)
q3.setValue(random()-0.5)
q4.setValue(random()-0.5)

## By parts in the box
## A component
#qa0.setValue(random(), where=(x > L/2.))
#qa1.setValue(random(), where=(x > L/2.))
#qa2.setValue(random(), where=(x > L/2.))
#qa3.setValue(random(), where=(x > L/2.))
#qa4.setValue(random(), where=(x > L/2.))

# Boundary conditions, not given (no--flux)
valueTop    = (3/2.)*sqrt(3/2.)
valueBottom = (3/2.)*sqrt(3/2.)
q0.constrain(valueTop, mesh.facesRight)
q0.constrain(valueBottom, mesh.facesLeft)
q1.constrain(valueTop, mesh.facesRight)
q1.constrain(valueBottom, mesh.facesLeft)
q2.constrain(valueTop, mesh.facesRight)
q2.constrain(valueBottom, mesh.facesLeft)
q3.constrain(valueTop, mesh.facesRight)
q3.constrain(valueBottom, mesh.facesLeft)
q4.constrain(valueTop, mesh.facesRight)
q4.constrain(valueBottom, mesh.facesLeft)

# View the initial conditions
#if __name__ == '__main__':
 #viewer = Viewer(vars=(q0,q1,q2,q3,q4), datamin = -0.1, datamax = 1.5)

#####################################################################
# Free energy density terms
q_sq     = q0**2 + q1**2 + q2**2 + q3**2 + q4**2

Phi0     = (Theta - 2.*q0 + 2.*q_sq)*q0 + 2.*(q1**2 + q2**2) - (q3**2 + q4**2)
Phi1     = (Theta + 4.*q0 + 2.*q_sq)*q1 - sqrt(3.)*(q3**2 - q4**2)
Phi2     = (Theta + 4.*q0 + 2.*q_sq)*q2 - (2.)*sqrt(3.)*q3*q4
Phi3     = (Theta - 2.*q0 + 2.*q_sq)*q3 - 2.*sqrt(3.)*(q1*q3 + q2*q4)
Phi4     = (Theta - 2.*q0 + 2.*q_sq)*q4 + 2.*sqrt(3.)*(q1*q4 - q2*q3)

Sq0      = -Phi0 - sqrt(1./3.)*shear*sigma*q2
Sq1      = -Phi1 + shear*q2
Sq2      = -Phi2 - sqrt(1./3.)*shear*sigma*q0 - shear*q1 + tumble*shear
Sq3      = -Phi3 + (1./2.)*shear*(sigma + 1.)*q4
Sq4      = -Phi4 + (1./2.)*shear*(sigma - 1.)*q3

#####################################################################
# Homogeneous equations
eqq0      = TransientTerm() == Sq0
eqq1      = TransientTerm() == Sq1
eqq2      = TransientTerm() == Sq2
eqq3      = TransientTerm() == Sq3
eqq4      = TransientTerm() == Sq4

### Inhomogeneous equations
#eqq0      = TransientTerm() == DiffusionTerm(coeff=D) + Sq0
#eqq1      = TransientTerm() == DiffusionTerm(coeff=D) + Sq1
#eqq2      = TransientTerm() == DiffusionTerm(coeff=D) + Sq2
#eqq3      = TransientTerm() == DiffusionTerm(coeff=D) + Sq3
#eqq4      = TransientTerm() == DiffusionTerm(coeff=D) + Sq4

#The actual integration of the equation is done here
read_files = []

for step in range(steps):
  
    eqq0.solve(var=q0,dt=dt, solver=DummySolver())
    eqq1.solve(var=q1,dt=dt, solver=DummySolver())
    eqq2.solve(var=q2,dt=dt, solver=DummySolver())
    eqq3.solve(var=q3,dt=dt, solver=DummySolver())
    eqq4.solve(var=q4,dt=dt, solver=DummySolver())

    name = 'data' + str('%05d' % step) + '.dat'

    TSVViewer(vars = (q0, q1, q2, q3, q4)).plot(filename=name)

    read_files.append(name)


# Two big files

FileFinalName = 'One_Tumbling_' + str(tumble) + '_Shear_' + str(shear) + '_Theta_' + str(Theta) + '.txt'

with open (FileFinalName,"wb") as outfile:
    for f in read_files:
        file = open(f,"r")
        lines = file.readlines()
        for i in range(1,len(lines)):
            outfile.write(f[4:9] + "\t" + str(lines[i]) )
        outfile.write("\n")
        

FileFinalName2 = "One_Alignment_Shear" + str(shear) + "_Theta_" + str(Theta) + "_Tumbling_" + str(tumble) + ".txt" 

# Calculate and print q = sqrt(q0**2 + q1**2 + q2**2 + q3**2 + q4**2 )

with open (FileFinalName2,"wb") as outfile:
    for f in read_files:
        file = open(f,"r")
        lines = file.readlines()
        for i in range(1,len(lines)):
            Split = lines[i].split()
            float_list = [float(i) for i in Split]
            
            y  = float_list[0]
            qa0 = float_list[1]
            qa1 = float_list[2]
            qa2 = float_list[3]
            qa3 = float_list[4]
            qa4 = float_list[5]
            
            align = sqrt(qa0**2 + qa1**2 + qa2**2 + qa3**2 + qa4**2)
            
            outfile.write(f[4:9] + "\t" + str(y) + "\t" + str(align) + "\n")
        outfile.write("\n")


os.system("rm *.dat")
