# Shear flow inhomogeneities (Two Components: R. Lugo-Frias, H. Reinken)
# This file generates raw data

from fipy   import *
from random import *
from math   import *
import os

L        = 1.
nx       = 200
dx       = L/nx
mesh     = Grid1D(nx=nx, dx=dx,overlap=2)

# Important parameters
D1         = 0.001  # 1/Er number
D2         = 0.002  # 1/Er number
shear     = 7.5
sigma     = 0.0
ThetaA    = -0.25
ThetaB    = -0.25
ThetaAB   = 1.5
tumbleA   = 1.2
tumbleB   = 3.0

dt        = 0.9*dx**2/(2*0.5*10*D1)
#dt        = 0.01
steps     = 2000


# Random initial conditions for qi0--qi4
x         = mesh.cellCenters[0]

##################
# Define (A component)

qa0        = CellVariable(name=r"$qa0$", mesh=mesh)
qa1        = CellVariable(name=r"$qa1$", mesh=mesh)
qa2        = CellVariable(name=r"$qa2$", mesh=mesh)
qa3        = CellVariable(name=r"$qa3$", mesh=mesh)
qa4        = CellVariable(name=r"$qa4$", mesh=mesh)

###################
## Define (B component)

qb0        = CellVariable(name=r"$qb0$", mesh=mesh)
qb1        = CellVariable(name=r"$qb1$", mesh=mesh)
qb2        = CellVariable(name=r"$qb2$", mesh=mesh)
qb3        = CellVariable(name=r"$qb3$", mesh=mesh)
qb4        = CellVariable(name=r"$qb4$", mesh=mesh)

##################
# Set initial conditions 

### Zero
# A component
#qa0.value  = 0.0
#qa1.value  = 0.0
#qa2.value  = 0.0
#qa3.value  = 0.0
#qa4.value  = 0.0
# B component
#qb0.value  = 0.0
#qb1.value  = 0.0
#qb2.value  = 0.0
#qb3.value  = 0.0
#qb4.value  = 0.0

### Random
## A component
#qa0.setValue(random())
#qa1.setValue(random())
#qa2.setValue(random())
#qa3.setValue(random())
#qa4.setValue(random())
## B Component
#qb0.setValue(random())
#qb1.setValue(random())
#qb2.setValue(random())
#qb3.setValue(random())
#qb4.setValue(random())

# By parts in the box
# A component
qa0.setValue(random(), where=(x > L/2.))
qa1.setValue(random(), where=(x > L/2.))
qa2.setValue(random(), where=(x > L/2.))
qa3.setValue(random(), where=(x > L/2.))
qa4.setValue(random(), where=(x > L/2.))
# B component
qb0.setValue(random(), where=(x < L/2.))
qb1.setValue(random(), where=(x < L/2.))
qb2.setValue(random(), where=(x < L/2.))
qb3.setValue(random(), where=(x < L/2.))
qb4.setValue(random(), where=(x < L/2.))

# Boundary conditions, not given (no--flux)
valueTop    = 1.5*sqrt(3./2.)
valueBottom = 1.5*sqrt(3./2.)

qa0.constrain(valueTop, mesh.facesRight)
qa0.constrain(valueBottom, mesh.facesLeft)
qa1.constrain(valueTop, mesh.facesRight)
qa1.constrain(valueBottom, mesh.facesLeft)
qa2.constrain(valueTop, mesh.facesRight)
qa2.constrain(valueBottom, mesh.facesLeft)
qa3.constrain(valueTop, mesh.facesRight)
qa3.constrain(valueBottom, mesh.facesLeft)
qa4.constrain(valueTop, mesh.facesRight)
qa4.constrain(valueBottom, mesh.facesLeft)

qb0.constrain(valueTop, mesh.facesRight)
qb0.constrain(valueBottom, mesh.facesLeft)
qb1.constrain(valueTop, mesh.facesRight)
qb1.constrain(valueBottom, mesh.facesLeft)
qb2.constrain(valueTop, mesh.facesRight)
qb2.constrain(valueBottom, mesh.facesLeft)
qb3.constrain(valueTop, mesh.facesRight)
qb3.constrain(valueBottom, mesh.facesLeft)
qb4.constrain(valueTop, mesh.facesRight)
qb4.constrain(valueBottom, mesh.facesLeft)

# View the initial conditions
#if __name__ == '__main__':
 #viewer = Viewer(vars=(qa0,qa1,qa2,qa3,qa4,qb0,qb1,qb2,qb3,qb4), datamin = -0.1, datamax = 1.5)

####################################################################
# Free energy density terms
qa_sq     = qa0**2 + qa1**2 + qa2**2 + qa3**2 + qa4**2
qb_sq     = qb0**2 + qb1**2 + qb2**2 + qb3**2 + qb4**2
qmix      = qa0*qb0 + qa1*qb1 + qa2*qb2 + qa3*qb3 + qa4*qb4

Phi0A     = (ThetaA - 3.*qa0 + 2.*qa_sq)*qa0 + 3.*(qa1**2 + qa2**2) - 1.5*(qa3**2 + qa4**2)
Phi1A     = (ThetaA + 6.*qa0 + 2.*qa_sq)*qa1 - 1.5*sqrt(3.)*(qa3**2 - qa4**2)
Phi2A     = (ThetaA + 6.*qa0 + 2.*qa_sq)*qa2 - 3.*sqrt(3.)*qa3*qa4
Phi3A     = (ThetaA - 3.*qa0 + 2.*qa_sq)*qa3 - 3.*sqrt(3.)*(qa1*qa3 + qa2*qa4)
Phi4A     = (ThetaA - 3.*qa0 + 2.*qa_sq)*qa4 + 3.*sqrt(3.)*(qa1*qa4 - qa2*qa3)

Phi0B     = (ThetaB - 3.*qb0 + 2.*qb_sq)*qb0 + 3.*(qb1**2 + qb2**2) - 1.5*(qb3**2 + qb4**2)
Phi1B     = (ThetaB + 6.*qb0 + 2.*qb_sq)*qb1 - 1.5*sqrt(3.)*(qb3**2 - qb4**2)
Phi2B     = (ThetaB + 6.*qb0 + 2.*qb_sq)*qb2 - 3.*sqrt(3.)*qb3*qb4
Phi3B     = (ThetaB - 3.*qb0 + 2.*qb_sq)*qb3 - 3.*sqrt(3.)*(qb1*qb3 + qb2*qb4)
Phi4B     = (ThetaB - 3.*qb0 + 2.*qb_sq)*qb4 + 3.*sqrt(3.)*(qb1*qb4 - qb2*qb3)

PhiMixA   = -qb0*(ThetaAB + qmix)
PhiMixA   = -qb1*(ThetaAB + qmix)
PhiMixA   = -qb2*(ThetaAB + qmix)
PhiMixA   = -qb3*(ThetaAB + qmix)
PhiMixA   = -qb4*(ThetaAB + qmix)

PhiMixB   = -qa0*(ThetaAB + qmix)
PhiMixB   = -qa1*(ThetaAB + qmix)
PhiMixB   = -qa2*(ThetaAB + qmix)
PhiMixB   = -qa3*(ThetaAB + qmix)
PhiMixB   = -qa4*(ThetaAB + qmix)

Sqa0      = -(Phi0A + PhiMixA) - sqrt(1./3.)*shear*sigma*qa2
Sqa1      = -(Phi1A + PhiMixA) + shear*qa2
Sqa2      = -(Phi2A + PhiMixA) - sqrt(1./3.)*shear*sigma*qa0 - shear*qa1 + tumbleA*shear
Sqa3      = -(Phi3A + PhiMixA) + (1./2.)*shear*(sigma + 1.)*qa4
Sqa4      = -(Phi4A + PhiMixA) + (1./2.)*shear*(sigma - 1.)*qa3

Sqb0      = -(Phi0B + PhiMixB) - sqrt(1./3.)*shear*sigma*qb2
Sqb1      = -(Phi1B + PhiMixB) + shear*qb2
Sqb2      = -(Phi2B + PhiMixB) - sqrt(1./3.)*shear*sigma*qb0 - shear*qb1 + tumbleB*shear
Sqb3      = -(Phi3B + PhiMixB) + (1./2.)*shear*(sigma + 1.)*qb4
Sqb4      = -(Phi4B + PhiMixB) + (1./2.)*shear*(sigma - 1.)*qb3

####################################################################
# Inhomogeneous equations
eqqa0      = TransientTerm() == DiffusionTerm(coeff=D1) + Sqa0
eqqa1      = TransientTerm() == DiffusionTerm(coeff=D1) + Sqa1
eqqa2      = TransientTerm() == DiffusionTerm(coeff=D1) + Sqa2
eqqa3      = TransientTerm() == DiffusionTerm(coeff=D1) + Sqa3
eqqa4      = TransientTerm() == DiffusionTerm(coeff=D1) + Sqa4

eqqb0      = TransientTerm() == DiffusionTerm(coeff=D2) + Sqb0
eqqb1      = TransientTerm() == DiffusionTerm(coeff=D2) + Sqb1
eqqb2      = TransientTerm() == DiffusionTerm(coeff=D2) + Sqb2
eqqb3      = TransientTerm() == DiffusionTerm(coeff=D2) + Sqb3
eqqb4      = TransientTerm() == DiffusionTerm(coeff=D2) + Sqb4

#The actual integration of the equation is done here
read_filesA = []
read_filesB = []

for step in range(steps):
  
    eqqa0.solve(var=qa0,dt=dt, solver=DummySolver())
    eqqa1.solve(var=qa1,dt=dt, solver=DummySolver())
    eqqa2.solve(var=qa2,dt=dt, solver=DummySolver())
    eqqa3.solve(var=qa3,dt=dt, solver=DummySolver())
    eqqa4.solve(var=qa4,dt=dt, solver=DummySolver())
    
    eqqb0.solve(var=qb0,dt=dt, solver=DummySolver())
    eqqb1.solve(var=qb1,dt=dt, solver=DummySolver())
    eqqb2.solve(var=qb2,dt=dt, solver=DummySolver())
    eqqb3.solve(var=qb3,dt=dt, solver=DummySolver())
    eqqb4.solve(var=qb4,dt=dt, solver=DummySolver())
        
    nameA = 'dataA' + str('%04d' % step) + '.dat'
    nameB = 'dataB' + str('%04d' % step) + '.dat'
        
    TSVViewer(vars = (qa0, qa1, qa2, qa3, qa4)).plot(filename=nameA)
    TSVViewer(vars = (qb0, qb1, qb2, qb3, qb4)).plot(filename=nameB)
    
    read_filesA.append(nameA)
    read_filesB.append(nameB)

# Two big files

FileFinalNameA = 'A_Shear_' + str(shear) + 'TumblingA_' + str(tumbleA) + 'TumblingB_' + str(tumbleB) + '_ThetaA_' + str(ThetaA) + '_ThetaB_' + str(ThetaB) + '_ThetaAB_' + str(ThetaAB) + '.txt'
FileFinalNameB = 'B_Shear_' + str(shear) + 'TumblingA_' + str(tumbleA) + 'TumblingB_' + str(tumbleB) + '_ThetaA_' + str(ThetaA) + '_ThetaB_' + str(ThetaB) + '_ThetaAB_' + str(ThetaAB) + '.txt'

with open (FileFinalNameA,"wb") as outfile:
    for f in read_filesA:
        file = open(f,"r")
        lines = file.readlines()
        for i in range(1,len(lines)):
            outfile.write(f[5:9] + "\t" + str(lines[i]) )
        outfile.write("\n")
        
with open (FileFinalNameB,"wb") as outfile:
    for f in read_filesB:
        file = open(f,"r")
        lines = file.readlines()
        for i in range(1,len(lines)):
            outfile.write(f[5:9] + "\t" + str(lines[i]) )
        outfile.write("\n")

FileFinalName2A = 'A_Alignment_Shear_' + str(shear) + 'TumblingA_' + str(tumbleA) + 'TumblingB_' + str(tumbleB) + '_ThetaA_' + str(ThetaA) + '_ThetaB_' + str(ThetaB) + '_ThetaAB_' + str(ThetaAB) + '.txt'
FileFinalName2B = 'B_Alignment_Shear_' + str(shear) + 'TumblingA_' + str(tumbleA) + 'TumblingB_' + str(tumbleB) + '_ThetaA_' + str(ThetaA) + '_ThetaB_' + str(ThetaB) + '_ThetaAB_' + str(ThetaAB) + '.txt'

# Calculate and print q = sqrt(q0**2 + q1**2 + q2**2 + q3**2 + q4**2 )

with open (FileFinalName2A,"wb") as outfile:
    for f in read_filesA:
        file = open(f,"r")
        lines = file.readlines()
        for i in range(1,len(lines)):
            Split = lines[i].split()
            float_list = [float(i) for i in Split]
            
            y   = float_list[0]
            qa0 = float_list[1]
            qa1 = float_list[2]
            qa2 = float_list[3]
            qa3 = float_list[4]
            qa4 = float_list[5]
            
            align = sqrt(qa0**2 + qa1**2 + qa2**2 + qa3**2 + qa4**2)
            
            outfile.write(f[5:9] + "\t" + str(y) + "\t" + str(align) + "\n")
        outfile.write("\n")
        

with open (FileFinalName2B,"wb") as outfile:
    for f in read_filesB:
        file = open(f,"r")
        lines = file.readlines()
        for i in range(1,len(lines)):
            Split = lines[i].split()
            float_list = [float(i) for i in Split]
            
            y  = float_list[0]
            qb0 = float_list[1]
            qb1 = float_list[2]
            qb2 = float_list[3]
            qb3 = float_list[4]
            qb4 = float_list[5]
            
            align = sqrt(qb0**2 + qb1**2 + qb2**2 + qb3**2 + qb4**2)
            
            outfile.write(f[5:9] + "\t" + str(y) + "\t" + str(align) + "\n")
        outfile.write("\n")

os.system("rm *.dat")
