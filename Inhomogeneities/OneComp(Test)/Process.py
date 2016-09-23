from math   import *
import glob
import os

shear     = 2.5
Theta     = 0.0
tumble    = 0.5

FileFinalName2 = "Alignment(" + str(shear) + "," + str(Theta) + "," + str(tumble) + ").txt" 

read_files = glob.glob("*.dat")

# Calculate and print q = sqrt(q0**2 + q1**2 + q2**2 + q3**2 + q4**2 )

with open (FileFinalName2,"wb") as outfile:
    for f in read_files:
        file = open(f,"r")
        lines = file.readlines()
        for i in range(1,len(lines)):
            Split = lines[i].split()
            float_list = [float(i) for i in Split]
            
            y  = float_list[0]
            q0 = float_list[1]
            q1 = float_list[2]
            q2 = float_list[3]
            q3 = float_list[4]
            q4 = float_list[5]
            
            align = sqrt(q0**2 + q1**2 + q2**2 + q3**2 + q4**2)
            
            outfile.write(f[4:8] + "\t" + str(y) + "\t" + str(align) + "\n")
        outfile.write("\n")

os.system("rm *.dat")
