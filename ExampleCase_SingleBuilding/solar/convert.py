import numpy as np
import sys
import os

cwd = os.getcwd()

if "/" not in sys.argv[1]:
	filename = cwd + "/" + sys.argv[1] 
else:
	filename = sys.argv[1]
	
if "/" not in sys.argv[2]:	
	outputfile = cwd + "/" + sys.argv[2]
else:
	outputfile = sys.argv[2]

data = np.genfromtxt(filename,skip_header=8)
npts = len(data)

oF = outputfile
oFW = open(oF,"w")

oFW.write("# vtk DataFile Version 2.0\n")
oFW.write("Radiation Building\n")
oFW.write("ASCII\n")
oFW.write("DATASET POLYDATA\n")

oFW.write("POINTS " + str(npts) + " float\n")

for i in range(0,npts):
	x = data[i][0]
	y = data[i][1]
	z = data[i][2]
	oFW.write(str(x) + " " + str(y) + " " + str(z) + "\n")

oFW.write("VERTICES " + str(npts) + " " + str(2*npts) + "\n")
for i in range(0,npts):
	oFW.write("1 " + str(i) + "\n")

oFW.write("POINT_DATA " + str(npts) + "\n")
oFW.write("SCALARS radiation float 1\n")
oFW.write("LOOKUP_TABLE default\n")
for i in range(0,npts):
	rad = data[i][3]
	oFW.write(str(rad) + "\n")
oFW.close()
