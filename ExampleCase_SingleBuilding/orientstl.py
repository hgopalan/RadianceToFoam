import stl 
import numpy as np 
import sys
from stl.mesh import Mesh
print(sys.argv)
filename=str(sys.argv[1])
stl_mesh = Mesh.from_file(filename)
data = np.zeros(stl_mesh.data.size, dtype=Mesh.dtype)
new_mesh = Mesh(data)
new_mesh.normals[:] = stl_mesh.normals
new_mesh.vectors[:] = stl_mesh.vectors
for i in range(0,new_mesh.normals.shape[0]):
	nsum=0
	x1=new_mesh.vectors[i,0,0]
	x2=new_mesh.vectors[i,1,0]
	x3=new_mesh.vectors[i,2,0]
	y1=new_mesh.vectors[i,0,1]
	y2=new_mesh.vectors[i,1,1]
	y3=new_mesh.vectors[i,2,1]
	z1=new_mesh.vectors[i,0,2]
	z2=new_mesh.vectors[i,1,2]
	z3=new_mesh.vectors[i,2,2]
	# Shoe-Lace
	if(abs(z1-z2)<0.5 and abs(z2-z3)<0.5 and abs(z3-z1)<0.5):
		nsum=nsum+(x2-x1)*(y1+y2)
		nsum=nsum+(x3-x2)*(y3+y2)
		nsum=nsum+(x1-x3)*(y1+y3)
		if(nsum<=0):
			continue
		else:
			new_mesh.vectors[i,1,:]=stl_mesh.vectors[i,2,:]
			new_mesh.vectors[i,2,:]=stl_mesh.vectors[i,1,:]
new_mesh.save(sys.argv[1],mode=stl.Mode.ASCII)
