import stl 
import numpy as np 
import sys
from stl.mesh import Mesh
print((sys.argv))
stl_mesh = Mesh.from_file(sys.argv[1])
data = np.zeros(stl_mesh.data.size, dtype=Mesh.dtype)
new_mesh = Mesh(data)
new_mesh.normals[:] = stl_mesh.normals
new_mesh.vectors[:] = stl_mesh.vectors
target=open(sys.argv[2],"w")
for i in range(0,new_mesh.normals.shape[0]):
	x1=new_mesh.vectors[i,0,0]
	x2=new_mesh.vectors[i,1,0]
	x3=new_mesh.vectors[i,2,0]
	y1=new_mesh.vectors[i,0,1]
	y2=new_mesh.vectors[i,1,1]
	y3=new_mesh.vectors[i,2,1]
	z1=new_mesh.vectors[i,0,2]
	z2=new_mesh.vectors[i,1,2]
	z3=new_mesh.vectors[i,2,2]
	xc=1.0/3.0*(x1+x2+x3)
	yc=1.0/3.0*(y1+y2+y3)
	zc=1.0/3.0*(z1+z2+z3)
	nx=new_mesh.normals[i,0]
	ny=new_mesh.normals[i,1]
	nz=new_mesh.normals[i,2]
	normalmag=np.sqrt(pow(nx,2)+pow(ny,2)+pow(nz,2))+1e-15
	target.write("%f %f %f %f %f %f\n"%(xc,yc,zc,nx/normalmag,ny/normalmag,nz/normalmag))
target.close()
