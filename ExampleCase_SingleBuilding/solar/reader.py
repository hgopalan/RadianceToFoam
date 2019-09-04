import os 
import numpy as np
from subprocess import call
print("Creating boundary data for OpenFOAM from Radiance heat-flux")
print("If OpenFOAM crashes, change the sampling variable from 1 to 10..20.. and rerun")


sampling=1
call(["rm","-rf","boundaryData/"])
source=open("boundary.info",'r')
patchnum=int(source.readline())
print("Number of patches:",patchnum)
for i in range (0,patchnum):
    patchname="boundaryData/"+source.readline()
    if "\n" in patchname:
#        print "Reading endline. So removing it."
        patchname=patchname[:-1]    
    print("Reading patch",patchname)
    pointsname=source.readline()
    if "\n" in pointsname:
#        print "Reading endline. So removing it."
        pointsname=pointsname[:-1]
#    print "Points File",pointsname
    data = np.loadtxt(pointsname,skiprows=8)
    x=data[:,0]
    y=data[:,1]
    z=data[:,2]
#    if(np.size(x)<=50000):
#        sampling=1
#    elif (np.size(x)>5e4 and np.size(x)<=1e5):
 #       sampling=int(1+9.0*(np.size(x)-5e4)/(1e5-5e4))
 #   else:
 #       print "OpenFOAM may crash. Sample size is huge"
 #       sampling=int(10+10.0*(np.size(x)-1e5)/(5e6-1e5))
    call(["mkdir","-p",patchname]);
    target = open(patchname+"/points", 'w')    
    openfoam_header_string1="/*--------------------------------*- C++ -*----------------------------------*\n"
    target.write(openfoam_header_string1)
    openfoam_header_string1="| =========                 |                                                 |\n"
    target.write(openfoam_header_string1)
    openfoam_header_string1="| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox\n"           
    target.write(openfoam_header_string1)
    openfoam_header_string1="|  \\    /   O peration     | Version:  2.3.0                                 |\n"
    target.write(openfoam_header_string1)
    openfoam_header_string1="|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |\n"
    target.write(openfoam_header_string1)
    openfoam_header_string1="|    \\/     M anipulation  |                                                 |\n"
    target.write(openfoam_header_string1)
    openfoam_header_string1="\*---------------------------------------------------------------------------*/\n"
    target.write(openfoam_header_string1)
    openfoam_header_string1="FoamFile\n"
    target.write(openfoam_header_string1)
    openfoam_header_string1="{\n"
    target.write(openfoam_header_string1)
    openfoam_header_string1="    version     2.0;\n"
    target.write(openfoam_header_string1)
    openfoam_header_string1="    format      ascii;\n"
    target.write(openfoam_header_string1)
    openfoam_header_string1="    class       vectorField;\n"
    target.write(openfoam_header_string1)
    openfoam_header_string1="    object      points;\n"
    target.write(openfoam_header_string1)
    openfoam_header_string1="}\n"
    target.write(openfoam_header_string1)
    openfoam_header_string1="// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n"
    target.write(openfoam_header_string1)
    target.write("\n\n\n")
    target.write("(\n")
    counter=0
    for j in range(0,np.size(x),sampling):
        counter=counter+1
        target.write('%s %f %f %f %s\n' %("( ",x[j],y[j],z[j],")"));
        #        target.write('%s %f %f %f %s\n' %("( ",x[i],y[i],0.0,")"));
    target.write(")\n")
    target.close()
    # Adjust for Albedo values .....
    if (patchname=="boundaryData/lower" or patchname=="boundaryData/terrain"):
        u1=(1-0.213)*data[:,3]
    elif (patchname=="boundaryData/road"):
        u1=(1-0.1587)*data[:,3]
        for i in range(0,len(u1)):
            if(u1[i]<1):
                u1[i]=np.average(u1)
    elif (patchname=="boundaryData/waterbody"):
        u1=1.0/25000*data[:,3]
        for i in range(0,len(u1)):
            if(u1[i]<1):
                u1[i]=np.average(u1)
    elif (patchname=="building" or patchname=="wall"):
        u1=(1-0.64)*data[:,3]        
    elif(patchname=="greenery"):
        u1=(1-0.24)*data[:,3]
        if(u1[i]<1):
            u1[i]=np.average(u1)
    else:
        u1=data[:,3]
    call(["mkdir",patchname+"/0"]);
    target = open(patchname+"/0/T", 'w')
    openfoam_header_string1="/*--------------------------------*- C++ -*----------------------------------*\n"
    target.write(openfoam_header_string1)
    openfoam_header_string1="| =========                 |                                                 |\n"
    target.write(openfoam_header_string1)
    openfoam_header_string1="| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox\n"           
    target.write(openfoam_header_string1)
    openfoam_header_string1="|  \\    /   O peration     | Version:  2.3.0                                 |\n"
    target.write(openfoam_header_string1)
    openfoam_header_string1="|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |\n"
    target.write(openfoam_header_string1)
    openfoam_header_string1="|    \\/     M anipulation  |                                                 |\n"
    target.write(openfoam_header_string1)
    openfoam_header_string1="\*---------------------------------------------------------------------------*/\n"
    target.write(openfoam_header_string1)
    openfoam_header_string1="FoamFile\n"
    target.write(openfoam_header_string1)
    openfoam_header_string1="{\n"
    target.write(openfoam_header_string1)
    openfoam_header_string1="    version     2.0;\n"
    target.write(openfoam_header_string1)
    openfoam_header_string1="    format      ascii;\n"
    target.write(openfoam_header_string1)
    openfoam_header_string1="    class       scalarAverageField;\n"
    target.write(openfoam_header_string1)
    openfoam_header_string1="    object      T;\n"
    target.write(openfoam_header_string1)
    openfoam_header_string1="}\n"
    target.write(openfoam_header_string1)
    openfoam_header_string1="// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n"
    target.write(openfoam_header_string1)
    target.write("\n\n\n")
    target.write('%i\n' % (counter))
    target.write("(\n\n\n")
    for j in range(0,np.size(x),sampling):
        target.write('%f\n' % (u1[j]))
    target.write(")\n")    
    target.close()
        
os.system("rm -rf ../constant/boundaryData")
os.system("mv boundaryData/ ../constant/")
source.close()

