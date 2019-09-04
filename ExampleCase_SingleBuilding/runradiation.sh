 #!/bin/bash 
declare -a file=(  building.stl   lower.stl )
declare -a patch=(  building   lower )
arraylength=2
# OpenFOAM Conversion Part 
for (( i=1; i<${arraylength}+1; i++ ));
do
    echo ${file[$i-1]}
    PATCHNAME=${patch[$i-1]}
    FILENAME=${file[$i-1]}
    surfaceMeshTriangulate -patches "($PATCHNAME)" $FILENAME > log.tmp
    surfaceOrient $FILENAME "( -5000 -5000 10000)" $FILENAME > log.tmp
    python orientstl.py $FILENAME > log.log
    gmsh -2 $FILENAME -format stl > log.log
    grep -- "0 -1" log.log | wc -l
    python stlpts.py $FILENAME $PATCHNAME.pts 
    surfaceConvert $FILENAME $PATCHNAME.obj > log.tmp
done
# Delete existing data 
rm -rf PTS
mv *.stl solar/
mv *.obj solar/
mkdir PTS
mv *.pts PTS/
rm -rf solar/PTS
rm -rf solar/OUT
mv PTS/ solar/
# Run Radiance 
echo "Running Radiance"
cd solar/
# Radiance output file generated here ......
# Reader expects output file name patchname.out 
# Convert to OpenFOAM
echo "Converting Radiance to OpenFOAM"
python reader.py


