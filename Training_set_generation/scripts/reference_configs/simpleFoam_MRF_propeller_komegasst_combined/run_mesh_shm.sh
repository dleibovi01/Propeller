#!/bin/bash

procs=32



# surfaceTransformPoints -rotate-y 180 constant/triSurface/aero_suv_1.stl constant/triSurface/aero_suv_2.stl

foamCleanTutorials
rm -rf 0 > /dev/null 2>&1
cp -r 0_org 0 > /dev/null 2>&1

surfaceFeatureExtract | tee log.surfaceFeatures
blockMesh | tee log.blockMesh
topoSet | tee



decomposePar
mpirun -np $procs --allow-run-as-root snappyHexMesh -overwrite -parallel |& grep -v "Read -1" | tee log.shm
mpirun -np $procs --allow-run-as-root snappyHexMesh -overwrite -parallel -dict system/snappyHexMeshDict_layering |& grep -v "Read -1" | tee log.shm_layering
mpirun -np $procs --allow-run-as-root checkMesh -parallel |& grep -v "Read -1" | tee log.checkmesh



reconstructParMesh -constant
rm -rf 0 > /dev/null 2>&1
cp -r 0_org 0 > /dev/null 2>&1

rm -r processor*


checkMesh | tee log.checkmesh


