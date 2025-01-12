#!/bin/bash

procs=16

#cp system/coarse/* system/
#cp system/fine/* system/


decomposePar
mpirun -np $procs --allow-run-as-root snappyHexMesh -overwrite -parallel -dict system/snappyHexMeshDict_layering |& grep -v "Read -1" | tee log.shm_layering
mpirun -np $procs --allow-run-as-root checkMesh -overwrite -parallel |& grep -v "Read -1" | tee log.checkmesh
reconstructParMesh -constant

# rm -rf 0 > /dev/null 2>&1
# cp -r 0_org 0 > /dev/null 2>&1

# rm -r processor*


# checkMesh | tee log.checkmesh
