#!/bin/bash

procs=32

# transformPoints -scale 0.001

topoSet
decomposePar -force

mpirun -np $procs --oversubscribe --allow-run-as-root simpleFoam -parallel |& grep -v "Read -1" | tee log.solver

reconstructPar

foamToVTK -latestTime
