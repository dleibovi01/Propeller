#!/bin/bash 


source $SCRATCH/miniconda3/bin/activate
pip install numpy

REF_GEOMETRY_DIR="$SCRATCH/datasets/driveSim/geometries" #"/scratch/10238/dleibovi/datasets/driveSim/geometries" 
REF_CONFIGS="$SCRATCH/automotive/Training_set_generation/dataset-generation/baseline_runs/scripts/reference_configs"  #"/work/10238/dleibovi/vista/automotive/Training_set_generation/dataset-generation/baseline_runs/scripts/reference_configs"
SCRIPTS="$SCRATCH/automotive/Training_set_generation/dataset-generation/baseline_runs/scripts" #"/work/10238/dleibovi/vista/automotive/Training_set_generation/dataset-generation/baseline_runs/scripts"
OPENFOAM_DIR="simpleFoam_steady_divided_floor_separate_wheels"

c=$1

workdir=${PWD}

echo $workdir 

cd $REF_GEOMETRY_DIR 


# echo $case_list

cd $workdir 

echo $PWD
# for c in $case_list; do


dir_prefix=`echo $c | sed "s/\.stl//g"`
echo $dir_prefix 

rm -r detailed_car_$dir_prefix
mkdir detailed_car_$dir_prefix
cp -rfp  $REF_CONFIGS/GPU_23M/$OPENFOAM_DIR detailed_car_$dir_prefix
mv detailed_car_$dir_prefix/$OPENFOAM_DIR detailed_car_$dir_prefix/simpleFoam_steady

# bash $REF_GEOMETRY_DIR/add_case_info.sh $c
rm $REF_GEOMETRY_DIR/$c/case_info.txt
rm $REF_GEOMETRY_DIR/$c/log.solve
cp $REF_GEOMETRY_DIR/case_info.txt $REF_GEOMETRY_DIR/$c
cp $REF_GEOMETRY_DIR/parea_script.sh $REF_GEOMETRY_DIR/$c
cd $REF_GEOMETRY_DIR/$c
surfaceMeshConvert body.obj body.stl
surfaceMeshConvert front_wheels.obj front_wheels.stl
surfaceMeshConvert back_wheels.obj back_wheels.stl
bash parea_script.sh
rm $REF_GEOMETRY_DIR/$c/log.solve
cd $SCRIPTS



cp -rfp  $REF_GEOMETRY_DIR/$c/*.stl detailed_car_$dir_prefix/simpleFoam_steady/constant/triSurface/
cp -rfp  $REF_GEOMETRY_DIR/$c/case_info.txt detailed_car_$dir_prefix

cp create_mesh_inputs_case.py detailed_car_$c
cp mesh_input.sh detailed_car_$c
cd detailed_car_$c
bash mesh_input.sh
rm mesh_input.sh
rm create_mesh_inputs_case.py

