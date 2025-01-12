#!/bin/bash 

source ../personal_variables.sh


REF_GEOMETRY_DIR="$GEN_ROOT_DIR/geometries" 
REF_CONFIGS="$GEN_ROOT_DIR/scripts/reference_configs" 
SCRIPTS="$GEN_ROOT_DIR/scripts"
OPENFOAM_DIR="simpleFoam_MRF_propeller_komegasst_combined"

c=$1

workdir=${PWD}
echo $workdir 

dir_prefix=`echo $c | sed "s/\.stl//g"`
echo $dir_prefix 

rm -r $dir_prefix
mkdir $dir_prefix
cp -rfp  $REF_CONFIGS/$OPENFOAM_DIR $dir_prefix
mv $dir_prefix/$OPENFOAM_DIR $dir_prefix/simpleFoam_steady

cd $SCRIPTS

cp -rfp  $REF_GEOMETRY_DIR/$c/*.stl $dir_prefix/simpleFoam_steady/constant/triSurface/
cp -rfp  $REF_GEOMETRY_DIR/$c/case_info.txt $dir_prefix

cp create_mesh_inputs_case.py $c
cp mesh_input.sh $c
cd $c
bash mesh_input.sh
rm mesh_input.sh
rm create_mesh_inputs_case.py

