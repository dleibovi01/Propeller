export case1=$1
export case2=$2

for case in $( eval echo {$case1..$case2} )
do
        cp cases/$case/case_info.csv .
        python propeller.py --case $case
        mkdir -p propeller_geometries/${case}
        cp propeller.stl case_info.csv case_info.txt propeller_geometries/${case}
        rm *stl
        rm case_info.csv
done

