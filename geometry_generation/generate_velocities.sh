export case1=$1
export case2=$2

for case in $( eval echo {$case1..$case2} )
do
        rm case_info.txt 
        python write_speeds.py
        cp case_info.txt propeller_geometries/$case/
        rm case_info.txt        
done
