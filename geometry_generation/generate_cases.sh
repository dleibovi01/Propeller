export case1=$1
export case2=$2

for case in $( eval echo {$case1..$case2} )
do
        mkdir -p cases/$case
        python csv_generation.py
        python write_speeds.py
        cp case_info.csv cases/$case/
        rm case_info.csv
        cp case_info.txt cases/$case/
        rm case_info.txt        
done