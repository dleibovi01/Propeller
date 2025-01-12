
export case1=$1
export case2=$2

for case in $( eval echo {$case1..$case2} )
do
        bash launch_setup_and_job ${case}
done

