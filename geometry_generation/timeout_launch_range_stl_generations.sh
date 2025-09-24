export case1=$1
export case2=$2

for case in $( eval echo {$case1..$case2} )
do
	# Define the timeout in seconds
	TIMEOUT=120
	
	TASK="
	    cp cases/$case/case_info.csv .
	    cp cases/$case/case_info.txt .;
	    python propeller.py --case $case;
	    mkdir -p propeller_geometries/${case};
	    cp propeller.stl case_info.csv case_info.txt propeller_geometries/${case};
	    rm *stl;
	    rm case_info.csv;
	    rm case_info.txt;
	"


	# Start the task in the background
	bash -c "$TASK" &
	TASK_PID=$!

	# Function to kill all child processes of the TASK_PID
	kill_children() {
	    # Find all child processes of the TASK_PID
	    CHILD_PIDS=$(ps --ppid $TASK_PID -o pid=)
	    
	    # Kill all child processes
	    for pid in $CHILD_PIDS; do
		kill -9 $pid
	    done
	}

	# Wait for the task to complete or timeout
	(
	    sleep $TIMEOUT
	    if ps -p $TASK_PID > /dev/null; then
		echo "Task exceeded timeout. Killing the task and all its children..."
		kill -9 $TASK_PID
		kill $(pgrep -f 'python propeller.py')
		# kill_children  # Kill all child processes of the main task
	    fi
	) 


	# Wait for the task to finish normally
	wait $TASK_PID

done

