#!/bin/bash

export case1=$1
export case2=$2

case0=$case1


# Function to check if an .stl file exists in the specified directory
check_stl_file() {
    directory=$1  # Directory to check for .stl file
    case0=$2
    # Check if any .stl files are present in the directory
    if ls "$directory"/*.stl 1> /dev/null 2>&1; then
	echo "true"
    else
	echo "false"
    fi
}

for case in $( eval echo {$case1..$case2} )
do
	directory="propeller_geometries/$case"  # Replace with your directory path
	echo ${case}

	# result=$(check_stl_file "$directory")
	
	if ls "$directory"/*.stl 1> /dev/null 2>&1; then
		echo "true"
		mv propeller_geometries/$case propeller_geometries/$case0
		case0=$((case0+1))
	else
		echo "false"
		rm -r propeller_geometries/$case
	fi

done
