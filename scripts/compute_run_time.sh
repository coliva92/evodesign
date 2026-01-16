#!/bin/bash

# --- Input Handling ---
# Check if an argument was provided
if [ -z "$1" ]; then
    echo "Usage: $0 <parent_directory>"
    exit 1
fi

PARENT_DIR="$1"
FILE_CREATION="initial_rng_state.txt"   # File to check creation date
FILE_MOD="pymoo_algorithm.bin"        # File to check modification date

# Arrays and variables to store data
declare -a differences
count=0
sum=0

# --- 1. Loop through subdirectories ---
# We verify that PARENT_DIR exists before looping
if [ ! -d "$PARENT_DIR" ]; then
    echo "Error: Directory '$PARENT_DIR' does not exist."
    exit 1
fi

echo "Processing subfolders in: $PARENT_DIR"

for dir in "$PARENT_DIR"/*/; do
    
    # Check if the specific files exist in this subdirectory
    if [[ -f "$dir$FILE_CREATION" && -f "$dir$FILE_MOD" ]]; then
        
        # Get timestamps (Linux syntax)
        # %W = Creation time (0 if unknown), %Y = Modification time
        t_create=$(stat -c %W "$dir$FILE_CREATION")
        t_mod=$(stat -c %Y "$dir$FILE_MOD")

        # Validation: Ensure creation time is supported/found (not 0)
        if [ "$t_create" -ne 0 ]; then
            # Calculate difference
            diff=$((t_mod - t_create))
            
            # Store in array and add to sum
            differences+=($diff)
            sum=$((sum + diff))
            ((count++))
        fi
    fi
done

# --- 2. Statistical Calculation ---

if [ "$count" -eq 0 ]; then
    echo "No valid subfolders found containing both files with readable creation dates."
    exit 0
fi

# Calculate Average (Mean) using bc for floating point math
mean=$(echo "scale=4; $sum / $count" | bc)

# Calculate Sum of Squared Differences for Standard Deviation
sum_sq_diff=0
for val in "${differences[@]}"; do
    # (val - mean)^2
    sq_diff=$(echo "scale=4; ($val - $mean)^2" | bc)
    sum_sq_diff=$(echo "scale=4; $sum_sq_diff + $sq_diff" | bc)
done

# Calculate Standard Deviation (Sample SD)
if [ "$count" -gt 1 ]; then
    std_dev=$(echo "scale=4; sqrt($sum_sq_diff / ($count - 1))" | bc)
else
    std_dev=0
fi

# --- 3. Output Results ---
echo "-----------------------------------"
echo "Subfolders processed: $count"
echo "Average time difference: $mean seconds"
echo "Standard Deviation:      $std_dev seconds"
echo "-----------------------------------"
