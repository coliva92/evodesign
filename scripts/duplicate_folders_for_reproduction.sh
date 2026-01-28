#!/bin/bash

# Check if correct number of arguments provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <source_parent_folder> <destination_parent_folder>"
    exit 1
fi

SOURCE_PARENT=$1
DEST_PARENT=$2

# Check if source directory exists
if [ ! -d "$SOURCE_PARENT" ]; then
    echo "Error: Source directory '$SOURCE_PARENT' does not exist."
    exit 1
fi

# Loop through each item in the source directory
# The */ pattern restricts the loop to directories only
for dir_path in "$SOURCE_PARENT"/*/; do
    # Check if the glob expansion actually found a directory
    # (Handles the case where the folder is empty)
    if [ -d "$dir_path" ]; then
        
        # 1. Extract the folder name (e.g., /path/to/data/run1 -> run1)
        # We strip the trailing slash first, then get basename
        folder_name=$(basename "${dir_path%/}")

        # 2. Define the new destination path for this specific subfolder
        target_dest="$DEST_PARENT/$folder_name"
        
        # 3. Call the Python script
        # Assuming duplicate.py is in the current directory
        python3 duplicate_for_reproduction.py "$dir_path" "$target_dest"
    fi
done
