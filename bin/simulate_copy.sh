#! /bin/bash

source_dir="/root/input"
destination_dir="/data/bff-test/data/input"
minutes_between_copies=20
seconds_between_copies=$((minutes_between_copies * 60))

# Find all the files in the source directory that end in .pod5
# Sort files numerically by the number between last _ and .pod5
files=($(find $source_dir -type f -name "*.pod5" | sort -V))

# Alternative with more explicit numeric extraction:
# files=($(find $source_dir -type f -name "*.pod5" | awk -F'[_.]' '{print $(NF-1) " " $0}' | sort -n | cut -d' ' -f2-))

echo "files: ${files[*]}"

# Every $minutes_between_copies minutes, copy the files to the destination directory
for file in "${files[@]}"; do
    echo "Copying $file to $destination_dir"
    cp "$file" "$destination_dir"
    
    # Wait for specified time or until Enter is pressed
    echo "Waiting $minutes_between_copies minutes before next copy... (press Enter to skip)"
    
    elapsed=0
    while [ $elapsed -lt $seconds_between_copies ]; do
        read -t 1 -n 1 input
        if [ $? -eq 0 ]; then
            echo "Skipping wait time..."
            break
        fi
        ((elapsed++))
        echo -ne "\rTime remaining: $((seconds_between_copies - elapsed)) seconds     "
    done
    echo # Add a new line after the countdown
done
