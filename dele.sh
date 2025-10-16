#!/bin/bash
# Usage: ./delete_sz200.sh /path/to/folder



target_dir="data/evolve_beta"

# Check if the directory exists
if [ ! -d "$target_dir" ]; then
    echo "Error: '$target_dir' is not a valid directory."
    exit 1
fi

# Find and delete files containing "sz_200" in the name
echo "Searching for files with 'sz_200' in their name under: $target_dir"
find "$target_dir" -type f -name '*sz_200*' -print -delete

echo "Deletion complete."
