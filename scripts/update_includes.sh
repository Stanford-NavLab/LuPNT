#!/bin/bash

# Navigate to the lupnt directory
cd "$(dirname "$0")"/../lupnt

# Initialize or clear the lupnt.h file content before appending
: > lupnt.h

# Find all subdirectories and process them
find . -type d ! -path . | sort | while read -r subdir; do
    # Format the directory name for a comment
    subdir_name=$(echo "$subdir" | sed 's|./||')
    echo "// $subdir_name" >> lupnt.h

    # Find all .h files in the current subdirectory
    find "$subdir" -maxdepth 1 -type f -name "*.h" | sort | while read -r file; do
        # Create relative include path
        include_path=$(echo "$file" | sed 's|./||')
        # Append include directive to the lupnt.h file
        echo "#include \"$include_path\"" >> lupnt.h
    done

    # Add an empty line for readability between subdirectory sections
    echo "" >> lupnt.h
done

echo "A single header file lupnt.h has been created with all includes."
