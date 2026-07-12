#!/bin/bash

# 1. Gather all physically existing .cpp files (except init.cpp) into a clean array
current_files=()
while IFS= read -r file; do
    # Remove leading ./
    clean_file="${file#./}"
    [ "$clean_file" == "init.cpp" ] && continue
    current_files+=("$clean_file")
done < <(find . -type f -name "*.cpp")

# 2. Parse the existing sources.list file to figure out what to keep and track deletions
retained_files=()
deleted_files=()

if [ -f sources.list ]; then
    while IFS= read -r line; do
        # Clean hidden carriage returns and whitespaces
        clean_line=$(echo "$line" | tr -d '\r' | xargs)
        [ -z "$clean_line" ] && continue
        
        # Check if this tracked file still exists on disk
        if [[ " ${current_files[*]} " =~ " ${clean_line} " ]]; then
            retained_files+=("$clean_line")
        else
            deleted_files+=("$clean_line")
        fi
    done < sources.list
fi

# 3. Find brand new files that aren't in the retained list yet
added_files=()
for file in "${current_files[@]}"; do
    if [[ ! " ${retained_files[*]} " =~ " ${file} " ]]; then
        added_files+=("$file")
    fi
done

# 4. Alphabetically sort ONLY the newly added files before appending them
if [ ${#added_files[@]} -gt 0 ]; then
    IFS=$'\n' sorted_added=($(sort <<<"${added_files[*]}"))
    unset IFS
else
    sorted_added=()
fi

# 5. Write everything back to sources.list
# First, write retained files in their EXACT original order
: > sources.list
for file in "${retained_files[@]}"; do
    echo "$file" >> sources.list
done

# Next, append only the new files to the bottom
for file in "${sorted_added[@]}"; do
    echo "$file" >> sources.list
done

# 6. Output Report
if [ ${#added_files[@]} -gt 0 ]; then
    echo "The following files were APPENDED to the end of sources.list:"
    for file in "${sorted_added[@]}"; do echo "  + $file"; done
fi

if [ ${#deleted_files[@]} -gt 0 ]; then
    echo "The following files were REMOVED from their positions in sources.list:"
    for file in "${deleted_files[@]}"; do echo "  - $file"; done
fi

if [ ${#added_files[@]} -eq 0 ] && [ ${#deleted_files[@]} -eq 0 ]; then
    echo "sources.list is already up to date. No changes made."
fi
