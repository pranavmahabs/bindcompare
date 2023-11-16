#!/bin/bash

# Request the filepaths for all BED Files
read -a inputs -p "Enter BED File Paths (Space-Separated): "

# Check for Null Input:
if [[ ${#inputs[@]} -eq 0 ]]; then
    echo "Error: No BED Files provided."
    exit 1
fi

# Run the Explore Function using the inputs. 
python3 explore.py ${#inputs[@]} "${inputs[@]}"

