#!/bin/bash

# Request the filepaths for all BED Files
read -a inputs -p "Enter BED File Paths (Space-Separated): "

# Check for Null Input:
if [[ ${#inputs[@]} -eq 0 ]]; then
    echo "Error: No BED Files provided."
    exit 1
fi

# Display all input BED file paths
echo "Provided BED File Paths:"
for file in "${inputs[@]}"; do
    echo "$file"
done

# Prompt for the scope
read -p "Enter the scope: " scope

# Display confirmation message
echo "Scope: $scope"
read -p "Is everything okay? Enter 'yes' to continue or 'no' to cancel: " response

if [[ $response == "yes" ]]; then
    # Run the Explore Function using the inputs and scope.
    python3 explore.py "$scope" "${inputs[@]}"
else
    echo "Operation canceled."
fi

# /Users/pranavmahableshwarkar/BrownUniversity/LarschanLab/Analytics/PM_NewScripts/dna_bed/CLAMP_KC_DNA.bed /Users/pranavmahableshwarkar/BrownUniversity/LarschanLab/Analytics/PM_NewScripts/dna_bed/CLAMP_S2_DNA.bed /Users/pranavmahableshwarkar/BrownUniversity/LarschanLab/Analytics/PM_NewScripts/dna_bed/gaf_chip.bed /Users/pranavmahableshwarkar/BrownUniversity/LarschanLab/Analytics/PM_NewScripts/dna_bed/MLE_DNA.bed
