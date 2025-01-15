#!/bin/bash

# This script generates a toy example pairs file from 4DN data
# Source: https://data.4dnucleome.org/files-processed/4DNFI9DUBGCW/

# Check if input file exists
if [ ! -f "4DNFI9DUBGCW.pairs.gz" ]; then
    echo "Error: Input file 4DNFI9DUBGCW.pairs.gz not found"
    exit 1
fi

# Extract first 100k lines to create toy example
echo "Generating toy example pairs file..."
zcat 4DNFI9DUBGCW.pairs.gz | head -n 100000 > toyExample.pair

echo "Done! Toy example saved as toyExample.pair"