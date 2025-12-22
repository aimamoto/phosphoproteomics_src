#!/bin/bash

source ~/Documents/fastpace_alignment/.pyenv/bin/activate

# Set input and output directories
INPUT_DIR="/media/akira/argentee/proteome/251014_FragPipeAnalystR_1.1.0_MaxLFQ_250527/consensus/FASTA_DE"

OUTPUT_DIR="$INPUT_DIR/fastpace_out"  # You can change this if you want output elsewhere
mkdir $OUTPUT_DIR

# Loop through all files starting with the filename starting from "pY_" ending with '.fasta'
for input_file in "$INPUT_DIR"/ALLpY*.fasta; do
  # Extract the base filename without path
  base_name=$(basename "$input_file" .fasta)
  
  # Construct output file path
  output_file="$OUTPUT_DIR/output_${base_name}.fasta"
  
  # Run the Python script
  python ~/fastpace/test/fastpace_cmd.py --input_file "$input_file" --output_file "$output_file" --draw_logo
done

