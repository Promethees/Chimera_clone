#!/bin/bash
# Get the directory containing this script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Project root is one directory up from scripts/
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
SRC_DIR="$PROJECT_ROOT/src"
DATA_DIR="$PROJECT_ROOT/data"

# Initialize Conda
source /opt/conda/etc/profile.d/conda.sh
# Activate conda environment
conda activate chimera_clone || { echo "Failed to activate chimera_clone"; exit 1; }

set -e  # Exit on error
echo "Starting pipeline..."
INPUT=${1:-"2CaM-BLA 41/197 Chimera"}
echo "Processing input: $INPUT"

# Run Python scripts with absolute paths, logging output
echo "Running parse_description.py..."
python "$SRC_DIR/parse_description.py" --input "$INPUT" --output "$DATA_DIR/parsed.json" 2>&1 | tee "$DATA_DIR/parse_description.log" && echo "Parsed description" || { echo "Parsing failed"; exit 1; }

echo "Running fetch_structures.py..."
python "$SRC_DIR/fetch_structures.py" --input "$DATA_DIR/parsed.json" --output "$DATA_DIR/structures" 2>&1 | tee "$DATA_DIR/fetch_structures.log" && echo "Fetched structures" || { echo "Structure fetch failed"; exit 1; }

echo "Running model_chimera.py..."
python "$SRC_DIR/model_chimera.py" --input "$DATA_DIR/structures" --fusion-points 41 197 --output "$DATA_DIR/chimera.pdb" 2>&1 | tee "$DATA_DIR/model_chimera.log" && echo "Modeled chimera" || { echo "Modeling failed"; exit 1; }

echo "Running design_clone.py..."
python "$SRC_DIR/design_clone.py" --input "$DATA_DIR/chimera.pdb" --output "$DATA_DIR/clone_design.json" 2>&1 | tee "$DATA_DIR/design_clone.json" && echo "Designed clone" || { echo "Clone design failed"; exit 1; }

echo "Pipeline completed"