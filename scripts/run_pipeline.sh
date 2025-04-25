#!/bin/bash
# Get the directory containing this script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Project root is one directory up from scripts/
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
SRC_DIR="$PROJECT_ROOT/src"
DATA_DIR="/app/data"  # Explicitly set to mounted directory

# Clear DATA_DIR contents
echo "Clearing $DATA_DIR..."
rm -rf "$DATA_DIR"/* || { echo "Failed to clear $DATA_DIR"; exit 1; }
mkdir -p "$DATA_DIR" || { echo "Failed to create $DATA_DIR"; exit 1; }
chmod -R u+rwx "$DATA_DIR" || { echo "Failed to set permissions on $DATA_DIR"; exit 1; }
echo "Cleared and initialized $DATA_DIR"

# Verify DATA_DIR is writable
if [ ! -d "$DATA_DIR" ]; then
    echo "Error: $DATA_DIR does not exist"
    exit 1
fi
if [ ! -w "$DATA_DIR" ]; then
    echo "Error: $DATA_DIR is not writable"
    exit 1
fi

# Initialize Conda
source /opt/conda/etc/profile.d/conda.sh
conda activate chimera_clone || { echo "Failed to activate chimera_clone"; exit 1; }

set -e  # Exit on error
echo "Starting pipeline..."
INPUT=${1:-"2CaM-BLA 41/197 Chimera"}
echo "Processing input: $INPUT"

# Run parse_description.py
echo "Running parse_description.py with command:"
echo "python $SRC_DIR/parse_description.py --input '$INPUT' --output $DATA_DIR/parsed.json"
python "$SRC_DIR/parse_description.py" --input "$INPUT" --output "$DATA_DIR/parsed.json" 2>&1 | tee "$DATA_DIR/parse_description.log" && echo "Parsed description" || { echo "Parsing failed"; exit 1; }

# Run fetch_structures.py
echo "Running fetch_structures.py..."
echo "python $SRC_DIR/fetch_structures.py --input $DATA_DIR/parsed.json --output $DATA_DIR/structures"
python "$SRC_DIR/fetch_structures.py" --input "$DATA_DIR/parsed.json" --output "$DATA_DIR/structures" 2>&1 | tee "$DATA_DIR/fetch_structures.log" && echo "Fetched structures" || { echo "Structure fetch failed"; exit 1; }

# Run model_chimera.py
echo "Running model_chimera.py..."
echo "python $SRC_DIR/model_chimera.py --input $DATA_DIR/structures --input-json $DATA_DIR/parsed.json --output $DATA_DIR/chimera.pdb"
python "$SRC_DIR/model_chimera.py" --input "$DATA_DIR/structures" --input-json "$DATA_DIR/parsed.json" --output "$DATA_DIR/chimera.pdb" 2>&1 | tee "$DATA_DIR/model_chimera.log" && echo "Modeled chimera" || { echo "Modeling failed"; exit 1; }

# Run design_clone.py
echo "Running design_clone.py..."
echo "python $SRC_DIR/design_clone.py --input $DATA_DIR/chimera.pdb --output $DATA_DIR/clone_design.json"
python "$SRC_DIR/design_clone.py" --input "$DATA_DIR/chimera.pdb" --output "$DATA_DIR/clone_design.json" 2>&1 | tee "$DATA_DIR/design_clone.log" && echo "Designed clone" || { echo "Clone design failed"; exit 1; }

echo "Pipeline completed"