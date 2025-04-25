import os
from typing import List, Dict

# Directory paths
DATA_DIR: str = "/app/data"
STRUCTURES_DIR: str = os.path.join(DATA_DIR, "structures")

# PDB mappings
DOMAIN_TO_PDB: Dict[str, str] = {
    "CaM": "1CLL",
    "BLA": "1ERQ"
}

# Domain colors for PyMOL visualization
DOMAIN_COLORS: List[str] = ["cyan", "green", "red"]  # CaM1, CaM2, BLA

# File extensions
PDB_EXT: str = ".pdb"
JSON_EXT: str = ".json"
PML_EXT: str = ".pml"