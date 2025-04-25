import os
import json
import argparse
import logging
import urllib.request
from typing import List, Dict
from config import DOMAIN_TO_PDB, STRUCTURES_DIR, PDB_EXT

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

def load_parsed_data(input_path: str) -> Dict:
    """
    Load parsed data from JSON file.

    Args:
        input_path: Path to input JSON file.

    Returns:
        Dictionary containing parsed data.

    Raises:
        ValueError: If JSON file is invalid or missing required fields.
    """
    try:
        with open(input_path) as f:
            data = json.load(f)
        if "domains" not in data:
            raise ValueError("Missing 'domains' in input JSON")
        logger.info(f"Loaded parsed data from {input_path}")
        return data
    except Exception as e:
        logger.error(f"Error loading JSON: {e}")
        raise ValueError(f"Failed to load JSON: {e}")

def fetch_pdb_files(domains: List[str], output_dir: str) -> List[str]:
    """
    Fetch PDB files for given domains and save to output directory using urllib.

    Args:
        domains: List of domain names.
        output_dir: Directory to save PDB files.

    Returns:
        List of paths to saved PDB files.

    Raises:
        ValueError: If PDB fetch fails or domain is unknown.
    """
    try:
        os.makedirs(output_dir, exist_ok=True)
        pdb_files = []
        for domain in domains:
            pdb_id = DOMAIN_TO_PDB.get(domain)
            if not pdb_id:
                raise ValueError(f"Unknown domain: {domain}")
            output_path = os.path.join(output_dir, f"{pdb_id}{PDB_EXT}")
            url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
            urllib.request.urlretrieve(url, output_path)
            if not os.path.exists(output_path) or os.path.getsize(output_path) == 0:
                raise ValueError(f"Failed to save PDB file: {output_path}")
            pdb_files.append(output_path)
            logger.info(f"Fetched PDB {pdb_id} for {domain} to {output_path}")
        return pdb_files
    except Exception as e:
        logger.error(f"Error fetching PDB files: {e}")
        raise ValueError(f"Failed to fetch PDB files: {e}")

def main(args: argparse.Namespace) -> None:
    """Main function to fetch PDB structures."""
    parsed_data = load_parsed_data(args.input)
    fetch_pdb_files(parsed_data["domains"], args.output)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch PDB structures for domains")
    parser.add_argument("--input", required=True, help="Input JSON file (e.g., parsed.json)")
    parser.add_argument("--output", required=True, help="Output directory for PDB files")
    args = parser.parse_args()

    try:
        main(args)
        logger.info("Structure fetching completed successfully")
    except Exception as e:
        logger.error(f"Structure fetching failed: {e}")
        raise