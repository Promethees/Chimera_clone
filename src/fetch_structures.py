import json
import os
import argparse
from Bio.PDB import PDBList
from utils import write_json

def fetch_structures(input_json, output_dir):
    """
    Fetch PDB files for domains listed in input JSON and save to output directory.
    
    Args:
        input_json (str): Path to input JSON file (e.g., parsed.json)
        output_dir (str): Directory to save PDB files
    
    Returns:
        dict: Metadata about downloaded PDB files
    """
    try:
        # Read input JSON
        print(f"Reading input JSON: {input_json}")
        with open(input_json) as f:
            data = json.load(f)
        domains = data.get("domains", [])
        if not domains:
            raise ValueError("No domains found in input JSON")

        # Map domains to PDB IDs
        pdb_map = {"CaM": "1CLL", "BLA": "1BTS"}
        pdbl = PDBList()
        pdb_files = []

        # Create output directory
        os.makedirs(output_dir, exist_ok=True)

        # Download PDB files
        for domain in domains:
            pdb_id = pdb_map.get(domain)
            if not pdb_id:
                raise ValueError(f"No PDB ID mapped for domain: {domain}")
            print(f"Downloading PDB {pdb_id} for domain {domain}")
            # Bio.PDB saves to <output_dir>/pdb<xxxx>.ent; we rename to <xxxx>.pdb
            pdbl.retrieve_pdb_file(pdb_id, pdir=output_dir, file_format="pdb")
            src_file = os.path.join(output_dir, f"pdb{pdb_id.lower()}.ent")
            dest_file = os.path.join(output_dir, f"{pdb_id}.pdb")
            if os.path.exists(src_file):
                os.rename(src_file, dest_file)
                pdb_files.append(dest_file)
                print(f"Saved PDB file: {dest_file}")
            else:
                raise IOError(f"Failed to download PDB {pdb_id}")

        # Return metadata
        result = {
            "domains": domains,
            "pdb_files": [os.path.basename(f) for f in pdb_files]
        }
        print(f"Fetched structures: {result}")
        return result
    except Exception as e:
        print(f"Error fetching structures: {e}")
        raise

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch PDB structures for domains")
    parser.add_argument('--input', required=True, help="Input JSON file (e.g., parsed.json)")
    parser.add_argument('--output', required=True, help="Output directory for PDB files")
    args = parser.parse_args()

    try:
        result = fetch_structures(args.input, args.output)
        # Optionally save metadata to JSON (not required by pipeline)
        metadata_path = os.path.join(args.output, "structures_metadata.json")
        write_json(result, metadata_path)
        print("Fetch structures completed")
    except Exception as e:
        print(f"Error: {e}")
        raise