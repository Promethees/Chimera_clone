import os
import argparse
import pymol.cmd as cmd
from utils import write_json

def model_chimera(pdb_files, fusion_points, output_pdb):
    """
    Align PDB structures and save a chimeric protein structure.
    
    Args:
        pdb_files (list): List of PDB file paths
        fusion_points (list): List of fusion point indices
        output_pdb (str): Path to output PDB file
    
    Returns:
        dict: Metadata about the modeled structure
    """
    try:
        print("Initializing PyMOL")
        cmd.reinitialize()
        if not pdb_files:
            raise ValueError("No PDB files provided")

        # Load PDB files
        for i, pdb in enumerate(pdb_files):
            print(f"Loading PDB {pdb} as domain_{i}")
            cmd.load(pdb, f"domain_{i}")

        # Align domains if multiple
        if len(pdb_files) > 1:
            print("Aligning domains")
            cmd.align("domain_0", "domain_1")

        # Select all atoms for saving
        cmd.select("all")
        print(f"Saving chimeric structure to {output_pdb}")
        cmd.save(output_pdb)

        # Verify output file
        if not os.path.exists(output_pdb) or os.path.getsize(output_pdb) == 0:
            raise IOError(f"Failed to save PDB file: {output_pdb}")

        # Return metadata
        result = {
            "pdb_files": [os.path.basename(f) for f in pdb_files],
            "fusion_points": fusion_points,
            "output_pdb": os.path.basename(output_pdb)
        }
        print(f"Modeled chimera: {result}")
        return result
    except Exception as e:
        print(f"Error modeling chimera: {e}")
        raise

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Model chimeric protein structure")
    parser.add_argument('--input', required=True, help="Directory with PDB files")
    parser.add_argument('--fusion-points', nargs=2, type=int, required=True, help="Fusion point indices")
    parser.add_argument('--output', required=True, help="Output PDB file")
    args = parser.parse_args()

    try:
        print(f"Reading PDB files from: {args.input}")
        pdb_files = [os.path.join(args.input, f) for f in os.listdir(args.input) if f.endswith('.pdb')]
        if not pdb_files:
            raise ValueError(f"No PDB files found in {args.input}")
        print(f"Found PDB files: {pdb_files}")
        result = model_chimera(pdb_files, args.fusion_points, args.output)
        # Save metadata to JSON
        metadata_path = os.path.splitext(args.output)[0] + "_metadata.json"
        write_json(result, metadata_path)
        print("Modeling completed")
    except Exception as e:
        print(f"Error: {e}")
        raise