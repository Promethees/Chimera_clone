import os
import json
import argparse
import pymol.cmd as cmd
from utils import write_json

def model_chimera(pdb_files, fusion_points, output_pdb):
    """
    Align PDB structures, assign distinct colors and chain IDs to each domain, and save a chimeric protein structure.
    
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

        # Assign domain names and track CaM instances
        domain_names = []
        cam_count = 0
        for pdb in pdb_files:
            if "1CLL" in os.path.basename(pdb):
                cam_count += 1
                domain_names.append(f"CaM{cam_count}")
            else:
                domain_names.append("BLA")

        # Load PDB files and assign chain IDs
        for i, (pdb, domain) in enumerate(zip(pdb_files, domain_names)):
            print(f"Loading PDB {pdb} as domain_{i} ({domain})")
            cmd.load(pdb, f"domain_{i}")
            # Assign unique chain ID (A, B, C, ...)
            cmd.alter(f"domain_{i}", f"chain='{chr(65+i)}'")

        # Align domains if multiple
        if len(pdb_files) > 1:
            print("Aligning domains")
            for i in range(1, len(pdb_files)):
                cmd.align(f"domain_{i}", "domain_0")

        # Assign distinct colors to each domain
        colors = ["cyan", "green", "red"]  # CaM1: cyan, CaM2: green, BLA: red
        for i, domain in enumerate(domain_names):
            cmd.color(colors[i % len(colors)], f"domain_{i}")
            print(f"Colored domain_{i} ({domain}) as {colors[i % len(colors)]}")

        # Select all atoms for saving
        cmd.select("all")
        print(f"Saving chimeric structure to {output_pdb}")
        cmd.save(output_pdb, "all")

        # Verify output file
        if not os.path.exists(output_pdb) or os.path.getsize(output_pdb) == 0:
            raise IOError(f"Failed to save PDB file: {output_pdb}")

        # Generate PyMOL script for enhanced visualization
        pml_file = os.path.splitext(output_pdb)[0] + "_color.pml"
        with open(pml_file, "w") as f:
            f.write(f"load {os.path.basename(output_pdb)}\n")
            for i, (domain, color) in enumerate(zip(domain_names, colors)):
                chain = chr(65+i)
                f.write(f"select {domain}, chain {chain}\n")
                f.write(f"color {color}, {domain}\n")
                # Add label at the center of mass of the domain
                f.write(f"center {domain}\n")
                f.write(f"label {domain} and name CA, '{domain}'\n")
            f.write("show cartoon\n")
            f.write("set cartoon_fancy_helices, 1\n")
            f.write("set cartoon_transparency, 0.2\n")
            f.write("show surface\n")
            f.write("set transparency, 0.7\n")
            f.write("bg_color white\n")
            f.write("zoom\n")
        print(f"Saved PyMOL coloring script to {pml_file}")

        # Return metadata
        result = {
            "pdb_files": [os.path.basename(f) for f in pdb_files],
            "domains": domain_names,
            "fusion_points": fusion_points,
            "output_pdb": os.path.basename(output_pdb),
            "color_script": os.path.basename(pml_file)
        }
        print(f"Modeled chimera: {result}")
        return result
    except Exception as e:
        print(f"Error modeling chimera: {e}")
        raise

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Model chimeric protein structure")
    parser.add_argument('--input', required=True, help="Directory with PDB files")
    parser.add_argument('--input-json', required=True, help="Input JSON file with fusion points (e.g., parsed.json)")
    parser.add_argument('--output', required=True, help="Output PDB file")
    args = parser.parse_args()

    try:
        print(f"Reading PDB files from: {args.input}")
        pdb_files = [os.path.join(args.input, f) for f in os.listdir(args.input) if f.endswith('.pdb')]
        if not pdb_files:
            raise ValueError(f"No PDB files found in {args.input}")
        print(f"Found PDB files: {pdb_files}")

        print(f"Reading fusion points from: {args.input_json}")
        with open(args.input_json) as f:
            data = json.load(f)
        fusion_points = data.get("fusion_points", [])
        if not fusion_points or len(fusion_points) != 2:
            raise ValueError(f"Invalid or missing fusion_points in {args.input_json}")

        result = model_chimera(pdb_files, fusion_points, args.output)
        # Save metadata to JSON
        metadata_path = os.path.splitext(args.output)[0] + "_metadata.json"
        write_json(result, metadata_path)
        print("Modeling completed")
    except Exception as e:
        print(f"Error: {e}")
        raise