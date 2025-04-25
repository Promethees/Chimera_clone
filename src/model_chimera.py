import os
import json
import argparse
import logging
from typing import List, Dict, Tuple
import pymol.cmd as cmd
from utils import write_json
from config import DOMAIN_TO_PDB, DOMAIN_COLORS, PDB_EXT, PML_EXT

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

def load_pdb_files(input_dir: str) -> List[str]:
    """
    Load PDB files from the input directory.

    Args:
        input_dir: Directory containing PDB files.

    Returns:
        List of PDB file paths.

    Raises:
        ValueError: If no PDB files are found.
    """
    pdb_files = [os.path.join(input_dir, f) for f in os.listdir(input_dir) if f.endswith(PDB_EXT)]
    if not pdb_files:
        raise ValueError(f"No PDB files found in {input_dir}")
    logger.info(f"Found PDB files: {pdb_files}")
    return sorted(pdb_files)

def assign_domain_names(pdb_files: List[str]) -> List[str]:
    """
    Assign domain names based on PDB file names.

    Args:
        pdb_files: List of PDB file paths.

    Returns:
        List of domain names (e.g., ['CaM1', 'CaM2', 'BLA']).
    """
    domain_names = []
    cam_count = 0
    for pdb in pdb_files:
        pdb_name = os.path.basename(pdb)
        for domain, pdb_id in DOMAIN_TO_PDB.items():
            if pdb_id in pdb_name:
                if domain == "CaM":
                    cam_count += 1
                    domain_names.append(f"CaM{cam_count}")
                else:
                    domain_names.append(domain)
                break
    logger.info(f"Assigned domain names: {domain_names}")
    return domain_names

def load_and_align_domains(pdb_files: List[str], domain_names: List[str]) -> None:
    """
    Load PDB files into PyMOL, assign chain IDs, align domains, and color them.

    Args:
        pdb_files: List of PDB file paths.
        domain_names: List of domain names.
    """
    cmd.reinitialize()
    for i, (pdb, domain) in enumerate(zip(pdb_files, domain_names)):
        cmd.load(pdb, f"domain_{i}")
        cmd.alter(f"domain_{i}", f"chain='{chr(65+i)}'")
        color = DOMAIN_COLORS[i % len(DOMAIN_COLORS)]
        cmd.color(color, f"domain_{i}")
        logger.info(f"Loaded {pdb} as domain_{i} ({domain}), chain {chr(65+i)}, colored {color}")
    
    if len(pdb_files) > 1:
        for i in range(1, len(pdb_files)):
            cmd.align(f"domain_{i}", "domain_0")
        logger.info("Aligned domains")

def save_chimera_structure(output_pdb: str) -> None:
    """
    Save the chimeric structure to a PDB file.

    Args:
        output_pdb: Path to output PDB file.

    Raises:
        IOError: If saving fails.
    """
    cmd.select("all")
    cmd.save(output_pdb, "all")
    if not os.path.exists(output_pdb) or os.path.getsize(output_pdb) == 0:
        raise IOError(f"Failed to save PDB file: {output_pdb}")
    logger.info(f"Saved chimeric structure to {output_pdb}")

def generate_pymol_script(output_pdb: str, domain_names: List[str]) -> str:
    """
    Generate a PyMOL script for coloring and labeling domains.

    Args:
        output_pdb: Path to PDB file.
        domain_names: List of domain names.

    Returns:
        Path to the generated PyMOL script.
    """
    pml_file = os.path.splitext(output_pdb)[0] + PML_EXT
    with open(pml_file, "w") as f:
        f.write(f"load {os.path.basename(output_pdb)}\n")
        for i, domain in enumerate(domain_names):
            chain = chr(65+i)
            color = DOMAIN_COLORS[i % len(DOMAIN_COLORS)]
            f.write(f"select {domain}, chain {chain}\n")
            f.write(f"color {color}, {domain}\n")
            f.write(f"center {domain}\n")
            f.write(f"label {domain} and name CA, '{domain}'\n")
        f.write("show cartoon\n")
        f.write("set cartoon_fancy_helices, 1\n")
        f.write("set cartoon_transparency, 0.2\n")
        f.write("show surface\n")
        f.write("set transparency, 0.7\n")
        f.write("bg_color white\n")
        f.write("zoom\n")
    logger.info(f"Saved PyMOL script to {pml_file}")
    return pml_file

def model_chimera(pdb_files: List[str], fusion_points: List[int], output_pdb: str) -> Dict:
    """
    Model a chimeric protein structure using PyMOL.

    Args:
        pdb_files: List of PDB file paths.
        fusion_points: List of fusion point indices.
        output_pdb: Path to output PDB file.

    Returns:
        Metadata dictionary.
    """
    domain_names = assign_domain_names(pdb_files)
    load_and_align_domains(pdb_files, domain_names)
    save_chimera_structure(output_pdb)
    pml_file = generate_pymol_script(output_pdb, domain_names)
    
    return {
        "pdb_files": [os.path.basename(f) for f in pdb_files],
        "domains": domain_names,
        "fusion_points": fusion_points,
        "output_pdb": os.path.basename(output_pdb),
        "color_script": os.path.basename(pml_file)
    }

def load_fusion_points(input_json: str) -> List[int]:
    """
    Load fusion points from a JSON file.

    Args:
        input_json: Path to input JSON file.

    Returns:
        List of fusion points.

    Raises:
        ValueError: If fusion points are invalid.
    """
    with open(input_json) as f:
        data = json.load(f)
    fusion_points = data.get("fusion_points", [])
    if not fusion_points or len(fusion_points) != 2:
        raise ValueError(f"Invalid or missing fusion_points in {input_json}")
    logger.info(f"Loaded fusion points: {fusion_points}")
    return fusion_points

def main(args: argparse.Namespace) -> None:
    """Main function to model chimeric structure."""
    pdb_files = load_pdb_files(args.input)
    fusion_points = load_fusion_points(args.input_json)
    result = model_chimera(pdb_files, fusion_points, args.output)
    metadata_path = os.path.splitext(args.output)[0] + "_metadata.json"
    write_json(result, metadata_path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Model chimeric protein structure")
    parser.add_argument("--input", required=True, help="Directory with PDB files")
    parser.add_argument("--input-json", required=True, help="Input JSON file with fusion points")
    parser.add_argument("--output", required=True, help="Output PDB file")
    args = parser.parse_args()

    try:
        main(args)
        logger.info("Modeling completed successfully")
    except Exception as e:
        logger.error(f"Modeling failed: {e}")
        raise