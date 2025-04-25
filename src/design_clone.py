import json
import argparse
import logging
from typing import Dict, List, Tuple
from Bio import PDB
from Bio.Seq import Seq
from Bio.SeqUtils import seq1
from Bio.Data import CodonTable
from utils import write_json

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

def load_parsed_data(json_path: str) -> Tuple[List[str], List[int]]:
    """
    Load domains and fusion points from parsed.json.

    Args:
        json_path: Path to parsed JSON file.

    Returns:
        Tuple of (domains, fusion_points).

    Raises:
        ValueError: If JSON is invalid or missing required fields.
    """
    try:
        with open(json_path) as f:
            data = json.load(f)
        domains = data.get("domains", [])
        fusion_points = data.get("fusion_points", [])
        if not domains or not fusion_points:
            raise ValueError("Missing 'domains' or 'fusion_points' in JSON")
        logger.info(f"Loaded domains: {domains}, fusion_points: {fusion_points}")
        return domains, fusion_points
    except Exception as e:
        logger.error(f"Error loading JSON: {e}")
        raise ValueError(f"Failed to load JSON: {e}")

def extract_chain_sequences(pdb_file: str) -> Dict[str, str]:
    """
    Extract protein sequences from each chain in a PDB file.

    Args:
        pdb_file: Path to PDB file.

    Returns:
        Dictionary mapping chain IDs to sequences.

    Raises:
        ValueError: If no valid sequences are extracted.
    """
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)
    chain_sequences = {}
    
    for model in structure:
        for chain in model:
            chain_id = chain.id
            sequence = []
            for residue in chain:
                resname = residue.get_resname()
                try:
                    aa = seq1(resname, undef_code="X")
                    sequence.append(aa)
                except KeyError:
                    logger.warning(f"Skipping non-standard residue {resname} in chain {chain_id}")
                    continue
            sequence_str = "".join(sequence)
            if sequence_str:
                chain_sequences[chain_id] = sequence_str
                logger.info(f"Extracted sequence for chain {chain_id}: {sequence_str[:20]}...")
    
    if not chain_sequences:
        raise ValueError("No valid protein sequences extracted")
    return chain_sequences

def construct_chimeric_sequence(chain_sequences: Dict[str, str], domains: List[str], fusion_points: List[int]) -> str:
    """
    Construct the chimeric protein sequence using fusion points.

    Args:
        chain_sequences: Dictionary of chain IDs to sequences.
        domains: List of domain names (e.g., ['CaM', 'CaM', 'BLA']).
        fusion_points: List of fusion point indices (e.g., [41, 197]).

    Returns:
        Chimeric protein sequence.

    Raises:
        ValueError: If sequences or fusion points are invalid.
    """
    try:
        if len(domains) != 3 or len(fusion_points) != 2:
            raise ValueError(f"Expected 3 domains and 2 fusion points, got {domains}, {fusion_points}")
        
        # Map domains to chains (A: CaM1, B: CaM2, C: BLA)
        chain_ids = ["A", "B", "C"]
        segments = []
        
        # CaM1: residues 1 to fusion_point[0]
        cam1_seq = chain_sequences.get("A", "")
        if len(cam1_seq) < fusion_points[0]:
            raise ValueError(f"CaM1 sequence too short: {len(cam1_seq)} < {fusion_points[0]}")
        segments.append(cam1_seq[:fusion_points[0]])
        
        # CaM2: residues fusion_point[0]+1 to fusion_point[1]
        cam2_seq = chain_sequences.get("B", "")
        cam2_length = fusion_points[1] - fusion_points[0]
        if len(cam2_seq) < cam2_length:
            raise ValueError(f"CaM2 sequence too short: {len(cam2_seq)} < {cam2_length}")
        segments.append(cam2_seq[:cam2_length])
        
        # BLA: residues fusion_point[1]+1 to end
        bla_seq = chain_sequences.get("C", "")
        segments.append(bla_seq)
        
        chimeric_seq = "".join(segments)
        if not chimeric_seq or "X" in chimeric_seq:
            raise ValueError("Invalid chimeric sequence generated")
        
        logger.info(f"Constructed chimeric sequence: {chimeric_seq[:20]}... (length: {len(chimeric_seq)})")
        return chimeric_seq
    except Exception as e:
        logger.error(f"Error constructing chimeric sequence: {e}")
        raise ValueError(f"Failed to construct chimeric sequence: {e}")

def back_translate_protein(protein_seq: str) -> str:
    """
    Back-translate protein sequence to DNA using a standard codon table.

    Args:
        protein_seq: Protein sequence in one-letter codes.

    Returns:
        DNA sequence.

    Raises:
        ValueError: If translation fails.
    """
    try:
        # Use standard genetic code (Table 1)
        codon_table = CodonTable.unambiguous_dna_by_id[1]
        aa_to_codon = {}
        for codon, aa in codon_table.forward_table.items():
            aa_to_codon.setdefault(aa, []).append(codon)
        
        dna_seq = []
        for aa in protein_seq:
            codons = aa_to_codon.get(aa)
            if not codons:
                raise ValueError(f"No codon for amino acid: {aa}")
            # Use the first codon (simplest approach)
            dna_seq.append(codons[0])
        
        dna_seq_str = "".join(dna_seq)
        logger.info(f"Back-translated DNA sequence: {dna_seq_str[:20]}... (length: {len(dna_seq_str)})")
        return dna_seq_str
    except Exception as e:
        logger.error(f"Error back-translating protein: {e}")
        raise ValueError(f"Failed to back-translate protein: {e}")

def design_primers(dna_seq: str, primer_length: int = 20) -> Dict[str, str]:
    """
    Design forward and reverse primers from DNA sequence.

    Args:
        dna_seq: DNA sequence.
        primer_length: Length of each primer.

    Returns:
        Dictionary with forward and reverse primers.
    """
    try:
        dna = Seq(dna_seq)
        forward = dna[:primer_length]
        reverse = dna[-primer_length:].reverse_complement()
        primers = {"forward": str(forward), "reverse": str(reverse)}
        logger.info(f"Designed primers: {primers}")
        return primers
    except Exception as e:
        logger.error(f"Error designing primers: {e}")
        raise ValueError(f"Failed to design primers: {e}")

def design_clone(pdb_file: str, json_path: str, vector: str = "pET-28a") -> Dict:
    """
    Design DNA sequence and primers for a chimeric protein.

    Args:
        pdb_file: Path to PDB file.
        json_path: Path to parsed JSON file.
        vector: Cloning vector name.

    Returns:
        Dictionary with protein sequence, DNA sequence, primers, and vector.

    Raises:
        ValueError: If design fails.
    """
    domains, fusion_points = load_parsed_data(json_path)
    chain_sequences = extract_chain_sequences(pdb_file)
    protein_seq = construct_chimeric_sequence(chain_sequences, domains, fusion_points)
    dna_seq = back_translate_protein(protein_seq)
    primers = design_primers(dna_seq)
    
    result = {
        "protein_sequence": protein_seq,
        "dna_sequence": dna_seq,
        "primers": primers,
        "vector": vector
    }
    return result

def main(args: argparse.Namespace) -> None:
    """Main function to design DNA sequence."""
    result = design_clone(args.input, args.input_json, vector="pET-28a")
    write_json(result, args.output)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Design DNA sequence for cloning")
    parser.add_argument("--input", required=True, help="Input PDB file")
    parser.add_argument("--input-json", required=True, help="Input JSON file with domains and fusion points")
    parser.add_argument("--output", required=True, help="Output JSON file")
    args = parser.parse_args()

    try:
        main(args)
        logger.info("Clone design completed successfully")
    except Exception as e:
        logger.error(f"Clone design failed: {e}")
        raise