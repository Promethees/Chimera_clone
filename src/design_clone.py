import argparse
from Bio.Seq import Seq
from Bio import PDB
from Bio.SeqUtils import seq1
from utils import write_json

def get_protein_sequence(pdb_file):
    """
    Extract protein sequence from a PDB file.
    
    Args:
        pdb_file (str): Path to PDB file
    
    Returns:
        str: Protein sequence (1-letter codes)
    """
    try:
        print(f"Reading PDB file: {pdb_file}")
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure('protein', pdb_file)
        sequence = ""
        for model in structure:
            for chain in model:
                for residue in chain:
                    resname = residue.get_resname()
                    try:
                        aa = seq1(resname)
                        sequence += aa
                    except KeyError:
                        print(f"Skipping non-standard residue: {resname}")
                        continue
        if not sequence:
            raise ValueError("No valid protein sequence extracted from PDB")
        print(f"Extracted sequence: {sequence}")
        return sequence
    except Exception as e:
        print(f"Error extracting sequence: {e}")
        raise

def design_clone(protein_seq, vector="pET-28a"):
    """
    Design DNA sequence and primers from protein sequence.
    
    Args:
        protein_seq (str): Protein sequence
        vector (str): Cloning vector (default: pET-28a)
    
    Returns:
        dict: DNA sequence and primers
    """
    try:
        print("Designing clone")
        dna_seq = Seq(protein_seq).back_transcribe()
        forward_primer = dna_seq[:20]
        reverse_primer = dna_seq[-20:].reverse_complement()
        result = {
            "dna_sequence": str(dna_seq),
            "primers": {
                "forward": str(forward_primer),
                "reverse": str(reverse_primer)
            },
            "vector": vector
        }
        print(f"Designed clone: {result}")
        return result
    except Exception as e:
        print(f"Error designing clone: {e}")
        raise

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Design DNA sequence for cloning")
    parser.add_argument('--input', required=True, help="Input PDB file")
    parser.add_argument('--output', required=True, help="Output JSON file")
    args = parser.parse_args()

    try:
        protein_seq = get_protein_sequence(args.input)
        result = design_clone(protein_seq)
        write_json(result, args.output)
        print("Clone design completed")
    except Exception as e:
        print(f"Error: {e}")
        raise