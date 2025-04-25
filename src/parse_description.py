import re
import json
import argparse
import logging
from typing import Dict, List, Tuple
from utils import write_json

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

def parse_input(input_str: str) -> Tuple[List[str], List[int]]:
    """
    Parse input string to extract domains and fusion points.

    Args:
        input_str: Input string, e.g., "2CaM-BLA 41/197 Chimera".

    Returns:
        Tuple of (domains, fusion_points), e.g., (["CaM", "CaM", "BLA"], [41, 197]).

    Raises:
        ValueError: If input format is invalid.
    """
    try:
        # Remove "Chimera" suffix if present
        input_str = input_str.replace("Chimera", "").strip()
        
        # Split into domains and fusion points
        match = re.match(r"(\S+)\s+(\d+)/(\d+)", input_str)
        if not match:
            raise ValueError(f"Invalid input format: {input_str}")

        domains_part, fp1, fp2 = match.groups()
        fusion_points = [int(fp1), int(fp2)]

        # Parse domains (e.g., "2CaM-BLA" -> ["CaM", "CaM", "BLA"])
        domains = []
        for part in domains_part.split("-"):
            if part.startswith("2"):
                domains.extend([part[1:]] * 2)
            else:
                domains.append(part)

        logger.info(f"Parsed domains: {domains}, fusion points: {fusion_points}")
        return domains, fusion_points
    except Exception as e:
        logger.error(f"Error parsing input: {e}")
        raise ValueError(f"Failed to parse input: {e}")

def save_parsed_data(domains: List[str], fusion_points: List[int], output_path: str) -> None:
    """
    Save parsed data to a JSON file.

    Args:
        domains: List of domain names.
        fusion_points: List of fusion point indices.
        output_path: Path to output JSON file.
    """
    data = {"domains": domains, "fusion_points": fusion_points}
    write_json(data, output_path)
    logger.info(f"Saved parsed data to {output_path}")

def main(args: argparse.Namespace) -> None:
    """Main function to parse input and save results."""
    domains, fusion_points = parse_input(args.input)
    save_parsed_data(domains, fusion_points, args.output)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse chimeric protein description")
    parser.add_argument("--input", required=True, help="Input description (e.g., '2CaM-BLA 41/197 Chimera')")
    parser.add_argument("--output", required=True, help="Output JSON file path")
    args = parser.parse_args()

    try:
        main(args)
        logger.info("Parsing completed successfully")
    except Exception as e:
        logger.error(f"Parsing failed: {e}")
        raise