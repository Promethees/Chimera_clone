import argparse
import re
from utils import write_json, default_parser

def parse_protein_description(description):
    try:
        pattern = r"(\d*)([A-Za-z]+)-([A-Za-z]+)\s(\d+)/(\d+)"
        match = re.match(pattern, description)
        if match:
            count, domain1, domain2, point1, point2 = match.groups()
            count = int(count) if count else 1
            domains = [domain1] * count + [domain2]
            result = {
                "domains": domains,
                "fusion_points": [int(point1), int(point2)]
            }
            print(f"Parsed result: {result}")
            return result
        raise ValueError(f"Invalid description format: {description}")
    except Exception as e:
        print(f"Error parsing description: {e}")
        raise

if __name__ == "__main__":
    args = default_parser(description="Parse protein description", arg1='--input', help1="Protein description (e.g., '2CaM-BLA 41/197 Chimera')",
                                                                arg2='--output', help2="Output JSON file")
    try:
        print(f"Parsing input: {args.input}")
        result = parse_protein_description(args.input)
        write_json(result, args.output)
        print("Parse completed")
    except Exception as e:
        print(f"Error: {e}")
        raise