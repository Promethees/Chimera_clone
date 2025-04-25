import json
import os
import argparse

def write_json(data, output_path, indent=2):
    """
    Write data to a JSON file with error handling and logging.
    
    Args:
        data: Data to write (e.g., dict, list)
        output_path: Path to output JSON file
        indent: JSON indentation level (default: 2)
    
    Returns:
        None
    
    Raises:
        ValueError: If data is None or empty
        IOError: If file writing fails
    """
    try:
        if data is None or (isinstance(data, (list, dict)) and not data):
            raise ValueError(f"Cannot write empty or None data to {output_path}")
        
        # Ensure the output directory exists
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        
        print(f"Writing JSON to: {output_path}")
        with open(output_path, 'w') as f:
            json.dump(data, f, indent=indent)
        print(f"Successfully wrote JSON to: {output_path}")
    except Exception as e:
        print(f"Error writing JSON to {output_path}: {e}")
        raise

def default_parser(description="Processing...", arg1='--input', help1="Processing Input",
                                                arg2='--output', help2="Processing output"):
    parser = argparse.ArgumentParser(description)
    parser.add_argument(arg1, required=True, help="Protein description (e.g., '2CaM-BLA 41/197 Chimera')")
    parser.add_argument(arg2, required=True, help="Output JSON file")
    return parser.parse_args()