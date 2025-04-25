import json
import logging
from typing import Any
import os 

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

def write_json(data: Any, output_path: str) -> None:
    """
    Write data to a JSON file.

    Args:
        data: Data to write (e.g., dict, list).
        output_path: Path to output JSON file.

    Raises:
        IOError: If writing fails.
    """
    try:
        with open(output_path, "w") as f:
            json.dump(data, f, indent=4)
        if not os.path.exists(output_path) or os.path.getsize(output_path) == 0:
            raise IOError(f"Failed to write JSON to {output_path}")
        logger.info(f"Wrote JSON to {output_path}")
    except Exception as e:
        logger.error(f"Error writing JSON: {e}")
        raise IOError(f"Failed to write JSON: {e}")