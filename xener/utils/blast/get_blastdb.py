from pathlib import Path
import os

from ..logger import logger

def get_blastdb(dir:Path) -> dict[str, str]:
    """
    Get the list of BLAST databases in a directory.

    Args:
        dir: Directory containing BLAST database files.

    Returns:
        Dictionary mapping database names to file paths.
    """
    name2path = {}
    for subdir, _, files in os.walk(dir):
        for file in files:
            if file.split('.')[-1] not in ['pin','nin']:
                continue
            file_name = file.split('.')[0]
            name2path[file_name] = os.path.join(subdir, file_name)
    logger.debug(f'Found {len(name2path)} BLAST databases in {dir}')
    return name2path

if __name__ == '__main__':
    import sys
    name2path = get_blastdb(Path(sys.argv[1]))
    logger.info(f'BLAST databases: {name2path}')
