import subprocess
import time
from typing import Literal
from pathlib import Path
import os

from ..seq import deduplicate_fasta_by_length
from ..logger import logger

cmd_makeblastdb_base = "makeblastdb -in {input_file} "\
                    "-out {outdir} "\
                    "-dbtype {dbtype} -parse_seqids -hash_index"

def makeblastdb(input_path:Path, output_path:Path, dbtype:Literal['nucl', 'prot']) -> Path:
    """
    Create BLAST database from a single file or all files in a directory.

    Args:
        input_path: Input fasta file or directory containing fasta files.
        output_path: BLAST database output directory.
        dbtype: Database type, 'nucl' or 'prot'.

    Returns:
        Path to the BLAST database.
    """
    input_path = Path(input_path)
    output_path = Path(output_path)
    os.makedirs(output_path, exist_ok=True)
    logger.info('makeblastdb: input=%s, output=%s, dbtype=%s', input_path, output_path, dbtype)
    if input_path.is_file():
        makeblastdb_singlefile(input_path, output_path, dbtype)
    elif input_path.is_dir():
        for subdir, _, files in os.walk(input_path):
            for file in files:
                if file.split('.')[-1] not in ['fasta', 'fa']:
                    continue
                makeblastdb_singlefile(os.path.join(subdir, file), output_path, dbtype)
    else:
        logger.error(f'Input path {input_path} does not exist')
        return
    if os.path.exists(output_path / 'tmp.fasta'):
        os.remove(output_path / 'tmp.fasta')
    return output_path

def makeblastdb_singlefile(input_file, output_path, dbtype) -> str:
    """
    Create BLAST database from a single file.

    Args:
        input_file: Input fasta file path.
        output_path: BLAST database output directory.
        dbtype: Database type.

    Returns:
        Path to the BLAST database.
    """
    file_name = os.path.basename(input_file).split('.')[0]
    deduplicate_fasta_by_length(output_path / 'tmp.fasta', input_file)
    cmd_makeblastdb = cmd_makeblastdb_base.format(
        input_file=output_path / 'tmp.fasta',
        outdir=output_path / file_name,
        dbtype=dbtype)
    logger.debug(f'makeblastdb command: {cmd_makeblastdb}')
    t0 = time.time()
    result = subprocess.run(cmd_makeblastdb, stderr=subprocess.PIPE, shell=True)
    elapsed = time.time() - t0
    if result.returncode != 0:
        logger.error('makeblastdb failed for %s (returncode=%s) after %.2fs: %s',
                     file_name, result.returncode, elapsed,
                     result.stderr.decode(errors='replace') if result.stderr else '(no stderr)')
        raise RuntimeError(f'makeblastdb failed for {file_name} (returncode={result.returncode})')
    if result.stderr:
        logger.warning('makeblastdb stderr for %s (returncode=0, %.2fs): %s',
                       file_name, elapsed, result.stderr.decode(errors='replace'))
    else:
        logger.info('makeblastdb done for %s in %.2fs', file_name, elapsed)
    return output_path / file_name

if __name__ == '__main__':
    import sys

    makeblastdb(Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3])
