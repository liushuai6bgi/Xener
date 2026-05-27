import sys
import subprocess
from typing import Literal
from pathlib import Path
import os

import pandas as pd

from ..logger import logger

cmd_blastp_base = "blastp -query {query_fasta} -db {db_path} -out {output_file} " +\
                    "{other_args} -task blastp-fast -outfmt 6 " +\
                    "-mt_mode 1 -num_threads {num_threads} "

def blastp(query_fasta:Path, db_path:Path, output_file:Path, num_threads:int=None, **kwargs) -> pd.DataFrame:
    """
    Execute BLASTP alignment.

    Args:
        query_fasta: Query fasta file path.
        db_path: BLAST database path.
        output_file: BLASTP result output file path.
        num_threads: Number of threads; if None, automatically set.
        kwargs: Additional BLASTP parameters, default evalue=1e-5.

    Returns:
        Processed BLASTP result DataFrame.
    """
    tmp_output_file = output_file.with_suffix('.csv')
    if not os.path.exists(output_file):
        if num_threads is None:
            num_threads = int(os.cpu_count() * 0.7)
        kwargs.setdefault('evalue', 1e-5)
        cmd_blastp = cmd_blastp_base.format(
            query_fasta=query_fasta, db_path=db_path,
            output_file=tmp_output_file, num_threads=num_threads,
            other_args=' '.join([f"-{k} {v}" for k, v in kwargs.items()])
        )
        logger.info(f'BLASTP command: {cmd_blastp}')
        result = subprocess.run(cmd_blastp, stderr=subprocess.PIPE, shell=True)
        if result.stderr:
            logger.error(f'BLASTP error: {result.stderr.decode()}')
        blast_result = pd.read_csv(tmp_output_file, sep='\t', header=None)
        blast_result.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
        blast_result = blast_result.loc[blast_result.groupby(['qseqid', 'sseqid'])['bitscore'].idxmax()]
        os.remove(tmp_output_file)
        blast_result.to_csv(output_file, index=False)
    else:
        blast_result = pd.read_csv(output_file)
        logger.info(f'{output_file} exists, loading from cache')

    return blast_result

if __name__ == '__main__':
    res = blastp(Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]))
    logger.info(f'BLASTP result: {res}')
