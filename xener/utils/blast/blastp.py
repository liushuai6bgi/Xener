import sys
import time
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
        logger.info(f'BLASTP starting: query={query_fasta}, db={db_path}, threads={num_threads}')
        logger.debug(f'BLASTP command: {cmd_blastp}')
        t0 = time.time()
        result = subprocess.run(cmd_blastp, stderr=subprocess.PIPE, shell=True)
        elapsed = time.time() - t0
        if result.returncode != 0:
            logger.error('BLASTP failed with returncode=%s after %.2fs: %s',
                         result.returncode, elapsed,
                         result.stderr.decode(errors='replace') if result.stderr else '(no stderr)')
            raise RuntimeError(f'BLASTP failed (returncode={result.returncode})')
        if result.stderr:
            logger.warning('BLASTP stderr (returncode=0, %.2fs): %s', elapsed, result.stderr.decode(errors='replace'))
        logger.info('BLASTP done in %.2fs: %s hits in %s', elapsed,
                    sum(1 for _ in open(tmp_output_file)) - 1 if os.path.exists(tmp_output_file) else 0,
                    tmp_output_file)
        blast_result = pd.read_csv(tmp_output_file, sep='\t', header=None)
        blast_result.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
        blast_result = blast_result.loc[blast_result.groupby(['qseqid', 'sseqid'])['bitscore'].idxmax()]
        logger.info('BLASTP best-hits per (qseqid,sseqid): %s rows', len(blast_result))
        os.remove(tmp_output_file)
        blast_result.to_csv(output_file, index=False)
    else:
        blast_result = pd.read_csv(output_file)
        logger.info(f'BLASTP cache hit: {output_file} ({len(blast_result)} rows, skipping alignment)')

    return blast_result

if __name__ == '__main__':
    res = blastp(Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]))
    logger.info(f'BLASTP result: {res}')
