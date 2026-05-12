from Bio import SeqIO
from Bio.Seq import Seq

from .logger import logger

def deduplicate_fasta_by_length(outfile:str, fasta_file:str) -> str:
    """
    Deduplicate a FASTA file, keeping only the longest sequence for each sequence ID.

    Args:
        outfile: Output file path.
        fasta_file: Input FASTA file path.

    Returns:
        Output file path.
    """
    all_sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        gene_id = record.id
        if gene_id in all_sequences and len(all_sequences[gene_id]) > len(record.seq):
            continue
        all_sequences[gene_id] = record.seq
    with open(outfile, "w") as output_handle:
        for gene_id, sequence in all_sequences.items():
            SeqIO.write(SeqIO.SeqRecord(sequence, id=gene_id, description=""), output_handle, "fasta")
    return outfile

def extract_fasta_by_name(outdir, name_list:list[str], fasta_file) -> str:
    """
    Extract sequences with specified names from a FASTA file.

    Args:
        outdir: Output directory.
        name_list: List of sequence names to extract.
        fasta_file: FASTA file path.

    Returns:
        Output file path.
    """
    name2original_name = {}
    extracted_sequences = {}
    name_set = set()
    for original_name in name_list:
        name = original_name
        name = name.replace('-', '_')
        name2original_name[name] = original_name
        name_set.add(name)
    for record in SeqIO.parse(fasta_file, "fasta"):
        gene_id = record.id
        gene_id = gene_id.replace('-', '_')
        if gene_id in name_set:
            gene_id = name2original_name[gene_id]
            extracted_sequences[gene_id] = record.seq
    query_fasta = f"{outdir}/query.fasta"
    with open(query_fasta, "w") as output_handle:
        for gene_id, sequence in extracted_sequences.items():
            SeqIO.write(SeqIO.SeqRecord(sequence, id=gene_id, description=""), output_handle, "fasta")
    logger.info(f"Extracted {len(extracted_sequences)}/{len(name_list)} sequences from {fasta_file}")
    return query_fasta

def translate2protein(fasta_file:str) -> str:
    """
    Translate nucleic acid sequences to protein sequences.

    Args:
        fasta_file: Input nucleic acid FASTA file path.

    Returns:
        Output protein FASTA file path.
    """
    outfile = fasta_file.split('.')[0] + '_prot.fasta'
    all_sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        gene_id = record.id
        if gene_id in all_sequences and len(all_sequences[gene_id]) > len(record.seq):
            continue
        all_sequences[gene_id] = Seq(record.seq).translate()
    with open(outfile, "w") as output_handle:
        for gene_id, sequence in all_sequences.items():
            SeqIO.write(SeqIO.SeqRecord(sequence, id=gene_id, description=""), output_handle, "fasta")
    return outfile
