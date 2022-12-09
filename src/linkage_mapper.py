#!/usr/bin/env python3
__description__ =\
"""
Purpose: Identify genomic positions of restriction map relative to assembled genome and create vcf.
"""
__author__ = "Erick Samera"
__version__ = "0.0.1"
__comments__ = "stable;"
# --------------------------------------------------
from argparse import (
    Namespace,
    ArgumentParser,
    RawTextHelpFormatter)
from pathlib import Path
# --------------------------------------------------
from Bio import SeqIO
from itertools import combinations
from multiprocessing import Pool
import pickle
import tempfile
import subprocess
import json
import collections
# --------------------------------------------------
def get_args() -> Namespace:
    """ Get command-line arguments """

    parser = ArgumentParser(
        description=__description__,
        epilog=f"v{__version__} : {__author__} | {__comments__}",
        formatter_class=RawTextHelpFormatter)
    parser.add_argument(
        'input_path',
        type=Path,
        help="path of input genomic file (.fna)")
    parser.add_argument(
        'csv_path',
        type=Path,
        help="path of input genomic file (.fna)")

    args = parser.parse_args()

    # parser errors and processing
    # --------------------------------------------------

    return args
# --------------------------------------------------
def _translate_snp_position(_allele_position: int, _match_position: int): return _match_position + _allele_position
def _identify_markers(_marker: dict, _chrom: str, _only_top: bool = False) -> list:
    """
    """
    if _marker['sequence'] == "#N/A": return None

    temp_fasta = tempfile.NamedTemporaryFile()
    with open(temp_fasta.name, 'w') as temp_fasta_file:
        temp_fasta_file.write("> temporary query\n")
        temp_fasta_file.write(f"{_marker['sequence']}")

    subprocess.run([
        "blastn",
        "-db", "data/Rubus_idaeus_JoanJ.fna",
        "-outfmt", "15",
        "-query", f"{temp_fasta.name}",
        "-out", f"blast_output.json"])
    
    potential_hits: list = []
    with open('blast_output.json', 'r') as json_file:
        blast_object = json.load(json_file)
        blast_hits = (blast_object['BlastOutput2'][0]['report']['results']['search']['hits'])
        for hit in blast_hits:
            if _chrom == hit['description'][0]['id']:
                for i in hit['hsps']:
                    potential_hits.append(i)
        if not potential_hits:
            for hit in blast_hits:
                if 'scaffold' in hit['description'][0]['id']:
                    for i in hit['hsps']:
                        potential_hits.append(i)
    
    potential_hits_sorted = [potential_hit for potential_hit in sorted(potential_hits, key=lambda x: x['evalue'])]
    if len(potential_hits_sorted) > 1 and not _only_top:
        potential_hits_filtered = [hit for hit in potential_hits_sorted if hit['evalue'] > 100]
    else:
        potential_hits_filtered = potential_hits_sorted
    
    marker_dict_list: list = []
    for hit in potential_hits_filtered:
            for i_snp_position, snp_position in enumerate(_marker['SNP_positions']):
                try:
                    processed_marker_name = f"{_marker['marker_name']}_{i_snp_position}" if len(_marker['SNP_positions']) > 1 else _marker['marker_name']
                    print('\t'.join([
                        f'{_chrom}', 
                        f"{_translate_snp_position(int(snp_position), int(hit['hit_from']))}", 
                        f"{processed_marker_name}", 
                        f"{hit['hseq'][snp_position]}", 
                        f"{','.join(_marker['alleles'][i_snp_position])}",
                        f".",
                        f"CHECK",
                        f"."
                        ]))
                except IndexError:
                    pass
    
def _identify_chromosome(_linkage_group_dict: dict) -> str:
    """
    """

    chromosome_candidates: list = []

    for marker in _linkage_group_dict.values():

        if marker['sequence'] == "#N/A": continue

        temp_fasta = tempfile.NamedTemporaryFile()
        with open(temp_fasta.name, 'w') as temp_fasta_file:
            temp_fasta_file.write("> temporary query\n")
            temp_fasta_file.write(f"{marker['sequence']}")

        subprocess.run([
            "blastn",
            "-db", "data/Rubus_idaeus_JoanJ.fna",
            "-outfmt", "15",
            "-query", f"{temp_fasta.name}",
            "-out", f"blast_output.json"])
        
        with open('blast_output.json', 'r') as json_file:
            blast_object = json.load(json_file)
            blast_hits = (blast_object['BlastOutput2'][0]['report']['results']['search']['hits'])
            matches = [str(hit['description'][0]['id']) for hit in blast_hits]

            chromosome_candidates.append(matches)
    
    chromosome_candidates = sum(chromosome_candidates, [])

    return collections.Counter(chromosome_candidates).most_common(1)[0][0]
def _parse_linkage_map_csv(_csv_path: Path) -> dict:
    """
    """

    linkage_map_dict: dict = {}

    with open(_csv_path) as csv_file:
        for line in csv_file.readlines()[1:]:
            line_info = [info.strip() for info in line.split('\t')]
            marker = line_info[0]
            linkage_group_name = line_info[2]
            marker_position = line_info[4]
            sequence = line_info[14]
            snp_positions = sorted(set([int(snp.split(',')[0]) for snp in line_info[15][1:-1].split(';') if 'N/' not in snp.split(',')]))

            if snp_positions:
                alleles = [allele for allele in line_info[16].split(';')]
                
                if any([True for allele in alleles if len(allele)>1]):

                    # gener
                    alleles_by_position = []
                    for _, allele_pairs in enumerate(alleles):
                        for i_allele, allele in enumerate(allele_pairs):
                            try:
                                alleles_by_position[i_allele].append(allele)
                            except IndexError:
                                alleles_by_position.append([])
                                alleles_by_position[i_allele].append(allele)
                    
                    alleles_by_position:list = [sorted(set(allele_position)) for allele_position in alleles_by_position]
                else: alleles = [sorted(set([allele for allele in alleles]))]
            else:
                alleles = []

            if linkage_group_name not in linkage_map_dict: linkage_map_dict[linkage_group_name] = {}
            linkage_map_dict[linkage_group_name][marker] = {
                'marker_name': marker,
                'linkage_group_name': linkage_group_name,
                'marker_position': marker_position,
                'sequence': sequence,
                'SNP_positions': snp_positions,
                'alleles': alleles
            }
    
    return linkage_map_dict
# --------------------------------------------------
def main() -> None:
    """ Insert docstring here """

    args = get_args()

    linkage_map_dict: dict = _parse_linkage_map_csv(args.csv_path)

    print('\t'.join([
        '#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO'
        ]))
    for linkage_group_name, linkage_group in linkage_map_dict.items():
        chromosome = _identify_chromosome(linkage_group)
        for marker in linkage_group.values():
            _identify_markers(marker, chromosome)

    return None
# --------------------------------------------------
if __name__ == '__main__':
    main()