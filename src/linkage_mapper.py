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
def _identify_chromosome(_key: str, _value: dict):
    """
    """

    temp_fasta = tempfile.NamedTemporaryFile()
    
    print(_value['sequence'])

    with open(temp_fasta.name, 'w') as temp_fasta_file:
        temp_fasta_file.write("> temporary query\n")
        temp_fasta_file.write(f"{_value['sequence']}")

    subprocess.run([
        "blastn",
        "-db", "data/Rubus_idaeus_JoanJ.fna",
        "-outfmt", "15",
        "-query", f"{temp_fasta.name}"])

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
            snp_positions = sorted(set([snp.split(',')[0] for snp in line_info[15][1:-1].split(';')]))

            if linkage_group_name in linkage_map_dict:
                linkage_map_dict[linkage_group_name][marker] = {
                    'marker_name': marker,
                    'linkage_group_name': linkage_group_name,
                    'marker_position': marker_position,
                    'sequence': sequence,
                    'SNP_positions': snp_positions
                }
            else:
                linkage_map_dict[linkage_group_name] = {}
                linkage_map_dict[linkage_group_name][marker] = {
                    'marker_name': marker,
                    'linkage_group_name': linkage_group_name,
                    'marker_position': marker_position,
                    'sequence': sequence,
                    'SNP_positions': snp_positions
                }
    
    return linkage_map_dict
# --------------------------------------------------
def main() -> None:
    """ Insert docstring here """

    args = get_args()

    linkage_map_dict: dict = _parse_linkage_map_csv(args.csv_path)

    for key, value in linkage_map_dict.items():
        for marker in value.values():
            _identify_chromosome(key, marker)

    return None
# --------------------------------------------------
if __name__ == '__main__':
    main()