#!/bin/env python3

import pysam
import sys
import typer
from typing import List
import re

# typical use for motus is : motus_genome_extractor.py ../motus_species_list.tsv ../db_mOTU_DB_CEN.fasta

MAX_FILE = 100

def get_genomes(path):
    with open(path,'r') as f:
        return list(map(lambda x: x.strip(),f.readlines()))
    
def main(genomes_file: str, catalog_files: List[str], reference_regexp:str=r'(.*?)\..*'):
    """
    --reference_regexp: a regular expression with a group (e.g. bound with parenthesis) which define where the genome name is within each sequence ref
    """
    genomes = get_genomes(genomes_file)
    reference_match=re.compile(reference_regexp)
    #print(repr(genomes))
    for catalog_file in catalog_files:
        print(f'\Processing catalog {catalog_file}',end='')
        reference_catalog=pysam.FastaFile(catalog_file)
        genomes_done = []
        file_cache = []
        remaining_genomes = list(genomes)
        for reference in reference_catalog.references:
            m=reference_match.match(reference)
            if m:
                ref_genome = m.groups()[0]
                if ref_genome in genomes:
                    if ref_genome in remaining_genomes:
                        #print(f'Found genome {ref_genome}')
                        print('.',end='')
                        remaining_genomes.remove(ref_genome)
                    for cache_genome, genome_file in file_cache:
                        if cache_genome == ref_genome:
                            break
                    else:
                        if len(file_cache)>=MAX_FILE:
                            _, genome_file = file_cache.pop(0)
                            genome_file.close()
                        genome_file = open(f'{ref_genome}.fa','a')
                        if ref_genome not in genomes_done:
                            genome_file.write(f'>{ref_genome}\n')
                        file_cache.append((ref_genome, genome_file))
                    if ref_genome not in genomes:
                        genomes_done.append(ref_genome)

                    genome_file.write(reference_catalog.fetch(reference)+'\n')     

    for _,genome_file in file_cache:
            genome_file.close()
    if remaining_genomes:
        print(f'Those genomes were not found: {", ".join(remaining_genomes)}')

if __name__ == '__main__':
    typer.run(main)
