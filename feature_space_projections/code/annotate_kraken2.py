import glob
import requests
import csv
import logging

# this small scripts should be run in the folder that is downloaded
# by the use of the scitq_kraken2.py script explained in the main README.md

print('specie\tgtdbTaxonomy')
for specie in glob.glob('*'):
    try:
        with open(f'{specie}/{specie}.report', 'r') as f:
            taxonomy = []
            for line in csv.reader(f, dialect='excel-tab'):
                rank = line[3].lower()
                value = line[5].strip()
                if rank in ['u','r']:
                    continue
                taxonomy.append(f'{rank}__{value}')
                if rank=='s': 
                    break
            else:
                taxonomy.append('unclassified')
            print(f'''{specie}\t{';'.join(taxonomy)}''')
    except Exception as e:
        logging.error(f'Could not parse {specie}')
        logging.exception(e)

