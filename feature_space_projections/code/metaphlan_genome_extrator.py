# all files are in http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/mpa_vOct22_CHOCOPhlAnSGB_202212.tar
import pysam

print('Loading marker genes')
with open('mpa_vOct22_CHOCOPhlAnSGB_202212_marker_info.txt','r',encoding='utf-8') as f:
    marker_genes_lines = [l for l in f.read().split('\n') if l]

marker_genes = {}
for line in marker_genes_lines:
    # the second part of the line is not proper json (single quote instead of double), but is a kind of python dump
    # UniRef90_UPI000E65AEFC|1__27|SGB32561   {'clade': 't__SGB32561', 'ext': [], 'len': 4050, 'taxon': 'k__Bacteria|p__Proteobacteria|c__Betaproteobacteria|o__Burkholderiales|f__Comamonadaceae|g__Calidifontimicrobium|s__Calidifontimicrobium_sediminis|t__SGB32561'}
    gene,description = line.split('\t')
    taxon=eval(description)['taxon']
    if taxon not in marker_genes:
        marker_genes[taxon]=[gene]
    else:
        marker_genes[taxon].append(gene)

print('Loading metaphlan genomes')
genomes = pysam.FastaFile('mpa_vOct22_CHOCOPhlAnSGB_202212_SGB.fna')

print('Loading species list to create')
with open('species.txt','r',encoding='utf-8') as f:
    species = [l for l in f.read().split('\n') if l and l!='clade_name']

print('Creating genomes')
msg = []
for specie in species:
    if specie in marker_genes:
        genes = marker_genes[specie]
        clade = specie.split('|')[-1]
        with open(f'{clade}.fa','w') as f:
            for gene in genes:
                try:
                    gene_seq = genomes.fetch(gene)
                except KeyError:
                    msg.append(f'Could not find gene {gene} for specie {specie}')
                    continue
                f.write(f'>{gene}\n')
                f.write(genomes.fetch(gene))
                f.write('\n')
    else:
        msg.append(f'Could not find specie {specie}')

print('\n'.join(msg))