#!/bin/env python3

from scitq.workflow import Workflow
from scitq.fetch import list_content, check_uri
import typer
import requests
import os
import subprocess
import io
import csv
from typing import Optional


BIOPROJECT = 'PRJNA987980'
SEED = 42

# this requires a recent version of scitq (>1.4) with protofilters activated, else providers, regions and flavors must be manually provided
DEFAULT_PROVIDER = 'auto'
DEFAULT_REGIONS = 'auto'
# you may want to adapt resources to your profiler needs here
DEFAULT_FLAVOR = 'auto:cpu>=32:ram>=125:disk>=400'

# YOU MUST ADAPT THIS PARAMETER TO A PROPER CLOUD LIKE URI LIKE azure://... OR s3://... SEE scitq documentation for details
BASE_PATH = 'azure://rnd/results'
# in case your custom profiler need a specific resource, it is generaly provided as a tar archive - that scitq will untar for
# you if you end the path with |untar (you may also use |unzip in case of a zip archive, see scitq documentation for details)
CUSTOM_RESOURCE = 'azure://path/to/my/catalogue.tgz|untar'
# you should also provide a docker containing your profiler
DOCKER_NAME = 'mydocker_name'

def count(items, ending):
    """Count how many items end with ending"""
    return len([item for item in items if item.endswith(ending)])

def custom_profiler(simulation:str, provider:str=DEFAULT_PROVIDER,
              region:Optional[str]=DEFAULT_REGIONS, scitq_project:Optional[str]=None,
              max_workflow_workers:int=10):
    """
    simulation must be refKrak, refMet4 or refBioms.

    --provider can be either azure or ovh

    --region should stay as auto unless in very specific cases

    --scitq-project is the name of the computational batch of scitq it defaults to f'custom_profiler-{simulation}'
    (if bioproject string contains some URI like pattern like s3://... only the last folder of the URI is retained)

    --max-workflow-workers is the maximal number of worker instances to deploy
   
    by default the script uses the ENA but --sra option means to use the SRA
    """
    assert simulation in ['refKrak','refMet4','refBioms']
    paired = None
    depth_string = depth
    if depth.startswith('2x'):
        paired = True
        depth = depth[2:]
    elif depth.startswith('1x'):
        paired = False
        depth = depth[2:]
    if depth=='None':
        depth=None
    else:
        depth = depth.lower().replace('k','000').replace('m','000000')
        try:
            int(depth)
        except ValueError:
            raise RuntimeError(f'Invalid depth form {depth_string}')




    ######################################################
    #                                                    #
    #    Collecting samples                              #
    #                                                    #
    ######################################################
    use_run_accessions = True

    p=subprocess.run(f'''docker run --rm -it ncbi/edirect sh -c "esearch -db sra -query '{BIOPROJECT}[bioproject]' | efetch -format runinfo"''',
                                    shell=True, check=True, capture_output=True, encoding='UTF-8')
    output=io.StringIO(p.stdout)
    samples = {}
    if paired is None:
        paired = True

    for item in csv.DictReader(output, delimiter=','):
        if item['LibraryStrategy']=='WGS':
            if paired != (item['LibraryLayout']=='PAIRED'):
                print(f'Rejecting run {item["Run"]} as it does not fit the main strategy')
                continue
            if simulation=='refMet4' and not item['LibraryName'].endswith('M'):
                continue
            if simulation=='refKrak' and not item['LibraryName'].endswith('K'):
                continue
            if simulation=='refBioms' and not item['LibraryName'].endswith('B'):
                continue
            
            if item['BioSample'] not in samples:
                samples[item['BioSample']]=[]
            samples[item['BioSample']].append(item['Run'])
    fetch_option='@sra'

    ######################################################
    #                                                    #
    #    Analytical Workflow                             #
    #                                                    #
    ######################################################

    if scitq_project is None:
        scitq_project = f'custom_profiler-{simulation}'

    azure_base = f'{BASE_PATH}/{scitq_project}' if output is None else output

    wf = Workflow(name=scitq_project, shell=True, 
                max_step_workers=5, retry=5, flavor=DEFAULT_FLAVOR, 
                provider=provider, region=region,
                max_workflow_workers=max_workflow_workers)

    for sample,runs in samples.items():

        # alignment step
        step1 = wf.step(
            batch='profiler',
            name=f'profiler:{sample}',
            command=f'cat /input/*.fastq |my_custom_profiler --input_type fastq \
                --no_map --offline --bowtie2db /resource/metaphlan/bowtie2 \
                --nproc $CPU -o /output/{sample}.metaphlan4_profile.txt',
            container=DOCKER_NAME,
            concurrency=6,
            input=[f'run+fastq{fetch_option}://{run}' for run in runs] if use_run_accessions else runs,
            output=f'{azure_base}/{sample}/metaphlan/',
            resource=CUSTOM_RESOURCE,
        )

    # final collection step
    step2 = wf.step(
        batch='compile',
        name='compile',
        shell=True,
        command=f'''cd /input && merge_metaphlan_tables.py *profile.txt > /output/merged_abundance_table.tsv''',
        container=DOCKER_NAME,
        concurrency=1,
        required_tasks=step1.gather(),
        input=step1.gather('output'),
        output=f'{azure_base}/compile/',
    )


    ######################################################
    #                                                    #
    #    Monitoring and post-treatment                   #
    #                                                    #
    ######################################################

    wf.run(refresh=10)
    os.makedirs(simulation,exist_ok=True)
    step2.download(destination=os.path.join(os.getcwd(),simulation)+'/')
    wf.clean()

    print(f'All done! Results are in {azure_base}/compile/ !')


if __name__=='__main__':
    typer.run(custom_profiler)