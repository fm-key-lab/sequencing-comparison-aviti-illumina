mlst_key = {
    'Escherichia coli': 'Escherichia_coli#1',
}

use rule srst2 from widevariant as sequence_typing with:
    input:
        pe=multiext(
            'results/{species}/fastp/{sample}',
            '_1.trim.fastq.gz',
            '_2.trim.fastq.gz',
        ),
        db=lambda wildcards: config['public_data']['mlst'][wildcards.species]['db'],
        prof=lambda wildcards: config['public_data']['mlst'][wildcards.species]['definitions'],
    params:
        extra='--forward "_1.trim" --reverse "_2.trim" --mlst_delimiter "_"',
        species_alias=lambda wildcards: config['public_data']['mlst'][wildcards.species]['alias'],
        prefix='results/{species}/mlst/{sample}'
    output:
        'results/{species}/mlst/{sample}__results.txt'
    envmodules:
        'srst2/0.2.0'


def sequence_typing_output(wildcards):
    import pandas as pd

    sample_ids = pd.read_csv(
        checkpoints.mapping_samplesheet.get(
            species=wildcards.species,
        ).output[0]
    )['sample'].astype(str)

    return expand(
        'results/{{species}}/mlst/{sample}__results.txt',
        sample=sample_ids
    )


rule:
    """Collect sequence typing output.
    
    Collect sequence typing output to resolve the sample wildcard.
    """
    input:
        sequence_typing_output
    output:
        touch('results/{species}/mlst/.done'),
    localrule: True
