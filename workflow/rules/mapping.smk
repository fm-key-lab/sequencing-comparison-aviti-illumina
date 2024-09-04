checkpoint mapping_samplesheet:
    input:
        'results/samplesheet.csv'
    output:
        'results/{species}/samplesheet.csv'
    localrule: True
    run:
        import pandas as pd

        samplesheet = pd.read_csv(input[0])

        if wildcards.species == 'Control':
            species_mask = samplesheet['species'].isna()
        else:
            species_mask = samplesheet['species'] == species_key.get(wildcards.species)

        (
            samplesheet[species_mask]
            .filter(['sample', 'fastq_1', 'fastq_2'])
            .to_csv(output[0], index=False)
        )


use rule bactmap from widevariant as mapping with:
    input:
        input='results/{species}/samplesheet.csv',
        reference=lambda wildcards: config['reference'][wildcards.species],
    params:
        pipeline='bactmap',
        profile='singularity',
        nxf='-work-dir results/{species}/work -config config/bactmap.config',
        outdir='results/{species}',
    output:
        'results/{species}/pipeline_info/pipeline_report.txt',
    localrule: True
    envmodules:
        'apptainer/1.3.2',
        'nextflow/21.10',
        'jdk/17.0.10'


rule:
    """Collect mapping output.
    
    Collect mapping output such that the sample wildcard can be
    resolved by downstream rules.
    """
    input:
        'results/{species}/pipeline_info/pipeline_report.txt'
    output:
        touch('results/{species}/fastp/{sample}_1.trim.fastq.gz'),
        touch('results/{species}/fastp/{sample}_2.trim.fastq.gz'),
    localrule: True