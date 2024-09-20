rule:
    output:
        'resources/raw_samplesheet.xlsx'
    params:
        url=config['sample_info_url']
    resources:
        slurm_partition='datatransfer'
    localrule: True
    shell:
        'wget -nc -O {output} {params.url}'


checkpoint samplesheet:
    input:
        ancient('resources/raw_samplesheet.xlsx')
    output:
        'results/samplesheet.csv'
    localrule: True
    run:
        from utils import create_samplesheet
        
        samplesheet = create_samplesheet(input[0], config['data'])
        
        samplesheet['species'] = samplesheet['species'].str.replace(' ', '_')

        samplesheet = samplesheet[
            samplesheet['donor'].isin(config['wildcards']['donors'].split('|')) &
            samplesheet['species'].isin(config['wildcards']['species'].split('|'))
        ]

        samplesheet.to_csv(output[0], index=False)