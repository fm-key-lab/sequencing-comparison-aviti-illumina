rule bcftools_fill_af_tag_query:
    input:
        ancient('results/{species}/variants/{sample}.vcf.gz'),
    output:
        'results/{species}/variants/{sample}_afreq.tsv'
    resources:
        cpus_per_task=2,
        runtime=5
    envmodules:
        'bcftools/1.20'
    shell:
        '''
        bcftools +fill-tags {input} -- -t AF | \
          bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE]\t%QUAL\t%INFO\t%DP\t%INFO/DP4\n' > {output}
        '''


rule genome_coverage_bed:
    input:
        ancient('results/{species}/samtools/{sample}.sorted.bam'),
    output:
        'results/{species}/samtools/{sample}_genomecov.tsv'
    resources:
        cpus_per_task=2,
        runtime=5
    envmodules:
        'bedtools/2.31.1'
    shell:
        '''
        touch {output}
        
        for strand in "+" "-"; do
          genomeCoverageBed -bga -strand $strand -ibam {input} | \
          awk -v strand=$strand '{{print $0, strand}}' OFS='\t' >> {output}
        done
        '''


def candidate_variant_tables(wildcards):
    import pandas as pd

    samplesheet = pd.read_csv(
        checkpoints.samplesheet.get(
            **wildcards,
        ).output[0]
    )

    samplesheet = samplesheet[
        ~samplesheet['sample'].astype(int).isin(EXCLUDE)
    ]

    return list(
        samplesheet
        [samplesheet['species'].isin(config['wildcards']['species'].split('|'))]
        .filter(['sample', 'species'])
        .drop_duplicates()
        .transpose()
        .apply(
            lambda df: [
                output_fn.format(**df.to_dict()) for output_fn in [
                    'results/{species}/variants/{sample}.filtered.vcf.gz',
                    'results/{species}/variants/{sample}_afreq.tsv',
                    'results/{species}/samtools/{sample}_genomecov.tsv'
                ]
            ]
        )
        .values
        .flatten()
    )


rule create_variants_db:
    input:
        ancient(candidate_variant_tables)
    params:
        vcfs="'results/*/variants/*.filtered.vcf.gz'"
    output:
        'results/candidate_variants.duckdb',
    resources:
        cpus_per_task=32,
        mem_mb=96_000,
        runtime=120
    envmodules:
        'duckdb/nightly'
    shell:
        '''
        export MEMORY_LIMIT="$(({resources.mem_mb} / 1200))GB" \
               FILTERED_VCFS={params.vcfs}
        
        duckdb {output} -c ".read workflow/scripts/create_variants_db.sql"
        '''