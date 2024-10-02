rule:
    input:
        ancient('results/{species}/variants/{sample}.vcf.gz'),
    output:
        multiext(
            'results/lake/species={species}/sample={sample}/',
            'indel.parquet',
            'snp.parquet',
            'nonindels.parquet',
        )
    resources:
        cpus_per_task=4,
        runtime=10
    envmodules:
        'bcftools/1.20',
        'vcf2parquet/0.4.1'
    shell:
        '''
        bcftools view -i 'TYPE="INDEL"' {input} |\
          vcf2parquet -i /dev/stdin convert -o {output[0]}

        bcftools view -i TYPE="SNP"' {input} |\
          vcf2parquet -i /dev/stdin convert -o {output[1]}

        bcftools view -e 'TYPE="INDEL"' {input} |\
          vcf2parquet -i /dev/stdin convert -o {output[2]}
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
        .apply(lambda df: 'results/lake/species={species}/sample={sample}/nonindels.parquet'.format(**df.to_dict()))
        .values
        .flatten()
    )


rule create_variants_db:
    input:
        ancient(candidate_variant_tables)
    params:
        vcfs="'results/lake/*/*/nonindels.parquet'"
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
        # export MEMORY_LIMIT="$(({resources.mem_mb} / 1200))GB" \
        #        VCFS={params.vcfs}

        duckdb {output} -c \
          "set memory_limit = '$(({resources.mem_mb} / 1200))GB';
          set threads = {resources.cpus_per_task};
          create table tmp_variants as 
          select species, "sample", chromosome, "position", reference, alternate, quality
          from read_parquet({params.vcfs}, hive_partitioning = true, hive_types = {'species': varchar, 'sample': usmallint});"

        # duckdb {output} -c ".read workflow/scripts/create_variants_db.sql"
        '''