rule:
    input:
        'results/{species}/variants/{sample}.vcf.gz',
    output:
        'results/{species}/variants/{sample}_af.tsv'
    resources:
        cpus_per_task=2,
        runtime=5
    localrule: False
    envmodules:
        'bcftools/1.20'
    shell:
        '''
        bcftools +fill-tags {input} -- -t AF | \
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE]\t%QUAL\t%DP\t%INFO/DP4\n' > {output}
        '''


rule:
    input:
        'results/{species}/samtools/{sample}.sorted.bam',
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

    sample_ids = pd.read_csv(
        checkpoints.mapping_samplesheet.get(
            species=wildcards.species,
        ).output[0]
    )['sample'].astype(str)

    return expand(
        [
            'results/{{species}}/variants/{sample}_af.tsv',
            'results/{{species}}/samtools/{sample}_genomecov.tsv',
        ],
        sample=sample_ids,
    )


rule:
    input:
        candidate_variant_tables
    params:
        af_glob="'results/{species}/variants/*_af.tsv'",
        gc_glob="'results/{species}/samtools/*_genomecov.tsv'",
    output:
        'results/{species}/variants/candidate_variants.duckdb',
    resources:
        cpus_per_task=32,
        mem_mb=48_000,
        runtime=15
    localrule: False
    envmodules:
        'duckdb/nightly'
    shell:
        '''
        export MEMORY_LIMIT="$(({resources.mem_mb} / 1000))GB" \
               BCFTOOLS_QUERY={params.af_glob} \
               BEDTOOLS_GENOMECOV={params.gc_glob}
        duckdb {output} -c ".read workflow/scripts/create_variants_db.sql"
        '''