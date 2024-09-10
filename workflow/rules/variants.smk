rule:
    input:
        ancient('results/{species}/variants/{sample}.vcf.gz')
    output:
        temp('results/{species}/variants/{sample}_af.tsv')
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


def variant_tables(wildcards):
    import pandas as pd

    sample_ids = pd.read_csv(
        checkpoints.mapping_samplesheet.get(
            species=wildcards.species,
        ).output[0]
    )['sample'].astype(str)

    return expand(
        'results/{{species}}/variants/{sample}_af.tsv',
        sample=sample_ids,
    )


rule:
    input:
        ancient(variant_stats_output)
    params:
        glob="'results/{species}/variants/*_af.tsv'",
    output:
        'results/{species}/variants/candidate_variants.duckdb',
    resources:
        cpus_per_task=48,
        mem_mb=64_000,
        runtime=30
    localrule: False
    envmodules:
        'duckdb/nightly'
    shell:
        '''
        export MEMORY_LIMIT="$(({resources.mem_mb} / 1000))GB" BCFTOOLS_QUERY={params.glob}
        duckdb {output} -c ".read workflow/scripts/create_variants_db.sql"
        '''


rule:
    input:
        'results/{species}/variants/candidate_variants.duckdb',
    params:
        maf=".85",
        qual=30,
    output:
        'results/{species}/aligned_pseudogenomes/aligned_pseudogenome.fas',
    localrule: True
    envmodules:
        'duckdb/nightly'
    shell:
        '''
        export MEMORY_LIMIT="$(({resources.mem_mb} / 1000))GB" MAF_THRESHOLD=".85" QUAL_THRESHOLD=30 \
        duckdb {output} -c ".read workflow/scripts/finalize_variants.sql" > {output}
        '''


rule:
    input:
        'results/{species}/aligned_pseudogenomes/aligned_pseudogenome.fas',
    output:
        'results/{species}/aligned_pseudogenomes/aligned_pseudogenome.phy',
    run:
        from Bio import AlignIO, SeqIO
        alignment = AlignIO.read(input[0], 'fasta')
        with open(output[0], 'w') as f:
            SeqIO.write(alignment, f, 'phylip')