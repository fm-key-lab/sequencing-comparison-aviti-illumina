rule:
    input:
        ancient('results/{species}/fastp/{sample}.fastp.json')
    output:
        'results/{species}/candidate_variant_table/base_freq/{sample}.csv'
    resources:
        cpus_per_task=2,
        runtime=5
    localrule: False
    envmodules:
        'yq/4.44.3'
    shell:
        '''
        touch {output}
        for orient in 1 2; do
          for base in A C T G N; do
            query=".read${{orient}}_before_filtering.content_curves.${{base}}"
            yq -o=csv $query {input} | awk -F',' -v orient=$orient -v base=$base '{{sum=0; for(i=1;i<=NF;i++) sum+=$i; avg=sum/NF; print avg","orient","base}}' >> {output}
          done
        done
        '''


rule:
    input:
        ancient('results/{species}/samtools/{sample}.sorted.bam')
    output:
        'results/{species}/candidate_variant_table/coverage/{sample}_strand.tsv',
        'results/{species}/candidate_variant_table/coverage/{sample}_total.tsv',
    resources:
        cpus_per_task=2,
        runtime=5
    localrule: False
    envmodules:
        'bedtools/2.31.1'
    shell:
        '''
        touch {output[0]}
        for strand in "+" "-" ""; do
          genomeCoverageBed -bg -strand $strand -ibam {input} | awk -v strand=$strand '{{print $0, strand}}' OFS='\t' >> {output[0]}
        done
        genomeCoverageBed -bg -ibam {input} > {output[1]}
        '''


rule:
    input:
        ancient('results/{species}/variants/{sample}.vcf.gz')
    output:
        'results/{species}/candidate_variant_table/allele_freq/{sample}.tsv'
    resources:
        cpus_per_task=2,
        runtime=5
    localrule: False
    envmodules:
        'bcftools/1.20'
    shell:
        '''
        bcftools +fill-tags {input} -- -t AF | \
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE]\t%QUAL\t%AF\t%DP\t%INFO/DP4\n' > {output}
        '''


def variant_stats_output(wildcards):
    import pandas as pd

    checkpoint_output = checkpoints.mapping_samplesheet.get(
        species=wildcards.species,
    ).output[0]

    sample_ids = pd.read_csv(checkpoint_output)['sample'].astype(str)

    return expand(
        [
            'results/{{species}}/candidate_variant_table/base_freq/{sample}.csv',
            'results/{{species}}/candidate_variant_table/coverage/{sample}_strand.tsv',
            'results/{{species}}/candidate_variant_table/coverage/{sample}_total.tsv',
            'results/{{species}}/candidate_variant_table/allele_freq/{sample}.tsv',
        ],
        sample=sample_ids,
    )


rule:
    input:
        ancient(variant_stats_output)
    params:
        # NOTE: Pretty sure nested quotes required here for DuckDB
        base_freq_glob="'results/{species}/candidate_variant_table/base_freq/*.csv'",
        coverage_strand_glob="'results/{species}/candidate_variant_table/coverage/*_strand.tsv'",
        coverage_total_glob="'results/{species}/candidate_variant_table/coverage/*_total.tsv'",
        allele_freq_glob="'results/{species}/candidate_variant_table/allele_freq/*.tsv'",
    output:
        'results/candidate_variant_table/{species}.duckdb',
    resources:
        cpus_per_task=32,
        mem_mb=64_000,
        runtime=120
    localrule: False
    envmodules:
        'duckdb/nightly'
    shell:
        '''
        export MEMORY_LIMIT="$(({resources.mem_mb} / 1000))GB"
        export BASE_FREQ={params.base_freq_glob}
        export COVERAGE={params.coverage_strand_glob}
        export COVERAGE_TOT={params.coverage_total_glob}
        export ALLELE_FREQ={params.allele_freq_glob}
        duckdb {output} -c ".read workflow/scripts/create_candidate_variant_tbls.sql"
        '''