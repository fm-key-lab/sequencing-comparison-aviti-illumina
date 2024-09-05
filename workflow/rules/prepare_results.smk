OUTPUT = [
    'samplesheet',
    'sequence_typing',
    'total_reads',
    'variant_quality_scores',
]


rule record_samples:
    input:
        'results/samplesheet.csv'
    output:
        temp('results/samplesheet.duckdb'),
    resources:
        cpus_per_task=1,
        mem_mb=1_024
    localrule: True
    envmodules:
        'duckdb/nightly'
    shell:
        '''
        duckdb {output} -s \
          "set memory_limit = '$(({resources.mem_mb} / 1000))GB';
          set threads = {resources.cpus_per_task};
          create table samplesheet as 
          select * 
          from read_csv('{input}');"
        '''


rule record_total_reads:
    input:
        expand(
            'results/{species}/multiqc/multiqc_data/multiqc_fastp.yaml',
            species=config['wildcards']['species_tmp'].split('|')
        )
    output:
        temp('results/total_reads.duckdb'),
    resources:
        cpus_per_task=8,
        mem_mb=16_000
    envmodules:
        'yq/4.44.3',
        'duckdb/nightly'
    shell:
        '''
        export MEMORY_LIMIT="$(({resources.mem_mb} / 1000))GB"
        yq -o=csv '.[] | [key, .summary.before_filtering.total_reads, .summary.after_filtering.total_reads]' {input} |\
          duckdb {output} -c ".read workflow/scripts/parse_fastp_multiqc_output.sql"
        '''


rule record_variant_quality_scores:
    input:
        expand(
            'results/{species}/multiqc/multiqc_data/mqc_bcftools_stats_vqc_Count_SNP.yaml',
            species=config['wildcards']['species_tmp'].split('|')
        )
    output:
        temp('results/variant_quality_scores.duckdb'),
    resources:
        cpus_per_task=8,
        mem_mb=16_000
    envmodules:
        'yq/4.44.3',
        'duckdb/nightly'
    shell:
        '''
        export MEMORY_LIMIT="$(({resources.mem_mb} / 1000))GB"
        yq -o=csv 'to_entries[] | .key as $sample | .value | to_entries[] | [$sample, .key, .value]' {input} |\
          duckdb {output} -c ".read workflow/scripts/parse_bcftools_variant_quality.sql"
        '''


rule record_sequence_typing_results:
    input:
        expand(
            'results/{species}/mlst/.done',
            species=config['wildcards']['species_tmp'].split('|')
        )
    params:
        mlst_glob='results/*/mlst/*[0-9][!a-zA-Z]__results.txt',
    output:
        temp('results/sequence_typing.duckdb'),
    resources:
        cpus_per_task=4,
        mem_mb=4_096
    envmodules:
        'duckdb/nightly'
    shell:
        '''
        export MEMORY_LIMIT="$(({resources.mem_mb} / 1000))GB"
        export MLST_RESULTS="{params.mlst_glob}"
        duckdb {output} -c ".read workflow/scripts/parse_srst2_output.sql"
        '''


rule:
    """Collect tables into one database.

    Avoids concurrency issues with multiple rules writing to the same database.
    """
    input:
        expand(
            'results/{output}.duckdb',
            output=OUTPUT
        )
    output:
        'results/results.duckdb',
    resources:
        cpus_per_task=1,
        mem_mb=1_024
    localrule: True
    envmodules:
        'duckdb/nightly'
    shell:
        '''
        for db in {input}; do
          duckdb -s "\
          set memory_limit = '$(({resources.mem_mb} / 1000))GB';
          set threads = {resources.cpus_per_task};
          attach '{output}' as results_db;
          attach '${{db}}' as output_db;
          copy from database output_db to results_db;"
        done
        '''