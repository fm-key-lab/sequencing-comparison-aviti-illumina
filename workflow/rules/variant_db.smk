rule:
    input:
        'results/{species}/fastp/{sample}.fastp.json'
    output:
        'results/{species}/candidate_variant_table/base_freq/{sample}.csv'
    envmodules:
        'yq/4.44.3'
    shell:
        '''
        touch {output}
        for orient in 1 2; do
          for base in A C T G N; do
            query=".read${{orient}}_before_filtering.content_curves.${{base}}"
            yq -o=csv $query {input} | \
              awk -F',' -v orient="$orient" -v base="$base" \
                '{sum=0; for(i=1;i<=NF;i++) sum+=$i; avg=sum/NF; print avg","orient","base}' >> {output}
          done
        done
        '''


rule:
    input:
        'results/{species}/samtools/{sample}.sorted.bam'
    output:
        'results/{species}/candidate_variant_table/coverage/{sample}.tsv'
    envmodules:
        'bedtools/2.31.1'
    shell:
        '''
        touch {output}
        for strand in "+" "-"; do
          genomeCoverageBed -bg -strand $strand -ibam {input} | \
            awk -v strand="$strand" '{print $0, strand}' OFS='\t' >> {output}
        done
        '''


rule:
    input:
        'results/{species}/variants/{sample}.vcf.gz'
    output:
        'results/{species}/candidate_variant_table/allele_freq/{sample}.tsv'
    envmodules:
        'bcftools/1.20'
    shell:
        '''
        bcftools +fill-tags {input} -- -t AF | \
          bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE]\t%QUAL\t%INFO/DP4\n' > {output}
        '''


rule:
    input:
        expand(
            [
                'results/{species}/candidate_variant_table/base_freq/{sample}.csv',
                'results/{species}/candidate_variant_table/coverage/{sample}.tsv',
                'results/{species}/candidate_variant_table/allele_freq/{sample}.tsv',
            ],
            species=['Ecoli'],
            sample=list(range(600, 620)),
        )
    output:
        touch('results/cmt_dev.done')