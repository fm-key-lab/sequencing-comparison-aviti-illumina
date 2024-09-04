rule vcf2pseudogenome:
    input:
        ancient('results/{species}/variants/{sample}.filtered.vcf.gz'),
    output:
        'results/{species}/pseudogenomes/{sample}.fas',
    params:
        reference=lambda wildcards: config['reference'][wildcards.species],
    envmodules:
        'sandbox'
    localrule: True
    shell:
        '''
        python workflow/scripts/vcf2pseudogenome.py  -r {params.reference} -b {input} -o {output}
        '''


rule:
    input:
        'results/{species}/pseudogenomes/{sample}.fas',
    output:
        'results/{species}/{group}/pseudogenomes/{sample}.fas'
    localrule: True
    shell:
        'cp {input} {output}'


def align_pseudogenomes_input(wildcards):
    """Align samples from the same cohorts, species."""
    import pandas as pd

    samplesheet = pd.read_csv(
        checkpoints.samplesheet.get(
            species=wildcards.species,
        ).output[0]
    )

    sample_ids = samplesheet[
        samplesheet['group'] == wildcards.group
    ]['sample'].astype(str)

    return expand(
        'results/{{species}}/{{group}}/pseudogenomes/{sample}.fas',
        sample=sample_ids
    )


rule align_pseudogenomes:
    input:
        align_pseudogenomes_input
    output:
        aligned='results/{species}/{group}/aligned_pseudogenomes/aligned_pseudogenome.fas',
        final_ref='results/{species}/{group}/aligned_pseudogenomes/final_reference.fas',
    params:
        reference=lambda wildcards: config['reference'][wildcards.species],
    envmodules:
        'sandbox'
    localrule: True
    shell:
        '''
        touch {output.aligned}
        for pseudogenome in {input}
        do
            cat $pseudogenome >> {output.aligned}
        done
        python workflow/scripts/reference2single_sequence.py -r {params.reference} -o {output.final_ref}
        cat {output.final_ref} >> {output.aligned}
        '''


use rule gubbins from widevariant as remove_recombination with:
    input:
        'results/{species}/{group}/aligned_pseudogenomes/aligned_pseudogenome.fas',
    params:
        f=config['gubbins']['filter_percentage'],
        tree_args=config['gubbins']['tree_args'],
        t=config['gubbins']['tree_builder']
    output:
        'results/{species}/{group}/gubbins/prefix.final_tree.tre',
    envmodules:
        'intel/21.2.0',
        'impi/2021.2',
        'gubbins/3.3.5'


rule:
    input:
        expand(
            'results/{{species}}/{group}/gubbins/prefix.final_tree.tre',
            group=['aviti', 'illumina'],
        )
    output:
        touch('results/{species}/phylogeny.done')