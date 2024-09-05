def pseudogenomes_input(wildcards):
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
        pseudogenomes_input
    output:
        aligned='results/{species}/{group}/aligned_pseudogenomes/aligned_pseudogenome.fas',
        reference_single_sequence='results/{species}/{group}/aligned_pseudogenomes/final_reference.fas',
    params:
        reference=lambda wildcards: config['public_data']['reference'][wildcards.species],
    envmodules:
        'sandbox'
    shell:
        '''
        touch {output.aligned}
        for pseudogenome in {input}
        do
            cat $pseudogenome >> {output.aligned}
        done
        python workflow/scripts/reference2single_sequence.py \
          -r {params.reference} \
          -o {output.reference_single_sequence}
        cat {output.reference_single_sequence} >> {output.aligned}
        '''


use rule gubbins from widevariant as build_tree with:
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
            group=config['wildcards']['sequencing'].split('|'),
        )
    output:
        touch('results/{species}/phylogeny.done')
    localrule: True