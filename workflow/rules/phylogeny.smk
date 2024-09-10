rule:
    input:
        ancient('results/{species}/variants/candidate_variants.duckdb'),
    params:
        maf=".85",
        qual=30,
    output:
        'results/{species}/aligned_pseudogenomes/{sequencing}.fas',
    localrule: True
    envmodules:
        'duckdb/nightly'
    shell:
        '''
        export MEMORY_LIMIT="$(({resources.mem_mb} / 1000))GB" \
               MAF_THRESHOLD=".85" \
               QUAL_THRESHOLD=30 \
               SEQUENCING={wildcards.sequencing}
        duckdb {input} -c ".read workflow/scripts/finalize_variants.sql" > {output}
        '''


rule veryfasttree:
    input:
        'results/{species}/aligned_pseudogenomes/{sequencing}.fas',
    output:
        'results/{species}/veryfasttree/{sequencing}.veryfasttree.phylogeny.nhx',
    resources:
        cpus_per_task=32,
        mem_mb=8_000,
        runtime=15,
    localrule: False
    envmodules:
        'veryfasttree/4.0.3.1'
    shell:
        '''
        export OMP_PLACES=threads
        veryfasttree {input} -nt -threads {resources.cpus_per_task} > {output}
        '''


rule raxml_ng:
    """Run RAxML-NG .

    Build a maximum-likelihood phylogeny from the multiple sequence alignment 
    using RAxML Next Generation.

    Args:

    Returns:
        {{ prefix }}.raxml.reduced.phy: Reduced alignment (with duplicates 
          and gap-only sites/taxa removed)
        {{ prefix }}.raxml.rba: Binary MSA file
        {{ prefix }}.raxml.bestTreeCollapsed: Best ML tree with collapsed near-zero branches
        {{ prefix }}.raxml.bestTree: Best ML tree
        {{ prefix }}.raxml.mlTrees: All ML trees
        {{ prefix }}.raxml.support: Best ML tree with Felsenstein bootstrap (FBP) support
        {{ prefix }}.raxml.bestModel: Optimized model
        {{ prefix }}.raxml.bootstraps: Bootstrap trees
        {{ prefix }}.raxml.log: Execution log
    """
    input:
        'results/{species}/aligned_pseudogenomes/{sequencing}.fas',
    params:
        extra='--all --model GTR+G --bs-trees 1000',
        prefix='results/{species}/raxml_ng/{sequencing}',
    output:
        multiext(
            'results/{species}/raxml_ng/{sequencing}.raxml',
            '.reduced.phy',
            '.rba',
            '.bestTreeCollapsed',
            '.bestTree',
            '.bestModel',
            '.log'
        )
    resources:
        cpus_per_task=48,
        mem_mb=64_000,
        runtime=60,
    localrule: False
    envmodules:
        'raxml-ng/1.2.2_MPI'
    shell:
        '''
        export OMP_PLACES=threads
        raxml-ng {params.extra} --msa {input} --threads {resources.cpus_per_task} --prefix {params.prefix}
        '''


rule:
    input:
        expand(
            [
                'results/{{species}}/veryfasttree/{sequencing}.veryfasttree.phylogeny.nhx',
                'results/{{species}}/raxml_ng/{sequencing}.raxml.bestTree',
            ],
            sequencing=config['wildcards']['sequencing'].split('|'),
        )
    output:
        touch('results/{species}/phylogenies.done')
    localrule: True


# def pseudogenomes_input(wildcards):
#     """Align samples from the same cohorts, species."""
#     import pandas as pd

#     samplesheet = pd.read_csv(
#         checkpoints.samplesheet.get(
#             species=wildcards.species,
#         ).output[0]
#     )

#     sample_ids = samplesheet[
#         samplesheet['group'] == wildcards.group
#     ]['sample'].astype(str)

#     return expand(
#         'results/{{species}}/{{group}}/pseudogenomes/{sample}.fas',
#         sample=sample_ids
#     )


# rule align_pseudogenomes:
#     input:
#         pseudogenomes_input
#     output:
#         aligned='results/{species}/{group}/aligned_pseudogenomes/aligned_pseudogenome.fas',
#         reference_single_sequence='results/{species}/{group}/aligned_pseudogenomes/final_reference.fas',
#     params:
#         reference=lambda wildcards: config['public_data']['reference'][wildcards.species],
#     envmodules:
#         'sandbox'
#     shell:
#         '''
#         touch {output.aligned}
#         for pseudogenome in {input}
#         do
#             cat $pseudogenome >> {output.aligned}
#         done
#         python workflow/scripts/reference2single_sequence.py \
#           -r {params.reference} \
#           -o {output.reference_single_sequence}
#         cat {output.reference_single_sequence} >> {output.aligned}
#         '''


# use rule gubbins from widevariant as build_tree with:
#     input:
#         'results/{species}/{group}/aligned_pseudogenomes/aligned_pseudogenome.fas',
#     params:
#         f=config['gubbins']['filter_percentage'],
#         tree_args=config['gubbins']['tree_args'],
#         t=config['gubbins']['tree_builder']
#     output:
#         'results/{species}/{group}/gubbins/prefix.final_tree.tre',
#     envmodules:
#         'intel/21.2.0',
#         'impi/2021.2',
#         'gubbins/3.3.5'


# rule:
#     input:
#         expand(
#             'results/{{species}}/{group}/gubbins/prefix.final_tree.tre',
#             group=config['wildcards']['sequencing'].split('|'),
#         )
#     output:
#         touch('results/{species}/phylogeny.done')
#     localrule: True