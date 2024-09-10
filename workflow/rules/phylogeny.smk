rule:
    input:
        ancient('results/{species}/variants/candidate_variants.duckdb'),
    params:
        maf=".85",
        qual=30,
    output:
        'results/{species}/aligned_pseudogenomes/{sequencing}.fas',
    resources:
        cpus_per_task=32,
        mem_mb=16_000,
        runtime=15,
    localrule: False
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
    
    Notes:
      - [GitHub](https://github.com/amkozlov/raxml-ng)
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