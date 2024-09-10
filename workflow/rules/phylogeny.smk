rule pseudogenome_alignment:
    input:
        samplesheet=ancient('results/samplesheet.csv'),
        db=ancient('results/{species}/variants/candidate_variants.duckdb'),
    params:
        covg=20,
        maf=".85",
        qual=30,
    output:
        'results/{species}/aligned_pseudogenomes/{sequencing}.fas',
    resources:
        cpus_per_task=1,
        mem_mb=64_000,
        runtime=30,
    envmodules:
        'duckdb/nightly'
    shell:
        '''
        export MEMORY_LIMIT="$(({resources.mem_mb} / 1100))GB" \
               COVERAGE_THRESHOLD={params.covg} \
               MAF_THRESHOLD={params.maf} \
               QUAL_THRESHOLD={params.qual} \
               SEQUENCING={wildcards.sequencing} \
               SAMPLESHEET={input.samplesheet}
        duckdb -readonly {input.db} -c ".read workflow/scripts/finalize_variants.sql" > {output}
        '''


rule extract_variant_sites:
    input:
        ancient('results/{species}/aligned_pseudogenomes/{sequencing}.fas'),
    output:
        'results/{species}/snpsites/{sequencing}.filtered_alignment.fas',
    localrule: True
    envmodules:
        'snp-sites/2.5.1'
    shell:
        '''
        snp-sites -o {output} {input}
        '''


rule veryfasttree:
    """Run VeryFastTree.

    Build a phylogeny from the multiple sequence alignment using the 
    VeryFastTree implementation of the FastTree-2 algorithm.

    Args:

    Returns:
    
    Notes:
      - [GitHub](https://github.com/citiususc/veryfasttree)
    """
    input:
        ancient('results/{species}/snpsites/{sequencing}.filtered_alignment.fas'),
    output:
        'results/{species}/veryfasttree/{sequencing}.veryfasttree.phylogeny.nhx',
    resources:
        cpus_per_task=32,
        mem_mb=64_000,
        runtime=30,
    envmodules:
        'veryfasttree/4.0.3.1'
    shell:
        '''
        export OMP_PLACES=threads
        veryfasttree {input} -nt -threads {resources.cpus_per_task} > {output}
        '''


rule raxml_ng:
    """Run RAxML-NG.

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
        ancient('results/{species}/snpsites/{sequencing}.filtered_alignment.fas'),
    params:
        extra='--all --model GTR+G --bs-trees 200',
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
        runtime=120,
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