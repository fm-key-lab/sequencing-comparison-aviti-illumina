rule pseudogenome_alignment:
    input:
        samplesheet='results/samplesheet.csv',
        db='results/{species}/variants/candidate_variants.duckdb',
    params:
        alt=config['variants_thresh']['reads_alt'],
        covg=config['variants_thresh']['coverage'],
        idist=config['variants_thresh']['interpos_dist'],
        maf=config['variants_thresh']['maf'],
        qual=config['variants_thresh']['qual'],
    output:
        'results/{species}/aligned_pseudogenomes/{sequencing}/{donor}.fas',
    resources:
        cpus_per_task=1,
        mem_mb=64_000,
        runtime=30,
    envmodules:
        'duckdb/nightly'
    shell:
        '''
        export MEMORY_LIMIT="$(({resources.mem_mb} / 1100))GB" \
               ALT_STRAND_DP_THRESHOLD={params.alt} \
               COVERAGE_THRESHOLD={params.covg} \
               INTERBASE_DISTANCE_THRESHOLD={params.idist} \
               MAF_THRESHOLD={params.maf} \
               QUAL_THRESHOLD={params.qual} \
               SAMPLESHEET={input.samplesheet} \
               DONOR={wildcards.donor} \
               SEQUENCING={wildcards.sequencing}
        duckdb -readonly {input.db} -c ".read workflow/scripts/finalize_variants.sql" > {output}
        '''


rule extract_variant_sites:
    input:
        'results/{species}/aligned_pseudogenomes/{sequencing}/{donor}.fas'
    output:
        'results/{species}/snpsites/{sequencing}/{donor}.filtered_alignment.fas',
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
        'results/{species}/snpsites/{sequencing}/{donor}.filtered_alignment.fas',
    params:
        extra='-double-precision -nt'
    output:
        'results/{species}/veryfasttree/{sequencing}/{donor}_5000.veryfasttree.phylogeny.nhx',
    resources:
        cpus_per_task=32,
        mem_mb=64_000,
        runtime=30,
    envmodules:
        'veryfasttree/4.0.3.1'
    shell:
        '''
        export OMP_PLACES=threads
        veryfasttree {input} {params.extra} -threads {resources.cpus_per_task} > {output}
        '''


# TODO: Add outgroup here

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
        'results/{species}/snpsites/{sequencing}/{donor}.filtered_alignment.fas',
    params:
        extra='--all --model GTR+G --bs-trees 200',
        prefix='results/{species}/raxml_ng/{sequencing}/{donor}_5000',
    output:
        multiext(
            'results/{species}/raxml_ng/{sequencing}/{donor}_5000.raxml',
            '.reduced.phy',
            '.rba',
            '.bestTreeCollapsed',
            '.bestTree',
            '.mlTrees',
            '.support',
            '.bestModel',
            '.bootstraps',
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
        raxml-ng {params.extra} --msa {input} --threads {resources.cpus_per_task} --prefix {params.prefix} --redo
        # `reduced.phy` may not be generated after SNP-Sites' `raxml.mlTrees`, etc. not generated unless bootstraping
        touch {output}
        '''


rule:
    input:
        expand(
            [
                'results/{{species}}/veryfasttree/{sequencing}/{donor}_5000.veryfasttree.phylogeny.nhx',
                'results/{{species}}/raxml_ng/{sequencing}/{donor}_5000.raxml.bestTree',
            ],
            sequencing=config['wildcards']['sequencing'].split('|'),
            donor=config['wildcards']['donors'].split('|'),
        )
    output:
        touch('results/{species}/phylogenies_5000.done')
    localrule: True