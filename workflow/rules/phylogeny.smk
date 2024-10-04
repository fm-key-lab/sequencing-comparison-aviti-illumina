rule pseudogenome_alignment:
    input:
        ancient('results/candidate_variants.duckdb'),
    output:
        'results/aligned_pseudogenomes/{species}.variants.csv',
        'results/aligned_pseudogenomes/{species}.fas',
    params:
        alt=config['variants_thresh']['reads_alt'],
        covg=config['variants_thresh']['coverage'],
        maf=config['variants_thresh']['maf'],
        qual=config['variants_thresh']['qual'],
    resources:
        cpus_per_task=32,
        mem_mb=96_000,
        runtime=15,
    envmodules:
        'duckdb/nightly'
    shell:
        '''
        export MEMORY_LIMIT="$(({resources.mem_mb} / 1100))GB" \
               DP={params.covg} \
               MAF={params.maf} \
               QUAL={params.qual} \
               SPECIES={wildcards.species} \
               STRAND_DP={params.alt}
        
        duckdb -readonly {input} -c ".read workflow/scripts/parse_variants.sql" > {output[0]}

        cat {output[0]} |\
          duckdb -c ".read workflow/scripts/to_msa.sql" > {output[1]}
        '''


rule filter_invariant_sites:
    input:
        'results/aligned_pseudogenomes/{species}.fas',
    output:
        'results/aligned_pseudogenomes/{species}-filtered.fas',
    localrule: True
    envmodules:
        'snp-sites/2.5.1'
    shell:
        '''
        snp-sites -o {output} {input}
        '''


rule raxml_ng:
    """Run RAxML-NG.

    Build a maximum-likelihood phylogeny from the multiple sequence alignment 
    using RAxML Next Generation.

    Args:

    Returns:
      {{ prefix }}.raxml.reduced.phy: Reduced alignment (with duplicates 
        and gap-only sites/taxa removed) [OPTIONAL]
      {{ prefix }}.raxml.rba: Binary MSA file
      {{ prefix }}.raxml.bestTreeCollapsed: Best ML tree with collapsed 
        near-zero branches [OPTIONAL]
      {{ prefix }}.raxml.bestTree: Best ML tree
      {{ prefix }}.raxml.mlTrees: All ML trees
      {{ prefix }}.raxml.support: Best ML tree with Felsenstein bootstrap 
        (FBP) support [OPTIONAL, with bootstrap]
      {{ prefix }}.raxml.bestModel: Optimized model
      {{ prefix }}.raxml.bootstraps: Bootstrap trees [OPTIONAL, with bootstrap]
      {{ prefix }}.raxml.log: Execution log
    
    Notes:
      - [GitHub](https://github.com/amkozlov/raxml-ng)
    """
    input:
        'results/aligned_pseudogenomes/{species}-filtered.fas',
    params:
        extra='--all --model GTR+G --bs-trees 200',
        # outgroup=lambda wildcards: '--outgroup ' + config['outgroup'][wildcards.donor][wildcards.species]['ID'],
        outgroup='',
        prefix='results/raxml_ng/{species}',
    output:
        multiext(
            'results/raxml_ng/{species}.raxml',
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
        raxml-ng \
          {params.extra} \
          {params.outgroup} \
          --msa {input} \
          --threads {resources.cpus_per_task} \
          --prefix {params.prefix} \
          --redo
        touch {output}
        '''


rule:
    input:
        expand(
            'results/raxml_ng/{species}.raxml.bestTree',
            species=config['wildcards']['species'].split('|'),
        )
    output:
        touch('results/phylogenies.done')
    localrule: True