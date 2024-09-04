use rule gubbins from widevariant as remove_recombination with:
    input:
        'results/{group}/{donor}/{species}/pseudogenomes/aligned_pseudogenomes.fas',
    params:
        f=config['gubbins']['filter_percentage'],
        tree_args=config['gubbins']['tree_args'],
        t=config['gubbins']['tree_builder']
    output:
        'results/{group}/{donor}/{species}/gubbins/prefix.final_tree.tre',
    envmodules:
        'intel/21.2.0',
        'impi/2021.2',
        'gubbins/3.3.5'