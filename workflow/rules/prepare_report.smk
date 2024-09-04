FIGURES = [
    'total_reads',
    'sequence_typing',
    'variant_quality_scores',
]


rule results_notebook:
    input:
        'results/results.duckdb'
    output:
        expand(
            'report/figures/{figure}.{ext}',
            figure=FIGURES,
            ext=['pdf', 'png']
        ),
        'report/figures/sequence_typing_table.tex',
    log:
        notebook='logs/notebooks/AVITI_Illumina_comparison.ipynb'
    notebook:
        'notebooks/AVITI_Illumina_comparison.ipynb'


rule create_report:
    input:
        expand(
            'report/figures/{figure}.{ext}',
            figure=FIGURES,
            ext=['pdf', 'png']
        ),
        'report/figures/sequence_typing_table.tex',
    output:
        'report/main.pdf'
    localrule: True
    envmodules:
        'tinytex/2024.07.03'
    shell:
        'pdflatex -output-directory=report report/main.tex'