.bail on

set enable_progress_bar = true;
set memory_limit = getenv('MEMORY_LIMIT');
set preserve_insertion_order = false;
set threads to getenv('SLURM_CPUS_PER_TASK');

-- See https://github.com/MultiQC/MultiQC/blob/a70cdd4e60eab452b9fabdd14f035e7483fe3a0b/multiqc/modules/bcftools/stats.py#L192-L202

create table bcftools_variant_quality_scores as
select * 
from read_csv(
    '/dev/stdin',
    header = false,
    columns = {
        'sample': 'varchar',
        'quality': 'bigint',
        'count': 'bigint',
    },
    auto_detect = false
)