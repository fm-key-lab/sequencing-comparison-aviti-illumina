.bail on

set enable_progress_bar = true;
set memory_limit = getenv('MEMORY_LIMIT');
set preserve_insertion_order = false;
set threads to getenv('SLURM_CPUS_PER_TASK');

create table fastp_read_counts as
select * 
from read_csv(
    '/dev/stdin',
    header = false,
    columns = {
        'sample': 'varchar',
        'before_filtering_total_reads': 'bigint',
        'after_filtering_total_reads': 'bigint',
    },
    auto_detect = false
)