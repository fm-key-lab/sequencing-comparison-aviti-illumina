set memory_limit = getenv('MEMORY_LIMIT');
set threads = getenv('SLURM_CPUS_PER_TASK');
set enable_progress_bar = false;

copy (
    select 
        '>' || ID || '_' || "group" || chr(10) || string_agg(ifnull(allele, 'N'), '')
    from read_csv('/dev/stdin')
    group by "group", ID
    order by ID
) to '/dev/stdout' (delimiter '', header false, quote '');"