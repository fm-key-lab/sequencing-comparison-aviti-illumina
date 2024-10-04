-- Environment variable arguments:
--  MEMORY_LIMIT='100G'
--  SLURM_CPUS_PER_TASK=32
--  QUAL=30
--  STRAND_DP=3
--  DP=8
--  MAF=".95"
--  POSITION=1586465
--
-- Usage:
--  $ export MEMORY_LIMIT='100G' NCORES=32 QUAL=30 STRAND_DP=3 DP=8 MAF=".95" POSITION=1586465
--  $ duckdb -c ".read workflow/scripts/parse_variants.sql" > /u/thosi/dev/tmp/variants2.csv

set memory_limit = getenv('MEMORY_LIMIT');
set threads = getenv('SLURM_CPUS_PER_TASK');
set enable_progress_bar = true;

attach 'results/results.duckdb' as sinfo (read_only);
attach 'results/candidate_variants.duckdb' as cmt (read_only);

select * from 