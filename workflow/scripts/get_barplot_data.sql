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
--  $ export MEMORY_LIMIT='100G' NCORES=32 QUAL=30 STRAND_DP=3 DP=8 MAF=".95" SPECIES=Escherichia_coli POSITION=1586465
--  $ mkdir -p /u/thosi/dev/tmp/barplots
--  $ duckdb -c ".read workflow/scripts/get_barplot_data.sql" > /u/thosi/dev/tmp/barplots/position=POSITION.csv

set memory_limit = getenv('MEMORY_LIMIT');
set threads = getenv('SLURM_CPUS_PER_TASK');
set enable_progress_bar = true;

attach 'results/results.duckdb' as sinfo (read_only);
attach 'results/candidate_variants.duckdb' as cmt (read_only);


-- 1. Filter samples
.print Filtering samples...

create temp table selected_samples as
select distinct(s.sample) 
from sinfo.samplesheet s
where s.species = getenv('SPECIES')
  and s.donor = 'B002'
  and s.sample in (
    select "sample"
    from sinfo.sequence_typing_results
    where st = 73
);

-- 2. Filter calls
.print Filtering genotype calls...

create temp table filtered_variants as
with filtered_calls as (
	select 
		v.sample
		, v.chromosome
		, v.position
		, v.reference
		, coalesce(v.alternate[1], v.reference) as allele 
	from variants v
	where v.chromosome = 'NC_002695.2'
	  and v.sample in (
		select "sample" from selected_samples
	  )
	  and v.quality >= cast(getenv('QUAL') as float)
	  and array_reduce(
	  	  v.info_DP4, (x, y) -> x + y
	  ) >= cast(getenv('DP') as int)
	  and list_aggregate(v.info_ADF, 'max') >= cast(getenv('STRAND_DP') as int)
	  and list_aggregate(v.info_ADR, 'max') >= cast(getenv('STRAND_DP') as int)
	  and abs(
		  array_reduce(
			  v.info_DP4[3:4], (x, y) -> x + y
		  ) / array_reduce(
			  v.info_DP4, (x, y) -> x + y
		  ) - .5
	  ) >= (cast(getenv('MAF') as float) / 2)
),
variable_sites as (
	select chromosome, "position"
	from filtered_calls
	group by chromosome, "position"
	having count(distinct allele) > 1
       and count(distinct "sample") > (
        select count(distinct "sample") / 2 from filtered_calls
       )
)
select v.* 
from variable_sites s
inner join filtered_calls v
on s.chromosome = v.chromosome
and s.position = v.position;