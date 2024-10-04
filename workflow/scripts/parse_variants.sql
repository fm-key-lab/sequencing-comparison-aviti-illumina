-- Environment variable arguments:
--  MEMORY_LIMIT='100G'
--  SLURM_CPUS_PER_TASK=32
--  QUAL=30
--  STRAND_DP=3
--  DP=8
--  MAF=".95"
--
-- Usage:
--  $ export export MEMORY_LIMIT='100G' NCORES=32 QUAL=30 STRAND_DP=10 DP=8 MAF=".95"
--  $ duckdb results/candidate_variants.duckdb -c ".read workflow/scripts/parse_variants.sql" > /u/thosi/dev/tmp/variants2.csv

set memory_limit = getenv('MEMORY_LIMIT');
set threads = getenv('SLURM_CPUS_PER_TASK');
set enable_progress_bar = true;

attach 'results/results.duckdb' as sinfo (read_only);

-- 1. Filter samples
.print Filtering samples...

create temp table selected_samples as
select distinct(s.sample) 
from sinfo.samplesheet s
where s.donor = 'B002'
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

-- 3. Annotate samples
.print Annotating sample data...

create temp table annotated_variants as
with cartesian_tmp as (
    select * from (
        select distinct on(chromosome, "position") 
            chromosome, "position", reference
        from filtered_variants
    )
    cross join (
        select distinct("sample")
        from filtered_variants
    )
),
annot_cartesian_tmp as (
    select s.group, s.ID, v.sample, v.chromosome, v.position, v.reference
    from sinfo.samplesheet s
    inner join cartesian_tmp v
    on s.sample = v.sample
)
select c.*, v.allele
from annot_cartesian_tmp c
left join filtered_variants v
on c.sample = v.sample
and c.chromosome = v.chromosome
and c.position = v.position
and c.reference = v.reference;

-- 4. Final cleaning
.print Preparing output...

create temp table parsed_variant_tbl as
select * exclude(allele)
    , ifnull(allele, 'N') as allele
from annotated_variants
where ID in (
    select ID from (
        select distinct on ("group", ID) "group", ID
        from annotated_variants
        where "sample" in (
            select "sample"
            from annotated_variants
            group by "sample"
            having count(allele) > (
                -- TODO:
                select count(distinct "position") * .8 from annotated_variants
            )
        )
    )
    group by "ID"
    having count(ID) = 2
);

-- 5. Output
.print Returning variants.

copy parsed_variant_tbl to '/dev/stdout' (format csv);