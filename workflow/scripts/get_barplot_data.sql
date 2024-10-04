-- Environment variable arguments:
--  MEMORY_LIMIT='100G'
--  SLURM_CPUS_PER_TASK=32
--  QUAL=30
--  STRAND_DP=3
--  DP=8
--  MAF=".95"
--  POSITION=1586465
--  SEQTYPE=73
--  DONOR=B002
--
-- Usage:
--  $ export MEMORY_LIMIT='100G' 
--  $ export QUAL=30 STRAND_DP=3 DP=8 MAF=".95"
--  $ export SPECIES=Escherichia_coli POSITION=1586465 SEQTYPE=73 DONOR=B002
--  $ mkdir -p /u/thosi/dev/tmp/barplots
--  $ duckdb -c ".read workflow/scripts/get_barplot_data.sql" > /u/thosi/dev/tmp/barplots/position=POSITION.csv

set memory_limit = getenv('MEMORY_LIMIT');
set threads = getenv('SLURM_CPUS_PER_TASK');
set enable_progress_bar = true;

attach 'results/results.duckdb' as sinfo (read_only);
attach 'results/candidate_variants.duckdb' as cmt (read_only);

create type orient as enum ('forward', 'reverse');

-- 1. Filter samples
.print Filtering samples...

create temp table selected_samples as
select distinct(s.sample) 
from sinfo.samplesheet s
where s.species = getenv('SPECIES')
  and s.donor = getenv('DONOR')
  and s.sample in (
    select "sample"
    from sinfo.sequence_typing_results
    where st = cast(getenv('SEQTYPE') as int)
);

-- 2. Filter calls
.print Filtering genotype calls...

create temp table filtered_variants as
with unnested_tmp as (
    select "sample", chromosome, "position"
        , quality
        , unnest([reference, alternate[1]]) as allele
        , unnest(info_ADF) as forward
        , unnest(info_ADR) as reverse
        , case when unnest(['reference', 'alternate']) = 'alternate' then true else false end as is_alternate
    from cmt.variants
    where "position" = cast(getenv('POSITION') as int)
    and chromosome = 'NC_002695.2'
    and "sample" in (select "sample" from selected_samples)
    and quality >= cast(getenv('QUAL') as float)
    and array_reduce(info_DP4, (x, y) -> x + y) >= cast(getenv('DP') as int)
    and list_aggregate(info_ADF, 'max') >= cast(getenv('STRAND_DP') as int)
    and list_aggregate(info_ADR, 'max') >= cast(getenv('STRAND_DP') as int)
    and abs(
        array_reduce(info_DP4[3:4], (x, y) -> x + y) / array_reduce(info_DP4, (x, y) -> x + y) - .5
    ) >= (cast(getenv('MAF') as float) / 2)
),
unpivoted_tmp as (
	unpivot (
		select * from unnested_tmp
		where allele is not null
	) 
	on forward, reverse
	into 
		name read_orientation
		value depth	
)
select
	* exclude(read_orientation, depth)
	, cast(read_orientation as orient) as read_orientation
	, cast(depth as usmallint) as depth
from unpivoted_tmp
where depth > 0;

-- 3. Annotate samples
.print Annotating sample data...

create temp table annotated_variants as
select 
    s.group, s.ID, v.sample, v.chromosome, v.position
    , v.quality, v.is_alternate, v.allele, v.read_orientation, v.depth
from sinfo.samplesheet s
inner join filtered_variants v
on s.sample = v.sample;

-- 4. Output
.print Returning variants.

copy annotated_variants to '/dev/stdout' (format csv);