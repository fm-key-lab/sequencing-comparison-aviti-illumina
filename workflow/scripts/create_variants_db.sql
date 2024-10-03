-- export MEMORY_LIMIT="100GB" VCFS='/ptmp/thosi/sequencing-comparison-aviti-illumina/results/lake/*/*65*/nonindels.parquet'

set memory_limit = getenv('MEMORY_LIMIT');
set threads = getenv('SLURM_CPUS_PER_TASK');

-- Create ENUM types

create type bases as enum ('A', 'C', 'T', 'G', 'N'); -- NOTE: Unused when indels

create type chroms as enum ('NC_002695.2', 'NC_002128.1', 'NC_002127.1');

create type samples as enum (
    select distinct(regexp_extract("file", '/sample=(\d+)/', 1))
	from glob(getenv('VCFS'))
);

create type taxon as enum ('Escherichia_coli');

create type orient as enum ('forward', 'reverse');

-- Parse VCF files

create temp table variants as
with unnested_tmp as (
	select	
		species, "sample", chromosome, "position"
		, unnest([reference, alternate]) as allele
		, unnest(info_ADF) as forward
		, unnest(info_ADR) as reverse
		, case when unnest(['reference', 'alternate']) = 'alternate' then true else false end as is_alternate
	from (
		select 
			species
			, "sample"
			, cast(chromosome as chroms) as chromosome
			, cast("position" as ubigint) as "position"
			, cast(reference as bases) as reference
			, cast(nullif(unnest(alternate), '') as bases) as alternate
			, quality
			, info_ADF
			, info_ADR
		from read_parquet(
			getenv('VCFS'), 
			hive_partitioning = true,
			hive_types = {
				'species': taxon,
				'sample': samples
			}
		)
	)
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
select * exclude(read_orientation, depth)
	, cast(read_orientation as orient) as read_orientation
	, cast(depth as usmallint) as depth
from unpivoted_tmp
where depth > 0;