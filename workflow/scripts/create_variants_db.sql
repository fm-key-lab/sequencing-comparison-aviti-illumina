set memory_limit = getenv('MEMORY_LIMIT');
set threads = getenv('SLURM_CPUS_PER_TASK');

create type bases as enum ('A', 'C', 'T', 'G', 'N'); -- NOTE: Unused when indels

create type chroms as enum ('NC_002695.2', 'NC_002128.1', 'NC_002127.1');

create type filters as enum ('PASS', 'LowQual');

create type samples as enum (
    select regexp_extract("file", 'variants/(\d+).filtered.vcf.gz', 1) 
	-- from glob('results/*/variants/*.filtered.vcf.gz')
	from glob(getenv('FILTERED_VCFS'))
);

drop table if exists raw_vcfs;
create temp table raw_vcfs as
select
	cast(regexp_extract(
        filename, 'variants/(\d+).filtered.vcf.gz', 1
    ) as samples) as "sample"
	, "#CHROM" as chrom
	, POS as chrom_pos
	, cast(REF as bases) as REF
	, cast(ALt as bases) as ALT
	, INFO
	, QUAL as qual
	, "FILTER" as "filter"
from read_csv(
	-- 'results/Escherichia_coli/variants/*.filtered.vcf.gz',
	getenv('FILTERED_VCFS'),
	auto_detect = false,
	columns = {
		'#CHROM': chroms,
		'POS': 'ubigint',
		'ID': 'varchar',
		'REF': 'varchar',
		'ALT': 'varchar',
		'QUAL': 'double',
		'FILTER': filters,
		'INFO': 'varchar',
		'FORMAT': 'varchar',
		'VALUES': 'varchar'
	},
	compression = gzip,
	delim = '\t',
    filename = true,
	header = true,
	-- TODO: With 'ignore errors', still need to pass `skip`?
	ignore_errors = true,
	nullstr = '.',
	parallel = true,
	skip = 43
)
-- NOTE: INDELs are filtered out
where INFO not ilike '%INDEL%';

drop table if exists variants;
create table variants as
with unpivot_allele as (
	unpivot raw_vcfs
	on REF, ALT
	into
		name is_ref
		value allele
),
unpivot_strand as (
	unpivot (
		select 
			"sample"
			, chrom
			, chrom_pos
			, qual
			, "filter"
			, is_ref
			, allele
			, regexp_extract(
				INFO,
				'DP4=(\d+),(\d+),(\d+),(\d+)',
				['ref-forward', 'ref-reverse', 'alt-forward', 'alt-reverse']
			) as dp4
			, case 
				when is_ref = 'REF' 
					then cast(dp4['ref-forward'] as usmallint)
				else cast(dp4['alt-forward'] as usmallint)
			end as forward
			, case 
				when is_ref = 'REF' 
					then cast(dp4['ref-reverse'] as usmallint)
				else cast(dp4['alt-reverse'] as usmallint)
			end as reverse
		from unpivot_allele
	)
	on forward, reverse
	into
		name read_orient
		value dp
) 
select
	"sample"
	, chrom
	, chrom_pos
	, qual
	, "filter"
	, case is_ref when 'ref' then true when 'alt' then false end as is_ref
	, read_orient
	, allele
	, dp
from unpivot_strand;




create table candidate_variant_tbl as 
select 
	"sample"
	, cast(chrom as chroms) as chrom
	, chrom_pos
	, ref
	, alt
	, cast(split_part(info_dp4, ',', 1) as usmallint) as ref_fwd_dp
	, cast(split_part(info_dp4, ',', 2) as usmallint) as ref_rev_dp
	, cast(split_part(info_dp4, ',', 3) as usmallint) as alt_fwd_dp
	, cast(split_part(info_dp4, ',', 4) as usmallint) as alt_rev_dp
	, (
		cast(split_part(info_dp4, ',', 1) as usmallint) + 
		cast(split_part(info_dp4, ',', 2) as usmallint)
	) / array_reduce(
		array_transform(
			regexp_split_to_array(info_dp4, ','), 
			x -> cast(x as usmallint)
		), 
		(x, y) -> x + y
	) as ref_af
	, (
		cast(split_part(info_dp4, ',', 3) as usmallint) + 
		cast(split_part(info_dp4, ',', 4) as usmallint)
	) / array_reduce(
		array_transform(
			regexp_split_to_array(info_dp4, ','), 
			x -> cast(x as usmallint)
		), 
		(x, y) -> x + y
	) as alt_af
	, cast(qual as decimal(4, 1)) as qual
	, dp
	, info_dp4
from read_csv(
	getenv('BCFTOOLS_QUERY'),
	header = false,
	delim = '\t',
	columns = {
		'chrom': 'varchar',
		'chrom_pos': 'ubigint',
		'ref': 'varchar',
		'alt': 'varchar',
		'sample': 'varchar',
		'qual': 'double',
		'info': 'varchar',
		'dp': 'bigint',
		'info_dp4': 'varchar'
	},
	nullstr = '.',
	auto_detect = false
)
-- NOTE: INDELs are filtered out
where info not ilike '%INDEL%';

-- TODO: Add metadata on calling, code, etc.

create type samples as enum (
    select distinct("sample") from candidate_variant_tbl
);

alter table candidate_variant_tbl alter "sample" type samples;