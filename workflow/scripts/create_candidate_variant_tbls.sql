set memory_limit = getenv('MEMORY_LIMIT');
set threads = getenv('SLURM_CPUS_PER_TASK');

-- Create enum types
create type bases as enum ('A', 'C', 'T', 'G', 'N');
create type chroms as enum ('NC_002695.2', 'NC_002128.1', 'NC_002127.1');
create type reads as enum ('1', '2');
create type strands as enum ('-', '+');

-- Base frequencies
create table base_freq as 
select
	regexp_extract(filename, 'base_freq/(\d+).csv', 1) as "sample"
	, cast(base as bases) as base
	, cast(read_pe as reads) as read_pe
	, cast(avg_perc as decimal(3, 2)) as avg_perc
from read_csv(
	getenv('BASE_FREQ'),
	header = false,
	delim = ',',
	columns = {
		'avg_perc': 'double',
		'read_pe': 'varchar',
		'base': 'varchar'
	},
	auto_detect = false,
	filename = true
);

-- Create and use sample enum
create type samples as enum (
    select distinct("sample") from base_freq
);

alter table base_freq alter "sample" type samples;

-- Coverage
create table coverage_strand as 
select 
	cast(regexp_extract(filename, 'coverage/(\d+)_strand.tsv', 1) as samples) as "sample"
	, cast(chrom as chroms) as chrom
	, cast(strand as strands) as strand
	, cast(chrom_start as ubigint) as chrom_start
	, cast(chrom_end as ubigint) as chrom_end
	, cast(coverage as usmallint) as coverage
from read_csv(
	getenv('COVERAGE'),
	header = false,
	delim = '\t',
	nullstr = '',
	columns = {
		'chrom': 'varchar',
		'chrom_start': 'bigint',
		'chrom_end': 'bigint',
		'coverage': 'bigint',
		'strand': 'varchar'
	},
	auto_detect = false,
	filename = true
);

create table coverage_total as 
select 
	cast(regexp_extract(filename, 'coverage/(\d+)_total.tsv', 1) as samples) as "sample"
	, cast(chrom as chroms) as chrom
	, cast(chrom_start as ubigint) as chrom_start
	, cast(chrom_end as ubigint) as chrom_end
	, cast(coverage as usmallint) as coverage
from read_csv(
	getenv('COVERAGE_TOT'),
	header = false,
	delim = '\t',
	nullstr = '',
	columns = {
		'chrom': 'varchar',
		'chrom_start': 'bigint',
		'chrom_end': 'bigint',
		'coverage': 'bigint'
	},
	auto_detect = false,
	filename = true
);

-- Estimate average coverage per sample
create table coverage_sample_avg as
select
	"sample"
	-- In case the covered interals << genome
	, intervals_size
	, cast((coverage_sum / intervals_size) as decimal(4, 1)) as coverage_sample
from (
	select
		"sample"
		, sum(chrom_end - chrom_start) as intervals_size
		, sum((chrom_end - chrom_start) * coverage) as coverage_sum
	from coverage_total
	group by "sample"
);

-- Allele frequencies
create table allele_freq as 
select 
	cast("sample" as samples) as "sample"
	, cast(chrom as chroms) as chrom
	, cast(chrom_pos as ubigint) as chrom_pos
	, cast(ref as bases) as ref
	, cast(alt as bases) as alt
	, cast(qual as decimal(4, 1)) as qual
	, af
	, dp
	, info_dp4
from read_csv(
	getenv('ALLELE_FREQ'),
	header = false,
	delim = '\t',
	columns = {
		'chrom': 'varchar',
		'chrom_pos': 'bigint',
		'ref': 'varchar',
		'alt': 'varchar',
		'sample': 'varchar',
		'qual': 'double',
		'af': 'double',
		'dp': 'bigint',
		'info_dp4': 'varchar'
	},
	nullstr = '.',
	auto_detect = false
)
-- NOTE: Use this filter for keeping indels. 
--		 Must also change type casting.
-- where strlen(alt) >= 1;
where
	strlen(ref) = 1 and 
	strlen(alt) = 1
;