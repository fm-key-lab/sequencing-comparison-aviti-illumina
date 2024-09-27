set memory_limit = getenv('MEMORY_LIMIT');
set threads = getenv('SLURM_CPUS_PER_TASK');

create type chroms as enum ('NC_002695.2', 'NC_002128.1', 'NC_002127.1');
create type strands as enum ('+', '-');

create table genome_coverage_bed as
select
	regexp_extract(
        "filename",
        'samtools/(\d+)_genomecov.tsv',
        1
    ) as "sample"
    , * exclude("filename")
from read_csv(
	getenv('BEDTOOLS_GENOMECOV'),
	header = false,
	delim = '\t',
	columns = {
		'chrom': chroms,
		'chrom_start': 'uinteger',
		'chrom_end': 'uinteger',
		'depth_of_coverage': 'usmallint',
		'strand': strands
	},
    filename = true,
	nullstr = '.',
	auto_detect = false
);

-- TODO: Add metadata on calling, code, etc.

create type samples as enum (
    select distinct("sample") from genome_coverage_bed
);

alter table genome_coverage_bed alter "sample" type samples;