set memory_limit = getenv('MEMORY_LIMIT');
set threads = getenv('SLURM_CPUS_PER_TASK');

-- Base frequencies
create table base_freq as 
select
	regexp_extract(filename, 'base_freq/(\d+).csv', 1) as sample
	, * exclude(filename)
from read_csv(
	getenv('BASE_FREQ'),
	header = false,
	delim = ',',
	columns = {
		'avg_perc': 'double',
		'read_pe': 'bigint',
		'base': 'varchar'
	},
	auto_detect = false,
	filename = true
);

-- Coverage
create table coverage as 
select 
	regexp_extract(filename, 'coverage/(\d+).tsv', 1) as sample
	, * exclude(filename)
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

-- Allele frequencies
create table allele_freq as 
select * from read_csv(
	getenv('ALLELE_FREQ'),
	header = false,
	delim = '\t',
	columns = {
		'chrom': 'varchar',
		'chrom_pos': 'bigint',
		'ref': 'varchar',
		'alt': 'varchar',
		'sample': 'bigint',
		'qual': 'double',
		'info_dp4': 'varchar'
	},
	nullstr = '.',
	auto_detect = false
)
where strlen(alt) >= 1;