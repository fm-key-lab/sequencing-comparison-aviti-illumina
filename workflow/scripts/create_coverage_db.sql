.bail on

-- export MEMORY_LIMIT='14G' BEDTOOLS_GENOMECOV="'results/Escherichia_coli/samtools/*_genomecov.tsv'"

set memory_limit = getenv('MEMORY_LIMIT');
set threads = getenv('SLURM_CPUS_PER_TASK');

create type strands as enum ('+', '-');

select *
from read_csv(
	getenv('BEDTOOLS_GENOMECOV'),
	header = false,
	delim = '\t',
	columns = {
		'chrom': 'varchar',
		'chrom_start': 'ubigint',
		'chrom_end': 'ubigint',
		'depth_of_coverage': 'ubigint',
		'strand': strands
	},
	nullstr = '.',
	auto_detect = false
);