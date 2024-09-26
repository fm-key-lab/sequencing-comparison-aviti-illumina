.bail on

-- export MEMORY_LIMIT='14G' BEDTOOLS_GENOMECOV="'results/Escherichia_coli/samtools/*_genomecov.tsv'"

set memory_limit = getenv('MEMORY_LIMIT');
set threads = getenv('SLURM_CPUS_PER_TASK');

select *
from read_csv(
	getenv('BEDTOOLS_GENOMECOV'),
	header = false,
	delim = '\t',
	columns = {
		'chrom': 'varchar',
		'chrom_pos': 'ubigint',
		'ref': 'varchar',
		'alt': 'varchar',
		'sample': 'varchar',
		'qual': 'double',
		'dp': 'bigint',
		'info_dp4': 'varchar'
	},
	nullstr = '.',
	auto_detect = false
)