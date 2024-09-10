set memory_limit = getenv('MEMORY_LIMIT');
set threads = getenv('SLURM_CPUS_PER_TASK');

create type chroms as enum ('NC_002695.2', 'NC_002128.1', 'NC_002127.1');

create table allele_freq as 
select 
	"sample"
	, cast(chrom as chroms) as chrom
	, cast(chrom_pos as ubigint) as chrom_pos
	, ref
	, alt
	, (
		cast(split_part(info_dp4, ',', 3) as usmallint) + 
		cast(split_part(info_dp4, ',', 4) as usmallint)
	) / array_reduce(
		array_transform(
			regexp_split_to_array(info_dp4, ','), 
			x -> cast(x as usmallint)
		), 
		(x, y) -> x + y
	) as maf
	, cast(qual as decimal(4, 1)) as qual
	, dp
	, info_dp4
from read_csv(
	getenv('BCFTOOLS_QUERY'),
	header = false,
	delim = '\t',
	columns = {
		'chrom': 'varchar',
		'chrom_pos': 'bigint',
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
where strlen(alt) >= 1;

-- NOTE: Since only dealing with SNVs
create type bases as enum ('A', 'C', 'T', 'G', 'N');
alter table allele_freq alter ref type bases;
alter table allele_freq alter alt type bases;

create type samples as enum (
    select distinct("sample") from allele_freq
);

alter table allele_freq alter "sample" type samples;

-- TODO: Add metadata.