-- ml duckdb/nightly
-- export MEMORY_LIMIT='100G' NCORES=32 QUAL=30 STRAND_DP=3 DP=8 MAF=".95"
-- duckdb

set memory_limit = getenv('MEMORY_LIMIT');
set threads = getenv('SLURM_CPUS_PER_TASK');
set enable_progress_bar = true;

attach 'results/results.duckdb' as sinfo (read_only);
attach 'results/candidate_variants.duckdb' as cmt (read_only);

-- 1. Filter samples
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
create temp table filtered_variants as
with filtered_calls as (
	select 
		v.sample
		, v.chromosome
		, v.position
		, v.reference
		, coalesce(v.alternate[1], v.reference) as allele 
	from cmt.variants v
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

-- X. Create fasta

-- Ensure stdout not polluted
set enable_progress_bar = false;






-- 1. Filter variants
create temp table variants as
select c.sample, c.position
from cmt.variants c 
where c.chromosome = 'NC_002695.2'
  and c.quality >= cast(getenv('QUAL') as float)
group by c.sample, c.position
having sum(c.depth) >= cast(getenv('DP') as int)
   and min(c.depth) >= cast(getenv('STRAND_DP') as int);

-- 2. Filter samples
create temp table selected_samples as
select distinct(s.sample) 
from sinfo.samplesheet s
where s.donor = 'B002'
  and s.sample in (
    select "sample"
    from sinfo.sequence_typing_results
    where st = 73
);

-- Update variants
create temp table variants as 
select c.* from cmt.variants c
where c.sample in (
    select "sample" from selected_samples
);


-- 2. Filter positions
create temp table selected_positions as

select * 
from cmt.variants c 
where c.chromosome = 'NC_002695.2'
  and c.quality >= cast(getenv('QUAL') as float)
  and c.sample in (select "sample" from selected_samples);

select distinct on(c.position)
    c.sample
    , c.position
    , list(distinct c.allele) as allele
from cmt.variants c 
where c.chromosome = 'NC_002695.2'
  and c.quality >= cast(getenv('QUAL') as float)
  and c.sample in (select "sample" from selected_samples)
group by c.sample, c.position
having sum(c.depth) >= cast(getenv('DP') as int)
   and min(c.depth) >= cast(getenv('STRAND_DP') as int);






create temp table selected_calls as
select c.sample, c.position
    , string_agg(distinct c.allele, ',') as allele
from cmt.variants c 
where c.chromosome = 'NC_002695.2'
  and c.quality >= cast(getenv('QUAL') as float)
  and c.sample in (select "sample" from selected_samples)
group by c.sample, c.position
having sum(c.depth) >= cast(getenv('DP') as int)
   and min(c.depth) >= cast(getenv('STRAND_DP') as int);

-- 3. Filter out invariant
with selected_calls_tmp as (
    select 
        "sample"
        , "position"
        , unnest(allele)
)

select c.*
from cmt.variants c 
where (c.sample, c.position) in (
    select ("sample", "position") from selected_calls
);




-- 2. Filter out invariant sites
create temp table variable_sites as
select 
    distinct on(c.chromosome, c.position) c.chromosome, c.position
from cmt.variants c 
group by c.chromosome, c.position
having count(distinct c.allele) > 1;


-- 2. Filter positions
create temp table selected_calls as
select c.sample, c.position
from cmt.variants c 
where c.chromosome = 'NC_002695.2'
  and c.quality >= cast(getenv('QUAL') as float)
  and c.sample in (
    select "sample" from selected_samples
)
group by c.sample, c.position
having sum(depth) >= cast(getenv('DP') as int)
   and min(depth) >= cast(getenv('STRAND_DP') as int)
;

-- 3. Update variants
with filtered_variants as (
    select "sample", "position", allele
    from cmt.variants
    where ("sample", "position") in (
        select ("sample", "position") from selected_calls
    )
),
invariant_positions as (
    select "position" from filtered_variants
    group by "position"
    having count(distinct allele) > 1
    and count(distinct "sample") > (
        select count(distinct "sample") / 2 from filtered_variants
    )
)


-- 3. Remove invariant positions
select "sample", "position"
from (
    select c.sample, c.position, c.allele
    from selected_calls s
    inner join cmt.variants c 
        on s.sample = c.sample
       and s.position = c.position
)
group by "position"
having count(distinct allele) > 1
   and count(distinct "sample") > (
    select count(distinct "sample") / 2 from selected_calls
);



create temp table variable_positions as
select distinct(c.position)
from cmt.variants c 
where c.position in (
    select "position" from selected_positions
)
  and c.sample in (
    select "sample" from selected_samples
)
group by c.position
having count(distinct c.allele) > 1
   and count(distinct c.sample) > (
    select count(distinct "sample") / 2 from selected_samples
);

-- 4. Filter variants table
create temp table filtered_variants as
select c.* from cmt.variants c 
where c.position in (
    select "position" from variable_positions
)
  and c.sample in (
    select "sample" from selected_samples
);




-- X. Create fasta

-- Ensure stdout not polluted
set enable_progress_bar = false;











set memory_limit = getenv('MEMORY_LIMIT');
set threads = getenv('SLURM_CPUS_PER_TASK');

-- Ensure stdout not polluted
set enable_progress_bar = false;

attach 'results/results.duckdb' as sinfo (read_only);
attach 'results/Escherichia_coli/all_variants.duckdb' as cmt (read_only);

-- 1. Filter variants by sequencing stats.
drop table if exists seq_filtered_cmt;
create table seq_filtered_cmt as
select 
    c.sample
    , c.chrom_pos
    , c.ref
    , coalesce(c.alt, c.ref) as allele
from cmt.candidate_variant_tbl c
-- NOTE: Limiting to this "chromosome"
where c.chrom = 'NC_002695.2'
  and c.qual >= cast(getenv('QUAL') as float)
  and c.dp >= cast(getenv('DP') as int)
  and greatest(ref_af, alt_af) >= cast(getenv('MAF') as float)
  and greatest(least(c.ref_fwd_dp, c.ref_rev_dp), least(c.alt_fwd_dp, c.alt_rev_dp)) >= cast(getenv('STRAND_DP') as int)
;

-- 2. Filter samples by sequence type.
drop table if exists st_seq_filtered_cmt;
create table st_seq_filtered_cmt as
with styping as (
    select "sample", st
    from sinfo.sequence_typing_results
    where st = 73
)
select c.*
from styping s
inner join seq_filtered_cmt c
on c.sample = s.sample;

-- 3. Parse donor metadata
create type groups as enum ('aviti', 'illumina'); -- NOTE: Unused when indels

create type ids as enum (
    select distinct(ID) from sinfo.samplesheet
);

drop table if exists idkey;
create table idkey as
select 
    try_cast("sample" as cmt.samples) as "sample"
    , try_cast("group" as groups) as "group"
    , try_cast(ID as ids) as ID
from sinfo.samplesheet;

-- 4. Add donor metadata
drop table if exists annot_filtered_cmt;
create table annot_filtered_cmt as
select
    s.group
    , s.ID
    , c.chrom_pos
    , c.ref
    , c.allele
from idkey s
inner join st_seq_filtered_cmt c
on c.sample = s.sample;

-- 5. Select variable sites with read support
drop table if exists selected_positions;
create table selected_positions as
with variable_sites as (
    select "group", chrom_pos, count(allele) as n
    from annot_filtered_cmt
    where allele is not null
    group by "group", chrom_pos
    having count(distinct allele) > 1
)
select chrom_pos 
from variable_sites 
-- TODO: Deprecate hard-coding of median (env var + add logic here)
where n >= 90 
group by chrom_pos 
having count(*) = 2;

-- 6. Create (cartesian) template
drop table if exists cartesian_template;
create table cartesian_template as
select * from (select distinct(ID) from annot_filtered_cmt)
cross join (select distinct(chrom_pos) from selected_positions)
cross join (select distinct("group") from annot_filtered_cmt);

-- 7. Parse output
drop table if exists fasta_output;
create table fasta_output as
select 
    template.group
    , template.ID
    , template.chrom_pos
    , geno.allele
from cartesian_template template
left join annot_filtered_cmt geno
on template.group = geno.group
and template.ID = geno.ID
and template.chrom_pos = geno.chrom_pos;

-- remove sample without
drop table if exists final_ids;
create table final_ids as
with minimum_samples as (
    select "group", ID, count(allele) as n 
    from fasta_output 
    group by "group", ID 
    -- TODO: Deprecate hard-coding of heuristic (env var + add logic here)
    having n > 50
)
select ID from minimum_samples group by ID having count("group") = 2;

drop table if exists final_fasta_output;
create table final_fasta_output as
select * from fasta_output 
where ID in (select ID from final_ids)
order by "group", ID, chrom_pos;

-- TODO: Error checking here
select 
    strlen(string_agg(ifnull(allele, 'N'), '')) as n
from (select * from final_fasta_output order by ID, chrom_pos)
group by "group", ID
order by n;

-- 8. Save combined output
copy (
    select 
        -- TODO: Refine imputation of missing sites
        '>' || ID || '_' || "group" || chr(10) || string_agg(ifnull(allele, 'N'), '')
    from final_fasta_output
    group by "group", ID
    order by ID
) to '/dev/stdout' (delimiter '', header false, quote '');