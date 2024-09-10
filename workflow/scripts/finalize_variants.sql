set memory_limit = getenv('MEMORY_LIMIT');
set threads = getenv('SLURM_CPUS_PER_TASK');

-- TODO: Must be a bit more dynamic, e.g., output fastas using sample info

-- TODO: SQL logic here can be simplified
create temp table phylogeny_input_candidates as
select 
    af.sample
    -- Collapse chromosome and position into single identifier
    , af.chrom || '_' || af.chrom_pos as chrom_chrom_pos
    , af.ref
    , af.alt
from allele_freq af, (
    select 
        chrom, chrom_pos, ref
    from (
        select 
            count("sample") as num_samples
            , chrom, chrom_pos, ref
        from (
            select *
            from allele_freq
            where
                -- ~ * ~ * ~ * ~
                -- TODO: Keep indels
                strlen(ref) = 1 and 
                strlen(alt) = 1 and
                -- ~ * ~ * ~ * ~
                maf >= cast(getenv('MAF_THRESHOLD') as double) and
                qual >= cast(getenv('QUAL_THRESHOLD') as usmallint)
        )
        group by chrom, chrom_pos, ref
    )
    -- TODO: How to pass variable here? `getenv(.)` doesn't work.
    using sample 100
) rand
where
    af.chrom = rand.chrom and
    af.chrom_pos = rand.chrom_pos
;

copy (
    with wide_tmp as (
        pivot phylogeny_input_candidates
        on chrom_chrom_pos
        using any_value(alt)
        group by "sample"
    )
    select '>' || "sample" || chr(10) || string_agg(ifnull(allele, '-'), '') as fasta_entry
    from wide_tmp unpivot include nulls (
        allele
        for chrom_chrom_pos in (columns(* exclude("sample")))
    )
    group by "sample"
    order by cast("sample" as ubigint)
) to '/dev/stdout' (delimiter '', header false, quote '');