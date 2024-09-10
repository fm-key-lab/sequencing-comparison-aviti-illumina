set memory_limit = getenv('MEMORY_LIMIT');
set threads = getenv('SLURM_CPUS_PER_TASK');

-- TODO: Must be a bit more dynamic, e.g., output fastas using sample info

-- TODO: SQL logic here can be simplified

-- TODO: Rewrite to deprecate create table (and thus allow read concurrency)

create temp table phylogeny_input_candidates as
select 
    variants.sample
    -- Collapse chromosome and position into single identifier
    , variants.chrom || '_' || variants.chrom_pos as chrom_chrom_pos
    , variants.ref
    , variants.alt
from candidate_variant_tbl variants, (
    select 
        chrom, chrom_pos, ref
    from (
        select 
            count("sample") as num_samples
            , chrom, chrom_pos, ref
        from (
            select *
            from candidate_variant_tbl
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
    using sample 1000
) selected
where
    variants.chrom = selected.chrom and
    variants.chrom_pos = selected.chrom_pos
;

copy (
    with wide_tmp as (
        pivot phylogeny_input_candidates
        on chrom_chrom_pos
        using any_value(alt)
        group by "sample"
    )
    select '>' || "sample" || chr(10) || string_agg(ifnull(allele, '-'), '-') as fasta_entry_gapped
    -- select '>' || "sample" || chr(10) || string_agg(ifnull(allele, '-'), '') as fasta_entry
    from wide_tmp unpivot include nulls (
        allele
        for chrom_chrom_pos in (columns(* exclude("sample")))
    )
    group by "sample"
    order by cast("sample" as ubigint)
) to '/dev/stdout' (delimiter '', header false, quote '');