set memory_limit = getenv('MEMORY_LIMIT');

-- Single-threading ensures deterministict sampling, given seed
set threads = 1;

-- Ensure stdout not polluted
set enable_progress_bar = false;

copy (
    -- Create FASTA entries
    select '>' || "sample" || chr(10) || string_agg(ifnull(allele, '-'), '')
    from (
        select 
            template.sample
            , selected.chrom
            , selected.chrom_pos
            -- Allele is a SNV (alt), else ref
            , coalesce(variant_data.alt, selected.ref) as allele
        from (
            with unpruned as (
                -- Create a random sample of SNV positions, excluding indels
                select 
                    distinct on (chrom, chrom_pos)
                    chrom, chrom_pos, ref
                    , lag(chrom_pos) over (partition by chrom order by chrom_pos) as last_pos
                    , lead(chrom_pos) over (partition by chrom order by chrom_pos) as next_pos
                from candidate_variant_tbl
                where 
                    -- ~ * ~ * ~ * ~
                    -- TODO: Keep indels
                    strlen(ref) = 1 
                    and strlen(alt) = 1 
                    -- ~ * ~ * ~ * ~
                -- TODO: Fix hard-coding of sample size and seed (`getenv(.)` fails)
                using sample reservoir(.01%) repeatable (10023)
            )
            select chrom, chrom_pos, ref
            from unpruned
            -- Filter SNVs by inter-position distance
            where
                (last_pos is null or chrom_pos - last_pos > cast(getenv('INTERBASE_DISTANCE_THRESHOLD') as ubigint))
                and (next_pos is null or next_pos - chrom_pos > cast(getenv('INTERBASE_DISTANCE_THRESHOLD') as ubigint))
        ) selected
        cross join (
            -- Subset using sample info
            select distinct "sample"
            from candidate_variant_tbl
            where 
                "sample" in (
                    select cast(cast("sample" as varchar) as samples)
                    from read_csv(getenv('SAMPLESHEET'))
                    where
                        "group" = getenv('SEQUENCING')
                        and "donor" = getenv('DONOR')
                )
        ) template
        left join (
            select "sample", chrom, chrom_pos, alt
            from candidate_variant_tbl
            where 
                -- ~ * ~ * ~ * ~
                -- TODO: Keep indels
                strlen(ref) = 1 
                and strlen(alt) = 1 
                -- ~ * ~ * ~ * ~
                and alt_fwd_dp >= cast(getenv('ALT_STRAND_DP_THRESHOLD') as usmallint)
                and alt_rev_dp >= cast(getenv('ALT_STRAND_DP_THRESHOLD') as usmallint)
                and dp >= cast(getenv('COVERAGE_THRESHOLD') as usmallint)
                and maf >= cast(getenv('MAF_THRESHOLD') as double) 
                and qual >= cast(getenv('QUAL_THRESHOLD') as usmallint) 
        ) variant_data
        on template.sample = variant_data.sample
            and selected.chrom = variant_data.chrom
            and selected.chrom_pos = variant_data.chrom_pos
        order by selected.chrom, selected.chrom_pos
    )
    group by "sample"
    order by cast("sample" as ubigint)
-- Pipe MSA FASTA to stdout
) to '/dev/stdout' (delimiter '', header false, quote '');