set memory_limit = getenv('MEMORY_LIMIT');
set threads = getenv('SLURM_CPUS_PER_TASK');

copy (
    select '>' || "sample" || chr(10) || string_agg(ifnull(allele, '-'), '')
    from (
        select 
            template.sample
            , selected.chrom
            , selected.chrom_pos
            , coalesce(variant_data.alt, selected.ref) as allele
        from (
            select distinct on (chrom, chrom_pos)
                chrom, chrom_pos, ref
            from candidate_variant_tbl
            using sample reservoir(5%) repeatable (10023)
        ) selected
        cross join (
            select distinct "sample"
            from candidate_variant_tbl
            where 
                "sample" in (
                    select cast(cast("sample" as varchar) as samples)
                    from read_csv(getenv('SAMPLESHEET'))
                    where
                        "group" = getenv('SEQUENCING')
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
) to '/dev/stdout' (delimiter '', header false, quote '');