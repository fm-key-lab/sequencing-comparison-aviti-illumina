-- Params:
--     MEMORY_LIMIT
--     SLURM_CPUS_PER_TASK
--     COVERAGE_THRESHOLD=2
--     RESOLUTION="-3"
--     SAMPLE_AVG_COVERAGE_THRESHOLD="1.5"
--     N_PERC_THRESHOLD=".1"
--     MAF_THRESHOLD=".85"
--     QUAL_THRESHOLD=30
--     N_POSITIONS=100
--     FASTA_OUTPUT="'test.fasta'"

-- Config
set memory_limit = getenv('MEMORY_LIMIT');
set threads = getenv('SLURM_CPUS_PER_TASK');

-- Create convenience table for removing low-coverage positions
create temp table low_coverage as
select * from coverage_total 
where coverage <= cast(getenv('COVERAGE_THRESHOLD') as usmallint);

create temp table coverage_bed as
select * 
from (
    select 
        chrom
        , chrom_start
        , lead(chrom_start, 1) over (partition by chrom order by chrom_start) as chrom_end
    from (
        with unpivot_tmp as (
            unpivot low_coverage
            on chrom_start, chrom_end
            into
                name feature
                value chrom_pos
        )
        select
            distinct on(chrom, chrom_start)
            chrom
            , cast(round(chrom_pos, cast(getenv('RESOLUTION') as smallint)) as ubigint) as chrom_start
        from unpivot_tmp
    )
    order by chrom, chrom_start
)
where chrom_end is not null;

create temp table coverage_blacklist as
select * from (
    select
        count(c.coverage) as count_low_coverage_samples
        , b.chrom
        , b.chrom_start
        , b.chrom_end
    from low_coverage c, coverage_bed b
    where
        b.chrom_start <= c.chrom_start and
        b.chrom_end >= c.chrom_end and
        b.chrom = c.chrom
    group by 
        "sample"
        , b.chrom
        , b.chrom_start
        , b.chrom_end
)
where 
    -- NOTE: Equivalent to a median coverage filter
    count_low_coverage_samples >= (
        select count(distinct "sample") / 2 from low_coverage
    )
order by count_low_coverage_samples;

-- TODO: Filter by blacklist (empty for current params)

-- Randomly sample variants
create temp table phylogeny_input_candidates as
select 
    af.sample
    -- Collapse chromosome and position into single unique identifier
    , af.chrom || '_' || af.chrom_pos as chrom_chrom_pos
    , af.ref
    , af.alt
from allele_freq af, (
    select 
        chrom
        , chrom_pos
        , ref
    from (
        select 
            count("sample") as num_samples
            , chrom
            , chrom_pos
            , ref
        from (
            select
                "sample"
                , chrom
                , chrom_pos
                , ref
                , alt
                -- Use AF == (DP4[2] + DP4[3]) / DP' where DP' /neq DP, but is sum(DP4)
                -- e.g., https://workflowhub.eu/workflows/354
                , (
                    cast(split_part(info_dp4, ',', 3) as usmallint) 
                    + cast(split_part(info_dp4, ',', 4) as usmallint)
                ) / array_reduce(
                    array_transform(
                        regexp_split_to_array(info_dp4, ','), 
                        x -> cast(x as usmallint)
                    ), 
                    (x, y) -> x + y
                ) as maf
            from allele_freq
            where
                "sample" in (
                    select "sample"
                    from coverage_sample_avg 
                    where coverage_sample >= cast(getenv('SAMPLE_AVG_COVERAGE_THRESHOLD') as double)
                ) and
                "sample" in (
                    select "sample"
                    from base_freq 
                    where 
                        base = 'N' and 
                        avg_perc < cast(getenv('N_PERC_THRESHOLD') as double)
                ) and
                maf >= cast(getenv('MAF_THRESHOLD') as double) and
                qual >= cast(getenv('QUAL_THRESHOLD') as usmallint)
        )
        group by chrom, chrom_pos, ref
    )
    -- NOTE: Source of bias: Implicit pre-down-sampling to have no private mutations
    where num_samples > 1
    -- NOTE: Does not work with `getenv('N_POSITIONS')`
    using sample 1000
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
-- NOTE: Does not work with `getenv('FASTA_OUTPUT')`
) to '/u/thosi/dev/test2.fasta' (delimiter '', header false, quote '');