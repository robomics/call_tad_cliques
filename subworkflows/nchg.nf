// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

// Workaround for optional input files: https://github.com/nextflow-io/nextflow/issues/1694
def make_optional_input(path) {
    if (path?.trim()) {
        return [file(path)]
    }
    return []
}

workflow NCHG {

    take:
        sample_sheet
        domains
        mad_max
        bad_bin_fraction

        cytoband
        gaps
        custom_mask

    main:

        GENERATE_MASK(
            make_optional_input(cytoband),
            make_optional_input(gaps),
            make_optional_input(custom_mask),
            params.zstd_compression_lvl
        )

        MASK_DOMAINS(
            domains,
            GENERATE_MASK.out.bed,
            params.zstd_compression_lvl
        )

        MASK_DOMAINS.out.set { filtered_domains }

        interaction_types = []
        if (params.call_cis_cliques) {
            interaction_types.add("cis")
        }
        if (params.call_trans_cliques) {
            interaction_types.add("trans")
        }

        sample_sheet
            .splitCsv(sep: "\t", header: true)
            .map { row -> tuple(row.sample,
                                file(row.hic_file, checkIfExists: true),
                                row.resolution)
            }
            .set { hic_files }


        EXPECTED(
            hic_files,
            interaction_types,
            mad_max
        )

        GENERATE_CHROMOSOME_PAIRS(
            hic_files.combine(interaction_types)
        )

        DUMP_CHROM_SIZES(
            hic_files
        )

        GENERATE_CHROMOSOME_PAIRS.out
            .splitCsv(sep: "\t", header: ["sample", "chrom1", "chrom2"])
            .map { tuple(it.sample, it.chrom1, it.chrom2) }
            .combine(hic_files, by: 0)
            .combine(filtered_domains, by: 0)
            .combine(EXPECTED.out.h5, by: 0)
            .set { nchg_compute_tasks }

        COMPUTE(
            nchg_compute_tasks,
            bad_bin_fraction,
        )

        COMPUTE.out.parquet
            .branch {
                cis: it[1] == it[2]
                trans: true
            }
            .set { significant_interactions }

        significant_interactions.cis
            .map { tuple(it[0], it[3]) }
            .groupTuple()
            .combine(["cis"])
            .join(DUMP_CHROM_SIZES.out.tsv)
            .set { nchg_merge_cis_tasks }

        significant_interactions.trans
            .map { tuple(it[0], it[3]) }
            .groupTuple()
            .combine(["trans"])
            .join(DUMP_CHROM_SIZES.out.tsv)
            .set { nchg_merge_trans_tasks }

        nchg_merge_cis_tasks
            .mix(nchg_merge_trans_tasks)
            .set { nchg_merge_tasks }

        MERGE(
            nchg_merge_tasks
        )


        MERGE.out.parquet
            .branch {
                cis: it[1] == "cis"
                trans: true
            }
            .set { nchg_merge_output }

        nchg_merge_output.cis
            .map { tuple(it[0],
                         it[1],
                         it[2],
                         params.nchg_fdr_cis,
                         params.nchg_log_ratio_cis)
            }
            .set { nchg_filter_cis_tasks }

        nchg_merge_output.trans
            .map { tuple(it[0],
                         it[1],
                         it[2],
                         params.nchg_fdr_trans,
                         params.nchg_log_ratio_trans)
            }
            .set { nchg_filter_trans_tasks }

        nchg_filter_cis_tasks
            .mix(nchg_filter_trans_tasks)
            .set { nchg_filter_tasks }

        FILTER(
            nchg_filter_tasks
        )

        VIEW(
            FILTER.out.parquet
        )

        CONCAT(
            VIEW.out.groupTuple()
        )

        if (!params.skip_expected_plots) {
            PLOT_EXPECTED(
                EXPECTED.out.h5,
                interaction_types
            )
        }

        if (!params.skip_sign_interaction_plots) {
            GET_HIC_PLOT_RESOLUTION(
                hic_files,
                params.hic_tgt_resolution_plots
            )

            GET_HIC_PLOT_RESOLUTION.out.tsv
                .splitCsv(header: ["sample", "resolution"],
                          sep: "\t")
                .map { tuple(it.sample, it.resolution) }
                .set { plotting_resolutions }

            GENERATE_CHROMOSOME_PAIRS.out.tsv
                .splitCsv(header: ["sample", "chrom1", "chrom2"],
                          sep: "\t")
                .map { tuple(it.sample, it.chrom1, it.chrom2) }
                .set { chrom_pairs }

            interaction_types = []
            if (params.call_cis_cliques) {
                interaction_types.push("cis")
            }
            if (params.call_trans_cliques) {
                interaction_types.push("trans")
            }

            hic_files
                .map { tuple(it[0], it[1]) }
                .join(plotting_resolutions)
                .join(CONCAT.out.tsv)
                .join(chrom_pairs.groupTuple())
                .set { plotting_tasks }

            if (!params.plot_sig_interactions_cmap_lb) {
                plot_sig_interactions_cmap_lb = Math.min(params.nchg_log_ratio_cis,
                                                         params.nchg_log_ratio_trans)
            } else {
                plot_sig_interactions_cmap_lb = params.plot_sig_interactions_cmap_lb
            }

            plot_sig_interactions_cmap_ub = Math.max(plot_sig_interactions_cmap_lb,
                                                     params.plot_sig_interactions_cmap_ub)

            PLOT_SIGNIFICANT(
               plotting_tasks,
               plot_sig_interactions_cmap_lb,
               plot_sig_interactions_cmap_ub
            )
        }

    emit:
        tsv = CONCAT.out.tsv

}

process GENERATE_MASK {
    label 'duration_very_short'

    cpus 1

    input:
        path cytoband
        path gaps
        path other
        val zstd_compression_lvl

    output:
        path "mask.bed.zst", emit: bed

    shell:
        '''
        set -o pipefail

        # Implement logic to support optional gaps/cytoband input files
        touch empty.tmp

        if [[ '!{gaps}' == "" ]]; then
            ln -s empty.tmp gaps
        else
            ln -s '!{gaps}' gaps
        fi

        if [[ '!{cytoband}' == "" ]]; then
            ln -s empty.tmp cytoband
        else
            ln -s '!{cytoband}' cytoband
        fi

        if [[ '!{other}' == "" ]]; then
            ln -s empty.tmp other
        else
            ln -s '!{other}' other
        fi

        # Concatenate and sort intervals from gaps and cytoband files,
        # then merge overlapping intervals

        mkfifo gaps.bed.tmp cytoband.bed.tmp other.bed.tmp
        trap 'rm -f mkfifo gaps.bed.tmp cytoband.bed.tmp other.bed.tmp empty.tmp' EXIT

        # Decompress BED files and select BED3 columns
        zcat -f gaps | cut -f 2-4 >> gaps.bed.tmp &
        zcat -f cytoband | grep 'acen$' | cut -f 1-3 >> cytoband.bed.tmp &
        zcat -f other | cut -f 1-3 >> other.bed.tmp &

        # Concatenate, sort and merge intervals
        cat *.bed.tmp |
        sort -k1,1V -k2,2n --parallel=!{task.cpus} |
        bedtools merge -i stdin |
        zstd -T!{task.cpus} \\
             -!{zstd_compression_lvl} \\
             -o mask.bed.zst
        '''
}

process MASK_DOMAINS {
    label 'duration_very_short'
    tag "$sample"

    cpus 1

    input:
        tuple val(sample),
              path(domains)

        path mask
        val zstd_compression_lvl

    output:
        tuple val(sample),
              path(outname)

    shell:
        outname=domains.toString().replace(".bed.gz", "_filtered.bed.zst")
        '''
        set -o pipefail

        bedtools subtract -A \\
                          -a <(zcat -f '!{domains}') \\
                          -b <(zstdcat -f '!{mask}') |
            zstd -T!{task.cpus} \\
                 -!{zstd_compression_lvl} \\
                 -o '!{outname}'
        '''
}


process DUMP_CHROM_SIZES {
    label 'process_very_short'
    tag "$sample"

    input:
        tuple val(sample),
              path(hic),
              val(resolution)

    output:
        tuple val(sample),
              path("*.chrom.sizes"),
        emit: tsv

    shell:
        '''
        #!/usr/bin/env python3

        import hictkpy

        chroms = hictkpy.File("!{hic}", int("!{resolution}")).chromosomes()

        with open("!{sample}.chrom.sizes", "w") as f:
            for chrom, size in chroms.items():
                print(f"{chrom}\\t{size}", file=f)
        '''
}

process GENERATE_CHROMOSOME_PAIRS {
    label 'process_very_short'
    tag "$sample"

    input:
        tuple val(sample),
              path(hic),
              val(resolution),
              val(interaction_type)

    output:
        stdout emit: tsv

    shell:
        '''
        #!/usr/bin/env python3

        import hictkpy

        chroms = list(
            hictkpy.File("!{hic}", int("!{resolution}")).chromosomes().keys()
        )

        sample = "!{sample}"
        interaction_type = "!{interaction_type}"

        for i, chrom1 in enumerate(chroms):
            for chrom2 in chroms[i:]:
                do_print = interaction_type == "cis" and chrom1 == chrom2
                do_print |= interaction_type == "trans" and chrom1 != chrom2

                if do_print:
                    print(f"{sample}\\t{chrom1}\\t{chrom2}")
        '''
}

// TODO optimize: trans expected values can be computed in parallel
process EXPECTED {
    publishDir "${params.publish_dir}/nchg/",
        enabled: !!params.publish_dir,
        mode: params.publish_dir_mode
    tag "$sample"

    input:
        tuple val(sample),
              path(hic),
              val(resolution)

        val interaction_types
        val mad_max

    output:
        tuple val(sample),
              path(outname),
        emit: h5

    shell:
        opts=[
            "--resolution='${resolution}'",
            "--mad-max='${mad_max}'"
        ]

        suffix=""
        if (interaction_types.size() == 1) {
            if ("cis" in interaction_types) {
                opts.push("--cis-only")
                suffix=".cis"
            } else if ("trans" in interaction_types) {
                opts.push("--trans-only")
                suffix=".trans"
            }
        }

        outname="expected_values_${sample}${suffix}.h5"
        opts=opts.join(" ")
        '''
        NCHG expected \\
            '!{hic}' \\
            --output='!{outname}' \\
            !{opts}
        '''
}

process COMPUTE{
    tag "$sample ($chrom1:$chrom2)"


    input:
        tuple val(sample),
              val(chrom1),
              val(chrom2),
              path(hic),
              val(resolution),
              path(domains),
              path(expected_values)

        val bad_bin_fraction

    output:
        tuple val(sample),
              val(chrom1),
              val(chrom2),
              path("*.parquet", optional: true),
        emit: parquet

    shell:
        outname="${sample}.${chrom1}.${chrom2}.parquet"
        '''
        set -o pipefail

        zstdcat -f '!{domains}' > domains.bed

        NCHG compute \\
            '!{hic}' \\
            '!{outname}' \\
            --resolution='!{resolution}' \\
            --chrom1='!{chrom1}' \\
            --chrom2='!{chrom2}' \\
            --expected-values='!{expected_values}' \\
            --domains=domains.bed \\
            --bad-bin-fraction='!{bad_bin_fraction}'
        '''
}

process MERGE {
    tag "$sample ($interaction_type)"

    cpus 2

    input:
        tuple val(sample),
              path(parquets),
              val(interaction_type),
              path(chrom_sizes)

    output:
        tuple val(sample),
              val(interaction_type),
              path(outname),
        emit: parquet

    shell:
        input_prefix="${sample}"
        outname="${sample}.${interaction_type}.parquet"
        '''
        NCHG merge '!{input_prefix}' '!{outname}' \\
            --threads='!{task.cpus}'
        '''
}

process FILTER {
    tag "$sample ($interaction_type)"

    cpus 2

    input:
        tuple val(sample),
              val(interaction_type),
              path(parquet),
              val(fdr),
              val(log_ratio)

    output:
        tuple val(sample),
              val(interaction_type),
              path(outname),
        emit: parquet

    shell:
        outname="${sample}.${interaction_type}.filtered.parquet"
        '''
        NCHG filter \\
            '!{parquet}' \\
            '!{outname}' \\
            --fdr='!{fdr}' \\
            --log-ratio='!{log_ratio}' \\
            --threads='!{task.cpus}'
        '''
}

process VIEW {
    tag "$sample ($interaction_type)"

    input:
        tuple val(sample),
              val(interaction_type),
              path(parquet)

    output:
        tuple val(sample),
              val(interaction_type),
              path(outname),
        emit: tsv

    shell:
        outname="${sample}.${interaction_type}.filtered.tsv.gz"
        '''
        set -o pipefail

        NCHG view '!{parquet}' |
            pigz -9 -p !{task.cpus} > '!{outname}'
        '''
}

process CONCAT {
    publishDir "${params.publish_dir}/nchg/",
        enabled: !!params.publish_dir,
        mode: params.publish_dir_mode

    tag "$sample"

    input:
        tuple val(sample),
              val(interaction_types),
              path(tsvs)

    output:
        tuple val(sample),
              path(outname),
        emit: tsv

    shell:
        outname="${sample}.filtered.tsv.gz"
        tsvs_str=tsvs.join(" ")
        '''

        # Write the file header
        zcat '!{tsvs[0]}' |
            head -n 1 |
            pigz -9 > '!{outname}'

        for f in !{tsvs_str}; do
            # Skip file headers
            zcat "$f" | tail -n +2
        done |
            sort -k1,1V -k2,2n -k4,4V -k5,5n |
            pigz -9 -p !{task.cpus} >> '!{outname}'
        '''
}

process PLOT_EXPECTED {
    publishDir "${params.publish_dir}/plots/nchg/${sample}",
        enabled: !!params.publish_dir,
        mode: params.publish_dir_mode

    tag "$sample"

    input:
        tuple val(sample),
              path(h5)

        val interaction_types

    output:
        tuple val(sample),
              path("*.${params.plot_format}")

    shell:
        plot_cis="cis" in interaction_types
        plot_trans="trans" in interaction_types
        '''
        if [[ !{plot_cis} == true ]]; then
            plot_nchg_expected_values.py \\
                '!{h5}' \\
                '!{sample}_cis.!{params.plot_format}' \\
                --yscale-log \\
                --plot-cis
        fi

        if [[ !{plot_trans} == true ]]; then
            plot_nchg_expected_values.py \\
                '!{h5}' \\
                '!{sample}_trans.!{params.plot_format}' \\
                --plot-trans
        fi
        '''
}

process GET_HIC_PLOT_RESOLUTION {
    label 'process_very_short'
    tag "$sample"

    input:
        tuple val(sample),
              path(hic),
              val(resolution)

        val tgt_resolution

    output:
        stdout emit: tsv

    shell:
        '''
        #!/usr/bin/env python3

        import hictkpy

        best_res = int("!{resolution}")

        try:
            resolutions = hictkpy.MultiResFile("!{hic}").resolutions()

            tgt_res = int("!{tgt_resolution}")

            for res in resolutions:
                delta1 = abs(res - tgt_res)
                delta2 = abs(best_res - tgt_res)

                if delta1 < delta2:
                    best_res = res

        except RuntimeError:
            pass
        finally:
            print(f"!{sample}\\t{best_res}")
        '''
}

process PLOT_SIGNIFICANT {
    publishDir "${params.publish_dir}/plots/nchg/${sample}",
        enabled: !!params.publish_dir,
        mode: params.publish_dir_mode

    label 'process_very_high'
    tag "$sample"

    input:
        tuple val(sample),
              path(hic),
              val(resolution),
              path(parquet),
              val(chroms1),
              val(chroms2)

        val cmap_lb
        val cmap_ub

    output:
        tuple val(sample),
              path("*.${params.plot_format}")

    shell:
        chroms1_str=chroms1.join(" ")
        chroms2_str=chroms2.join(" ")
        '''
        chroms1=(!{chroms1_str})
        chroms2=(!{chroms2_str})

        commands=()
        for i in "${!chroms1[@]}"; do
            chrom1="${chroms1[$i]}"
            chrom2="${chroms2[$i]}"

            outname="!{sample}.$chrom1.$chrom2.!{params.plot_format}"

            command=(
                plot_significant_interactions.py \\
                    '!{hic}' \\
                    '!{parquet}' \\
                    "$chrom1" \\
                    "$chrom2" \\
                    "$outname" \\
                    --resolution='!{resolution}' \\
                    --min-value='!{cmap_lb}' \\
                    --max-value='!{cmap_ub}'
            )

            commands+=("${command[*]}")
        done

        printf '%s\\0' "${commands[@]}" |
            xargs -0 -I '{}' -P '!{task.cpus}' bash -c '{}'
        '''
}
