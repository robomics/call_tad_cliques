// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

def strip_resolution_from_cooler_uri(uri) {
    file(uri.replaceFirst(/::\/resolutions\/\d+$/, ""), checkIfExists: true)
}

// Workaround for optional input files: https://github.com/nextflow-io/nextflow/issues/1694
def make_optional_input(path) {
    if (path?.trim()) {
        return [file(path)]
    }
    return []
}

def parse_sample_sheet_row(row) {
    hic_fname = strip_resolution_from_cooler_uri(row.hic_file)
    files = [hic_fname]

    tads = make_optional_input(row.tads)
    mask = make_optional_input(row.mask)

    tuple(row.sample,
          files,
          tads,
          mask)
}


workflow SAMPLESHEET {

    take:
        sample_sheet
        sample
        hic_file
        resolution
        tads
        mask

    main:
        if (sample_sheet) {
            if(sample) log.warn("'sample' parameter is ignored when 'sample_sheet' is defined")
            if(hic_file) log.warn("'hic_file' parameter is ignored when 'sample_sheet' is defined")
            if(mask) log.warn("'mask' parameter is ignored when 'sample_sheet' is defined")
            if(resolution) log.warn("'resolution' parameter is ignored when 'sample_sheet' is defined")

            CHECK_SYNTAX(file(sample_sheet, checkIfExists: true))
        } else {
            if(!sample || !hic_file || !resolution) {
                log.error("The following parameters are required when 'sample_sheet' is undefined! Required parameters: sample, hic_file, resolution")
            }

            GENERATE(
                sample,
                hic_file,
                resolution,
                tads ? tads : "",
                mask ? mask : ""
            )

            CHECK_SYNTAX(GENERATE.out.tsv)
        }

        CHECK_SYNTAX.out.tsv.set { sample_sheet }

        sample_sheet
            .splitCsv(sep: "\t", header: true)
            .map {
                    it = parse_sample_sheet_row(it)
                    // Concatenate path to coolers and optional files
                    it[1] + it[2] + it[3]
            }
            .flatten()
            .unique()
            .filter { !!it }
            .set { files_from_sample_sheet }

        CHECK_FILES(
            sample_sheet,
            files_from_sample_sheet.collect()
        )

    emit:
        tsv = CHECK_FILES.out.tsv

}

process GENERATE {
    label 'duration_very_short'

    cpus 1

    input:
        val sample
        val hic_file
        val resolution
        val tads
        val mask

    output:
        path "sample_sheet.tsv", emit: tsv

    shell:
        '''
        for param in '!{sample}' '!{hic_file}' '!{resolution}'; do
            if [[ "$param" == 'null' ]]; then
                2>&1 echo 'Parameters sample, hic_file is required when no samplesheet is provided!'
                2>&1 echo "sample='!{sample}'; hic_file='!{hic_file}'; resolution='!{resolution}'"
                exit 1
            fi
        done

        printf 'sample\\thic_file\\tresolution\\ttads\\tmask\\n' > sample_sheet.tsv
        printf '%s\\t%s\\t%s\\t%s\\t%s\\n' '!{sample}' \\
                                 '!{hic_file}' \\
                                 '!{resolution}' \\
                                 '!{tads}' \\
                                 '!{mask}' >> sample_sheet.tsv
        '''
}

process CHECK_SYNTAX {
    label 'duration_very_short'

    cpus 1

    input:
        path sample_sheet

    output:
        path "${sample_sheet}", includeInputs: true, emit: tsv

    shell:
        '''
        parse_samplesheet.py --detached '!{sample_sheet}' > /dev/null
        '''
}

process CHECK_FILES {
    label 'duration_very_short'

    cpus 1

    input:
        path sample_sheet
        path files

    output:
        path "*.ok", emit: tsv

    shell:
        '''
        parse_samplesheet.py '!{sample_sheet}' > '!{sample_sheet}.ok'
        '''
}
