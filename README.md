<!--
Copyright (C) 2022 Roberto Rossini <roberros@uio.no>

SPDX-License-Identifier: MIT
-->

# Nextflow workflow to call TAD cliques

[![CI](https://github.com/robomics/call_tad_cliques/actions/workflows/ci.yml/badge.svg)](https://github.com/robomics/call_tad_cliques/actions/workflows/ci.yml)

This repository hosts a Nextflow workflow to call [TAD cliques](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-021-07815-8).

The workflow is largely based on [this](https://github.com/Chrom3D/INC-tutorial) tutorial.

## Requirements

### Software requirements

- Nextflow (at least version: v20.07.1. Pipeline was developed using v22.10.8)
- Docker or Singularity/Apptainer

Running the pipeline without containers is technically possible, but it is not recommended.

<details>
<summary>If you absolutely cannot use containers...</summary>

Have a look at the `env.yml` for the list of dependencies to be installed.

To install the dependencies in a Conda environment named `myenv`, run the following:

```bash
conda env update --name myenv --file env.yml --prune
```

You will also need to compile `NCHG` from the source code available at [paulsengroup/NCHG](https://github.com/paulsengroup/NCHG).

</details>

### Required input files

The workflow can be run in two ways:
1. Using a sample sheet (recommended, supports processing multiple samples at once)
2. By specifying options directly on the CLI or using a config

#### Using a samplesheet

The samplesheet should be a TSV file with the following columns:

| sample       | hic_file                               | resolution | tads                      |
|--------------|----------------------------------------|------------|---------------------------|
| sample_name  | myfile.hic                             | 50000      | tads.bed                  |
| 4DNFI74YHN5W | 4DNFI74YHN5W.mcool::/resolutions/50000 | 50000      | 4DNFI74YHN5W_domains.bed  |


- __sample__: Sample names/ids. This field will be used as prefix to in the output file names (see [below](#running-the-workflow)).
- __hic_file__: Path to a file in .hic or Cooler format.
- __resolution__: Resolution to be used for the data analysis (50-100kbp are good starting points).
- __tads__ (optional) : path to a BED3+ file with the list of TADs. When not specified, the workflow will use [hicFindTADs](https://hicexplorer.readthedocs.io/en/latest/content/tools/hicFindTADs.html) from the [HiCExplorer](https://github.com/deeptools/HiCExplorer) suite to call TADs.

URI syntax for multi-resolution Cooler files is supported (e.g. `myfile.mcool::/resolutions/bin_size`).

Furthermore, all contact matrices (and TADs when provided) should use the same reference genome assembly.

<details>
<summary> <b>Without using a samplesheet</b> </summary>

To run the workflow without a samplesheet is not available, the following parameters are required:

- __sample__
- __hic_file__
- __tads__

Parameters have the same meaning as the header fields outlined in the [previous section](#using-a-samplesheet).

The above parameters can be passed directly through the CLI when calling `nextflow run`:

```bash
nextflow run --sample='4DNFI74YHN5W' \
             --hic_file='data/4DNFI74YHN5W.mcool' \
             --resolution=50000
             ...
```

Alternatively, parameters can be written to a `config` file:
```console
user@dev:/tmp$ cat myconfig.txt

sample       = '4DNFI74YHN5W'
hic_file     = 'data/4DNFI74YHN5W.mcool'
resolution   = 50000
```

and the `config` file is then passed to `nextflow run`:
``` bash
nextflow run -c myconfig.txt ...
```

</details>

### Optional files and parameters

In addition to the mandatory parameters, the pipeline accepts the following parameters:

- __cytoband__: path to a [cytoband](https://software.broadinstitute.org/software/igv/cytoband) file. Used to mask centromeric regions.
- __assembly_gaps__: path to a BED file with the list of assembly gaps/unmappable regions.
- __custom_mask__: path to a BED file with a list of custom regions to be masked out.

Note that NCHG by default uses the `MAD-max` filter to remove bins with suspiciously high or low marginals, so providing the above files is usually not requirerd.
One exception is when dealing with genomes affected by structural variants, in which case we reccommend masking out these regions using __custom_mask__.

- __hicexplorer_hic_norm__: normalization to use when calling TADs from .hic files.
- __hicexplorer_cool_norm__: normalization to use when calling TADs from .\[m\]cool files.
- __nchg_mad_max__: cutoff used by NCHG when performing the `MAD-max` filtering.
- __nchg_bad_bin_fraction__: bad bin fraction used by NCHG to discard domains overlapping with a high fraction of bad bins.
- __nchg_fdr_cis__: adjusted pvalue used by NCHG to filter significant cis interactions.
- __nchg_log_ratio_cis__: log ratio used by NCHG to filter significant cis interactions.
- __nchg_fdr_trans__: adjusted pvalue used by NCHG to filter significant trans interactions.
- __nchg_log_ratio_trans__: log ratio used by NCHG to filter significant trans interactions.
- __clique_size_thresh__: minimum clique size.
- __call_cis_cliques__: call cliques overlapping cis regions of the Hi-C matrix.
- __call_trans_cliques__: call cliques overlapping trans regions of the Hi-C matrix.

By default, the workflow results are published under `result/`. The output folder can be customized through the __outdir__ parameter.

For a complete list of parameters supported by the workflow refer to the workflow main [config](nextflow.config) file.

## Running the workflow

First, download the example datasets using script `utils/download_example_datasets.sh`.

```bash
# This will download files inside folder data/
utils/download_example_datasets.sh data/
```

Next, create a `samplesheet.tsv` file like the follwing (make sure you are using tabs, not spaces!)

```tsv
sample   hic_file      resolution    tads
example  data/4DNFI74YHN5W.mcool   50000   
```

Finally, run the workflow with:
```console
user@dev:/tmp$ nextflow run --max_cpus=8 \
                            --max_memory=16.GB \
                            --max_time=2.h \
                            --sample_sheet=samplesheet.tsv \
                            --outdir=data/results/ \
                            https://github.com/robomics/call_tad_cliques \
                            -r v0.4.0 \
                            -with-singularity  # Replace this with -with-docker to use Docker instead

 N E X T F L O W   ~  version 24.04.2

Launching `./main.nf` [irreverent_wescoff] DSL2 - revision: 0be6cbedc5

executor >  local (41)
[7b/fb313f] SAMPLESHEET:CHECK_SYNTAX                                   | 1 of 1 ✔
[b0/4f77c1] SAMPLESHEET:CHECK_FILES                                    | 1 of 1 ✔
[0b/0ba9ff] TADS:SELECT_NORMALIZATION_METHOD (example)                 | 1 of 1 ✔
[09/ca2567] TADS:APPLY_NORMALIZATION (example (50000; weight))         | 1 of 1 ✔
[da/a14762] TADS:HICEXPLORER_FIND_TADS (example)                       | 1 of 1 ✔
[-        ] TADS:COPY                                                  -
[92/f261fc] NCHG:GENERATE_MASK                                         | 1 of 1 ✔
[8e/40c02a] NCHG:MASK_DOMAINS (example)                                | 1 of 1 ✔
[36/f950a6] NCHG:EXPECTED (example)                                    | 1 of 1 ✔
[bf/6fbee1] NCHG:GENERATE_CHROMOSOME_PAIRS (example)                   | 1 of 1 ✔
[47/945cac] NCHG:DUMP_CHROM_SIZES (example)                            | 1 of 1 ✔
[16/bb5e75] NCHG:COMPUTE (example (chr1:chr1))                         | 21 of 21 ✔
[51/55dfbb] NCHG:MERGE (example (cis))                                 | 1 of 1 ✔
[f1/196069] NCHG:FILTER (example (cis))                                | 1 of 1 ✔
[01/c17c93] NCHG:VIEW (example (cis))                                  | 1 of 1 ✔
[a7/68bab8] NCHG:CONCAT (example)                                      | 1 of 1 ✔
[bd/e70a4d] NCHG:PLOT_EXPECTED (example)                               | 1 of 1 ✔
[a4/235495] NCHG:GET_HIC_PLOT_RESOLUTION (example)                     | 1 of 1 ✔
[78/19ba93] NCHG:PLOT_SIGNIFICANT (example)                            | 1 of 1 ✔
[76/b21803] CLIQUES:CALL (example)                                     | 1 of 1 ✔
[e7/812637] CLIQUES:PLOT_MAXIMAL_CLIQUE_SIZE_DISTRIBUTION_BY_TAD (cis) | 1 of 1 ✔
[b6/d1f15a] CLIQUES:PLOT_CLIQUE_SIZE_DISTRIBUTION (cis)                | 1 of 1 ✔
Completed at: 07-Jun-2024 16:10:27
Duration    : 1m 28s
CPU hours   : 0.1
Succeeded   : 41
```

This will create a `data/results/` folder with the following files:
- `cliques/example_cis_cliques.tsv.gz` - TSV with the list of cliques computed from cis significant interactions (i.e. both intra and inter-chromosomal interactions).
- `cliques/example_cis_domains.bed.gz` - BED file with the list of domains part of cliques computed from cis significant interactions. The last column encodes the domain ID.
- `nchg/example.filtered.tsv.gz` - TSV with the statistically significant interactions detected by NCHG.
- `nchg/expected_values_example.cis.h5` - HDF5 file with the expected values computed by NCHG.
- `plots/cliques/cis_clique_size_distribution*` - Plots showing the clique size distribution.
- `plots/cliques/cis_tad_max_clique_size_distribution*` - Plots showing the maximal clique size distribution.
- `plots/nchg/example/example.*.*.png` - Plots showing the log ratio computed by NCHG for each chromosome pair analyzed.
- `plots/nchg/example/example_cis.png` - Plot showing the expected value profile computed by NCHG.
- `tads/example_tads.bed.gz` - TADs used to generate the list of genomic coordinates to be tested for significance.

The list of pairs of interacting domains can be generated using `bin/generate_cliques_bedpe.py`

```console
user@dev:/tmp$ bin/generate_cliques_bedpe.py data/results/example_cis_domains.bed.gz data/results/example_cis_cliques.tsv.gz |
chr3	94300000	95850000	chr3	94300000	95850000	CLIQUE_#0
chr3	94300000	95850000	chr3	108150000	108950000	CLIQUE_#0
chr3	94300000	95850000	chr3	116000000	116900000	CLIQUE_#0
chr3	94300000	95850000	chr3	137750000	138500000	CLIQUE_#0
chr3	94300000	95850000	chr3	152100000	153900000	CLIQUE_#0
chr3	108150000	108950000	chr3	94300000	95850000	CLIQUE_#0
chr3	108150000	108950000	chr3	108150000	108950000	CLIQUE_#0
chr3	108150000	108950000	chr3	116000000	116900000	CLIQUE_#0
chr3	108150000	108950000	chr3	137750000	138500000	CLIQUE_#0
chr3	108150000	108950000	chr3	152100000	153900000	CLIQUE_#0
```

<details>
<summary>Troubleshooting</summary>

If you get permission errors when using `-with-docker`:
- Pass option `-process.containerOptions="--user root"` to `nextflow run`

If you get an error similar to:
```
Cannot find revision `v0.4.0` -- Make sure that it exists in the remote repository `https://github.com/robomics/call_tad_cliques`
```

try to remove folder `~/.nextflow/assets/robomics/call_tad_cliques` before running the workflow

</details>

## Getting help

If you are having trouble running the workflow feel free to reach out by starting a new discussion [here](https://github.com/robomics/call_tad_cliques/discussions).

Bug reports and feature requests can be submitted by opening an [issue](https://github.com/robomics/call_tad_cliques/issues).
