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

- Nextflow (at least version: v20.07.1. Pipeline was developed using v22.10.1)
- Docker or Singularity/Apptainer

Running the pipeline without containers is technically possible, but it is not recommended.

<details>
<summary>If you absolutely cannot use containers...</summary>

Have a look at the `env.yml` for the list of dependencies to be installed.

To install the dependencies in a Conda environment named `myenv`, run the following:

```bash
conda env update --name myenv --file env.yml --prune 
```

You will also need to compile `NCHG` from the source code available at [Chrom3D/preprocess_scripts](https://github.com/Chrom3D/preprocess_scripts).

Check out the `Dockerfile` from this repo for an example of how this can be done using Conda.

</details>

### Required input files

The workflow can be run in two ways:
1. Using a sample sheet (recommended, supports processing multiple samples at once)
2. By specifying options directly on the CLI or using a config

#### Using a samplesheet

The samplesheet should be a TSV file with the following columns:

| sample       | cooler_cis                                | cooler_trans                                | tads                     |
|--------------|-------------------------------------------|---------------------------------------------|--------------------------|
| sample_name  | high_resolution.cool                      | low_resolution.cool                         | tads.bed                 |
| 4DNFI74YHN5W | 4DNFI74YHN5W.mcool::/resolutions/50000    | 4DNFI74YHN5W.mcool::/resolutions/1000000    | 4DNFI74YHN5W_domains.bed |

- __sample__: Sample names/ids. This field will be used as prefix to in the output file names (see [below](#running-the-workflow)).
- __cooler_cis__: path to a cooler file with the contact matrix used to read intra-chromosomal/cis interactions (usually a matrix with resolution ~50kbp).
- __cooler_trans__: path to a cooler file with the contact matrix used to read intra-chromosomal/trans interactions (usually a matrix with resolution ~1Mbp).
- __tads__ (optional) : path to a BED3+ file with the list of TADs. When not specified, the workflow will use [hicFindTADs](https://hicexplorer.readthedocs.io/en/latest/content/tools/hicFindTADs.html) from the [HiCExplorer](https://github.com/deeptools/HiCExplorer) suite to call TADs.

URI syntax for multi-resolution Cooler files is supported (e.g. `myfile.mcool::/resolutions/bin_size`).

Furthermore, all contact matrices (and TADs when provided) should use the same reference genome assembly.

<details>
<summary> <b>Without using a samplesheet</b> </summary>

To run the workflow without a samplesheet is not available, the following parameters are required:

- __sample__
- __cooler_cis__
- __cooler_trans__

Parameters have the same meaning as the header fields outlined in the [previous section](#using-a-samplesheet).

The above parameters can be passed directly through the CLI when calling `nextflow run`:

```bash
nextflow run --sample='4DNFI74YHN5W' \
             --cooler_cis='data/4DNFI74YHN5W.mcool::/resolutions/100000' \
             --cooler_trans='data/4DNFI74YHN5W.mcool::/resolutions/1000000' \
             ...
```

Alternatively, parameters can be written to a `config` file:
```console
user@dev:/tmp$ cat myconfig.txt

sample       = '4DNFI74YHN5W'
cooler_cis   = 'data/4DNFI74YHN5W.mcool::/resolutions/100000'
cooler_trans = 'data/4DNFI74YHN5W.mcool::/resolutions/1000000'
```

and the `config` file is then passed to `nextflow run`:
``` bash
nextflow run -c myconfig.txt ...
```

</details>

### Optional files and parameters

In addition to the mandatory parameters, providing the following parameters is highly recommended:

- __cytoband__: path to a [cytoband](https://software.broadinstitute.org/software/igv/cytoband) file. Used to mask centromeric regions.
- __assembly_gaps__: path to a BED file with the list of assembly gaps/unmappable regions.

UCSC publishes the above files for common reference genomes: [link](https://hgdownload.cse.ucsc.edu/goldenPath/) (files are usually named `cytoBand.txt.gz` and `gaps.txt.gz` respectively).

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
sample   cooler_cis      cooler_trans    tads
example  data/4DNFI74YHN5W.mcool::/resolutions/50000   data/4DNFI74YHN5W.mcool::/resolutions/1000000
```

Finally, run the workflow with:
```console
user@dev:/tmp$ nextflow run --sample-sheet=samplesheet.tsv \
                            --cytoband=data/cytoBand.txt.gz \
                            --assembly_gaps=data/gaps.txt.gz \
                            --outdir=data/results/ \
                            https://github.com/robomics/call_tad_cliques \
                            -r v0.1.0 \
                            -with-singularity  # Replace this with -with-docker to use Docker instead

N E X T F L O W  ~  version 22.10.1
Launching `https://github.com/robomics/call_tad_cliques` [focused_bohr] DSL2 - revision: e3eb8d8633 [v0.0.10]
executor >  local (14)
[1a/d8a231] process > check_sample_sheet                             [100%] 1 of 1 ✔
[8d/6f96cc] process > process_sample_sheet                           [100%] 1 of 1 ✔
[0a/bb2a08] process > extract_chrom_sizes_from_cooler                [100%] 1 of 1 ✔
[c0/8ef7a2] process > generate_bed_mask                              [100%] 1 of 1 ✔
[46/68ecdb] process > process_tads (1)                               [100%] 1 of 1 ✔
[a7/538a68] process > fill_gaps_between_tads (1)                     [100%] 1 of 1 ✔
[b3/b678b9] process > bedtools_bed_setdiff (1)                       [100%] 1 of 1 ✔
[fc/6de0ad] process > map_intrachrom_interactions (1)                [100%] 1 of 1 ✔
[09/734100] process > select_significant_intrachrom_interactions (1) [100%] 1 of 1 ✔
[1e/cd7d31] process > collect_interchrom_interactions (1)            [100%] 1 of 1 ✔
[f5/c2feb2] process > select_significant_interchrom_interactions (1) [100%] 1 of 1 ✔
[be/c4afd0] process > map_interchrom_interactions_to_tads (1)        [100%] 1 of 1 ✔
[89/27c01f] process > merge_interactions (1)                         [100%] 1 of 1 ✔
[cc/399a12] process > call_cliques (1)                               [100%] 1 of 1 ✔
Completed at: 13-Jan-2023 20:59:06
Duration    : 6m 53s
CPU hours   : 0.1
Succeeded   : 14
```

This will create a `data/results/` folder with the following files:

- `example_clique_interactions.bedpe` - BEDPE with the list of interacting TADs part of cliques.
- `example_clique_sizes.bedGraph` - bedGraph file with the list of TADs with the size of the maximal clique to which they belong.


<details>
<summary> Sample output files </summary>

```console
user@dev:/tmp$ head data/results/*
==> data/results/example_clique_interactions.bedpe <==
chr1	33600000	34500000	chr1	33600000	34500000
chr1	33600000	34500000	chr1	63100000	64550000
chr1	33600000	34500000	chr1	73850000	75600000
chr1	33600000	34500000	chr1	87150000	88300000
chr1	33600000	34500000	chr1	92950000	94100000
chr1	33600000	34500000	chr1	131900000	132500000
chr1	33600000	34500000	chr1	133000000	133650000
chr1	63100000	64550000	chr1	63100000	64550000
chr1	63100000	64550000	chr1	73850000	75600000
chr1	63100000	64550000	chr1	87150000	88300000

==> data/results/example_clique_sizes.bedGraph <==
chr1	33600000	34500000	5
chr1	63100000	64550000	5
chr1	73850000	75600000	5
chr1	87150000	88300000	5
chr1	92950000	94100000	5
chr1	131900000	132500000	5
chr1	133000000	133650000	5
chr10	61700000	62850000	5
chr10	69350000	70300000	5
chr10	76250000	77500000	5
```

</details>

<details>
<summary>Troubleshooting</summary>

If you get permission errors when using `-with-docker`:
- Pass option `-process.containerOptions="--user root"` to `nextflow run`

If you get an error similar to:
```
Cannot find revision `v0.1.0` -- Make sure that it exists in the remote repository `https://github.com/robomics/call_tad_cliques`
```

try to remove folder `~/.nextflow/assets/robomics/call_tad_cliques` before running the workflow

</details>

## Getting help

If you are having trouble running the workflow feel free to reach out by starting a new discussion [here](https://github.com/robomics/call_tad_cliques/discussions).

Bug reports and feature requests can be submitted by opening an [issue](https://github.com/robomics/call_tad_cliques/issues).
