<!--
Copyright (C) 2022 Roberto Rossini <roberros@uio.no>

SPDX-License-Identifier: MIT
-->

# Nextflow workflow to call TAD cliques

[![CI](https://github.com/robomics/call_tad_cliques/actions/workflows/ci.yml/badge.svg)](https://github.com/robomics/call_tad_cliques/actions/workflows/ci.yml)

This repository hosts a Nextflow workflow to call TAD cliques.

The workflow is largely based on [this](https://github.com/Chrom3D/INC-tutorial) tutorial.

## Requirements

### Software requirements

- Nextflow (tested with v22.10.1, DSL2 support required)
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

The workflow requires the following input files:

- `--mcools` - one or more Hi-C matrices in .mcool format
- `--outdir` - directory where to store output files
- `--cytoband` - CytoBand file. Required to mask centromeric regions
- `--assembly_gaps` - Assembly gaps in BED3+ format

The complete list of supported parameters can be found in the `nextflow.config` file.

#### Additional requirements

- All .mcool files should use the same reference genome
- All .mcool files should have matrices for (at least) the resolutions specified by the `--cis_bin_size` and `--trans_bin_size` options (50kbp and 1Mbp by default)

<details>
<summary>Notes</summary>

The `--tads` option can be used to pass one or more BED3+ files with the list of TADs.

When this option is not passed, the workflow will call TADs
using [HiC-Explorer](https://hicexplorer.readthedocs.io/en/latest/content/tools/hicFindTADs.html#hicfindtads).

When provided, BED files should be in the same order as the .cool files, so that the first BED file is mapped to the first .cool etc...

Example:

| Hi-C matrix                        | BED                              |
|------------------------------------|----------------------------------|
| condition_001.mcool                | condition_001.bed                |
| mysamplehaveverycomplexnames.mcool | mysamplehaveverycomplexnames.bed |

</details>

## Running the workflow

Throughout this section we assume input files are stored inside a `data/` folder.

```bash
nextflow run --mcools='data/sample_001.mcool' \
             --cytoband=data/cytoBand.txt.gz \
             --assembly_gaps=data/gaps.txt.gz \
             --outdir=results \
             https://github.com/robomics/call_tad_cliques \
             -r v0.0.3 \
             -with-singularity  # Replace this with -with-docker to use Docker instead
```

This will create a `result/` folder with the following files:

- `sample_001_clique_interactions.bedpe`
- `sample_001_clique_sizes.bed`
- `sample_001_clique_stats.tsv`
- `sample_001_clique_tad_interactions.tsv`

<!-- TODO: describe output files -->

To test the workflow using publicly available mouse datasets, run the following:

```bash
nextflow run --mcools='https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/8b9d6836-0d6b-4c2f-9eaa-323f4fd7b6e4/4DNFI74YHN5W.mcool' \
             --cytoband='https://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/cytoBand.txt.gz' \
             --assembly_gaps='https://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/gap.txt.gz' \
             --outdir=results/ \
             https://github.com/robomics/call_tad_cliques \
             -r v0.0.3 \
             -with-singularity
```

This should take approximately 5-10 minutes (on slow internet connections it could take significantly longer).

<details>
<summary>Troubleshooting</summary>


If you get spurious errors about missing files, try one of the following:
- Re-run the workflow.
- Manually download files and pass the local file paths directly to the workflow.

If you get permission errors when using `-with-docker`:
- Pass option `-process.containerOptions="--user root"` to `nextflow run`

</details>
