# viral-asapseq-snakemake
Upstream processing for single cell ASAPseq for viral vector detection. This pipeline is based on a previous pipeline that was designed for processing of sequencing techniques with antibody derived tags (see [scc-proc](github.com/betts-lab/scc-proc))

## Necessary files
### 1) `config.yaml` file
Example configuration is shown below.

- If ATAC mode is enabled, the sample name MUST match with the sample name used to the demultiplex the fastq files during cellranger-atac mkfastq.
- Multiple samples can be added to the config assuming they follow the same antibody mixture and bead setup.
- All fastq files must be in the fastq_dir path
  - Specific filenames are provided in the `samples` attribute of the yaml file that correspond to the ADT files
  - Order of fastq files in the `samples` attribute corresponds to the `cbc`, `umi`, `tag` numbering system (refer to kallisto documentation).
- In the general settings:
  - `run_modes` section will toggle different modalities to be run.
  - `threads` and `memory` allow rule-specific thread/memory setting respectively. This pipeline is designed for LSF based submission but the lsf.yml file can be edited as needed (see below)

```
paths:
  atac_fastq_dir: ""
  cellranger_ref_dir: ""
  cellranger_exe: "~/pkg/cellranger-atac-2.1.0/bin/cellranger-atac"

  adt_fastq_dir: ""
  adt_catalog: ""
  allowlist: ""

out_dirs:
  atac: "." # due to cellranger settings, this cannot be currently changed.
  haystack: "haystack_out"
  adt: "adt_out"

general:
  run_modes:
    atac: True # boolean
    atac_level: "namesort" # must be namesort or haystack
    adt: False # boolean

  library:
    cbc: [1,0,16]
    umi: [0,0,10]
    tag: [2,0,15]

  threads:
    cellranger_count: 12
    namesort: 8
    make_bus: 4
    sort_bus: 8

  memory:
    cellranger_count: 120000
    namesort: 96000
    make_bus: 36000

samples:
  - name: "A"
    adt_fastqs: ["A_S2_R1_001.fastq.gz", "A_S2_R2_001.fastq.gz", "A_S2_R3_001.fastq.gz"]
```

### 2) `lsf.yaml`
This file allows for specific LSF submission settings based on rules.
```
__default__:
  - "-q normal"


cellranger_count:
  - "-q denovo"

namesort:
  - "-q denovo"
```

### 3) fastq files for ADT; directory of fastq files for ATAC
Refer to config above

### 4) Antibody catalog csv file (required if doing ADT portion)
Comma separted file with barcodes and descriptions. Example below:

```
CD40,CTCAGATGGAGTATG
CD44,AATCCTTCCGAATGT
CD48,CTACGACGTAGAAGA
CD21,AACCTAGTAGTTCGG
```

### 4) Allowlist cell barcode file (required if doing ADT portion)
Allowlist file containing accepted cell barcodes. One per line. Example below:

```
GTCTGCTATGTCTA
GATGATGCATAGAA
```

## Setup
1. Create conda env and set up LSF job submission for Snakemake (follow instructions [here](https://github.com/Snakemake-Profiles/lsf)).
```
conda create --name <env> --file env.yml
```

2. Run Snakemake
```
bsub -e snek.e -o snek.o snakemake --configfile=config.yaml --profile=lsf -s <path-to-Snakefile>
```

The overall process will look like this (to be added):

<!-- ![DAG image](dag.png) -->

3. Downstream analysis of your own choosing. Output files are as follows (to be added)...

## Credits
- [Conda](https://docs.conda.io/en/latest/)
- [Snakemake](https://snakemake.readthedocs.io/en/stable/)
- [Kallisto](https://pachterlab.github.io/kallisto/)
- [Kallisto KITE featureMap.py](https://github.com/pachterlab/kite/tree/master/featuremap)
- [Bustools](https://github.com/BUStools/bustools)
- [kallisto | bustools KITE protocol](https://bustools.github.io/BUS_notebooks_R/10xv3.html)
    - I made my Snakemake pipeline heavily based on their pipeline.
