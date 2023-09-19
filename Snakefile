import os
import sys
from collections import defaultdict

smplsDict = {x['name']: x for x in config["samples"]}
smplsNames = smplsDict.keys()

snakefileDir = os.path.dirname(workflow.snakefile)
scriptDir = snakefileDir + "/scripts/"

###############################################################################
# input paths dict
###############################################################################
inputPaths = defaultdict(str)
for inKey in config["paths"]:
  inPath = os.path.normpath(config["paths"][inKey])

  inputPaths[inKey] = inPath


###############################################################################
# make all out directories if needed
###############################################################################
outDirs = defaultdict(str)
for outDirKey in config["out_dirs"]:
  outDir = os.path.normpath(config["out_dirs"][outDirKey])

  if not os.path.exists(outDir):
    os.mkdir(outDir)

  outDirs[outDirKey] = outDir


###############################################################################
# resource allocation
###############################################################################
availThreads = config["general"]["threads"]
availMem = config["general"]["memory"]


###############################################################################
# determine what are the final outputs and check if req'd files exist
###############################################################################
def final_output(modes):
  outs = []

  # file check
  if not modes["atac"] and not modes["adt"]:
    raise Exception("ATAC and adt modes are both turned off")

  if modes["atac"]:
    if not os.path.exists(inputPaths["atac_fastq_dir"]):
      raise Exception("ATAC fastq directory does not exist.")

    if not os.path.exists(inputPaths["cellranger_ref_dir"]):
      raise Exception("No cellranger-atac genome reference exists.")

    if not os.path.exists(inputPaths["cellranger_exe"]):
      raise Exception("No cellranger-atac executable exists.")

  if modes["adt"]:
    if not os.path.exists(inputPaths["adt_fastq_dir"]):
      raise Exception("ADT fastq directory does not exist.")

    if not os.path.exists(inputPaths["adt_catalog"]):
      raise Exception("No catalog of ADT tags/oligos exists.")

    if not os.path.exists(inputPaths["allowlist"]):
      raise Exception("No cell barcode allowlist exists.")

  # atac outs
  if modes["atac"] and modes["atac_level"] == "namesort":
    outs.append(
      expand(
        f"cr_out_{{smpl}}/outs/namesorted.bam",
        smpl = smplsNames))

    outs.append(
      expand(
        f"{outDirs['amulet']}/{{smpl}}/MultipletBarcodes_01.txt",
        smpl = smplsNames))

  elif modes["atac"] and modes["atac_level"] == "haystack":
    outs.append(
      expand(
        f"{outDirs['haystack']}/{{smpl}}/viralFrags.tsv",
        smpl = smplsNames))

  else:
    sys.exit("ATAC mode is not set correctly. Currently allowed values: 'namesort', 'haystack'")

  # adt outs
  if modes["adt"]:
      outs.append(
        expand(
          f"{outDirs['adt']}/{{smpl}}_output_count/output.mtx",
          smpl = smplsNames))

  return(outs)


rule all: 
  input:
    *final_output(config["general"]["run_modes"])


###############################################################################
# atac rules
###############################################################################
rule cellranger_count:
  output:
    touch(".snek_cr/{smpl}_cr_count.done")
  params:
    cellrangerAtac = config["paths"]["cellranger_exe"],
    fastqDir = config["paths"]["atac_fastq_dir"],
    ref = config["paths"]["cellranger_ref_dir"]
  threads: availThreads["cellranger_count"]
  resources:
    mem_mb = availMem["cellranger_count"]
  shell:
    """
      {params.cellrangerAtac} count \
        --id=cr_out_{wildcards.smpl} \
        --reference={params.ref} \
        --fastqs={params.fastqDir} \
        --sample={wildcards.smpl} \
        --localcores={threads}
    """

rule namesort:
  input:
    rules.cellranger_count.output
  output:
    "cr_out_{smpl}/outs/namesorted.bam"
  threads: availThreads["namesort"]
  resources:
    mem_mb = availMem["namesort"]
  shell:
    "samtools sort -@ {threads} -n -o {output} -O bam cr_out_{wildcards.smpl}/outs/possorted_bam.bam"


rule amulet:
  input:
    rules.cellranger_count.output
  output:
    "amulet_out/{smpl}/" + "MultipletBarcodes_01.txt"
  params:
    amulet_dir = config["paths"]["amulet_dir"],
    autosomes = config["paths"]["autosomes"],
    denylist = config["paths"]["denylist"],
    overlap = config["amulet_settings"]["overlap"]
  threads: availThreads["amulet"]
  resources:
    mem_mb = availMem["amulet"]
  shell:
    """
    {params.amulet_dir}/AMULET.sh \
      --expectedoverlap={params.overlap} \
      cr_out_{wildcards.smpl}/outs/fragments.tsv.gz \
      cr_out_{wildcards.smpl}/outs/singlecell.csv \
      {params.autosomes} \
      {params.denylist} \
      amulet_out/{wildcards.smpl} \
      {params.amulet_dir}
    """


###############################################################################
# adt rules
###############################################################################
# python script from https://github.com/pachterlab/kite/tree/master/featuremap
rule build_preIndex:
  input:
    inputPaths["adt_catalog"]
  output:
    t2g = outDirs["adt"] + "/adt.t2g",
    fa = outDirs["adt"] + "/adt.fa"
  params:
    feat_map = scriptDir + "featuremap.py"
  shell:
    "python {params.feat_map} {input} --t2g {output.t2g} --fa {output.fa}"


rule make_index:
  input:
    rules.build_preIndex.output.fa
  output:
    outDirs["adt"] + "/adt.idx"
  params:
    kmer = config["general"]["library"]["tag"][2] - config["general"]["library"]["tag"][1]
  shell:
    "kallisto index -i {output} -k {params.kmer} {input}"


def get_sample_fastqs(wildcards):
  inputs = [inputPaths["adt_fastq_dir"] + "/" + x for x in smplsDict[wildcards.smpl]["adt_fastqs"]]

  return(inputs)

def convert_list_to_bus_params(l):
  return(",".join([str(x) for x in l]))


rule make_bus:
  input:
    fastqs = get_sample_fastqs,
    idx = rules.make_index.output
  output:
    bus = outDirs["adt"] + "/{smpl}_output_bus/output.bus"
  params:
    cbc = convert_list_to_bus_params(config["general"]["library"]["cbc"]),
    umi = convert_list_to_bus_params(config["general"]["library"]["umi"]),
    tag = convert_list_to_bus_params(config["general"]["library"]["tag"]),
    out_dir = outDirs["adt"] + "/{smpl}_output_bus/",
  threads: availThreads["make_bus"]
  resources:
    mem_mb = availMem["make_bus"]
  shell:
    "kallisto bus -x {params.cbc}:{params.umi}:{params.tag} -i {input.idx} -t {threads} -o {params.out_dir} {input.fastqs}"


rule correct_cbc:
  input:
    rules.make_bus.output.bus
  output:
    rules.make_bus.params.out_dir + "output.corrected.bus"
  params:
    allowlist = inputPaths["allowlist"]
  shell:
    "bustools correct -w {params.allowlist} -o {output} {input} "


rule sort_bus:
  input:
    rules.correct_cbc.output
  output:
    rules.make_bus.params.out_dir + "output.corrected.sorted.bus"
  threads: availThreads["sort_bus"]
  shell:
    "bustools sort -t {threads} -o {output} {input}"


rule make_count_matrix:
  input:
    rules.sort_bus.output
  output:
    outDirs["adt"] + "/{smpl}_output_count/output.mtx",
  params:
    count_output = outDirs["adt"] + "/{smpl}_output_count/",
    t2g = rules.build_preIndex.output.t2g,
    ec = rules.make_bus.params.out_dir + "matrix.ec",
    txt = rules.make_bus.params.out_dir + "transcripts.txt"
  shell:
    "bustools count -o {params.count_output} --genecounts -g {params.t2g} -e {params.ec} -t {params.txt} {input}"
