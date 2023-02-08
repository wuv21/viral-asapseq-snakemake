import yaml
import os
import re
import sys
from collections import defaultdict

smplsNames = config["samples"].keys()
snakefileDir = os.path.dirname(workflow.snakefile)
scriptDir = snakeFileDir + "scripts/"

###############################################################################
# input paths dict
###############################################################################
inputPaths = defaultdict(str)
for inKey in config["paths"]:
  inPath = os.path.normpath(config["paths"][inKey])

  inputPaths[inKey] = inPath


###############################################################################
# make all out directories
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

  # TODO run checks....

  if not modes["atac"] and not modes["adt"]:
    sys.exit("ATAC and adt modes are both turned off")

  # atac outs
  if modes["atac"] and modes["atac_level"] == "namesort":
    outs.append(
      expand(
        f"{outDirs['atac']}/{{smpl}}/outs/namesorted.bam",
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
     outDirs["atac"] + "/{smpl}/outs/possorted_bam.bam"
  params:
    outDir = outDirs["atac"],
    cellrangerAtac = config["paths"]["cellranger_exe"],
    fastqDir = config["paths"]["rna_fastq_dir"],
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
    outDirs["atac"] + "/{smpl}/outs/namesorted.bam"
  threads: availThreads["namesort"]
  resources:
    mem_mb = availMem["namesort"]
  shell:
    "samtools sort -@ {threads} -n -o {output} -O bam {input}"
  

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
    feat_map = scriptDir + "/featuremap.py"
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
  inputs = [inPath["adt_fastq_dir"] + "/" + x for x in config["samples"]["adt_fastqs"][wildcards.smpl]]

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
    cbc = convert_list_to_bus_params(config["general_settings"]["cbc"]),
    umi = convert_list_to_bus_params(config["general_settings"]["umi"]),
    tag = convert_list_to_bus_params(config["general_settings"]["tag"]),
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
    allowlist = inPath["allowlist"]
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