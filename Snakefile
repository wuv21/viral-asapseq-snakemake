import yaml
import os
import re
import sys
from collections import defaultdict
from math import floor

fastqDir = os.path.normpath(config["paths"]["fastq_dir"])
smplsNames = config["samples"]

# smplDict = {}
# for s in config["samples"]:
#   smplDict[s["name"]] = s

########################################
# make all out directories
########################################
outDirs = defaultdict(str)
for outDirKey in config["out_dirs"]:
  outDir = os.path.normpath(config["out_dirs"][outDirKey])

  if not os.path.exists(outDir):
    os.mkdir(outDir)

  outDirs[outDirKey] = outDir

########################################
# basic checking to make sure required directories exist and correct
########################################





########################################
# resource allocation
########################################
threads = config["general"]["threads"]
availMem = config["general"]["memory"]

if threads < len(config["samples"]):
  nThreadsPerSample = threads

else:
  nThreadsPerSample = floor(threads / len(config["samples"]))


if nThreadsPerSample == threads:
  nMemPerSample = availMem

else:
  nMemPerSample = floor(availMem / len(config["samples"]))


########################################
# rules
########################################
def final_output(mode):
  print(outDirs)
  print(smplsNames)

  if mode == "namesort":
    return expand(
      f"{outDirs['atac']}/{{smpl}}/outs/namesorted.bam",
      outDir = outDirs["atac"],
      smpl = smplsNames)

  elif mode == "haystack":
    return expand(f"{{outDir}}/{{smpl}}/viralFrags.tsv",
      outDir = outDirs[outDirKey]["haystack"],
      smpl = smplsNames)

  else:
    sys.exit("Mode is not set correctly. Currently allowed values: 'namesort', 'haystack'")


rule all: 
  input:
    final_output(config["general"]["mode"])


rule cellranger_count:
  output:
     outDirs["atac"] + "/{smpl}/outs/possorted_bam.bam"
  
  params:
    outDir = outDirs["atac"],
    cellrangerAtac = config["paths"]["cellranger_exe"],
    fastqDir = config["paths"]["fastq_dir"],
    ref = config["paths"]["cellranger_ref_dir"]

  threads: nThreadsPerSample
  resources:
    mem_mb = nMemPerSample

  shell:
    """
    (
      cd {params.outDir}
      {params.cellrangerAtac} count \
        --id={wildcards.smpl} \
        --reference={params.ref} \
        --fastqs={params.fastqDir} \
        --sample={wildcards.smpl} \
        --localcores={threads} \
        --localmem={resources.mem_mb}
    )
    """

rule namesort:
  input:
    rules.cellranger_count.output

  output:
    outDirs["atac"] + "/{smpl}/outs/namesorted.bam"
  
  threads: nThreadsPerSample
  resources:
    mem_mb = nMemPerSample
  shell:
    "samtools sort -@ {threads} -n -o {output} -O bam {input}"
  

