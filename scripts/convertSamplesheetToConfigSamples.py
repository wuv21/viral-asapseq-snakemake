import argparse
import csv
import os.path
import re


def generateSample(s_adt, s_hto):
  adt_sn = s_adt["sample"]
  hto_sn= s_hto["sample"]
  adt_n = s_adt["sample_n"]
  hto_n = s_hto["sample_n"]
  
  s = adt_sn.replace("_ADT", "") 

  sample = f'''\
    - name: {s}
      adt_fastqs: ["{adt_sn}_S{adt_n}_R1_001.fastq.gz", "{adt_sn}_S{adt_n}_R3_001.fastq.gz"]
      hto_fastqs: ["{hto_sn}_S{hto_n}_R1_001.fastq.gz", "{hto_sn}_S{hto_n}_R3_001.fastq.gz"]'''
  
  return(sample)

def main(args):
  samples = []
  skip = False
  foundSamples = False

  with open(args.samplesheet, "r") as csvfile:
    reader = csv.reader(csvfile, delimiter = ",")

    for row in reader:
      if row[0] == "[Data]":
        skip = True
        foundSamples = True
        continue

      elif skip:
        skip = False
        continue

      if foundSamples:
        samples.append(row[0])


  samples = [{"sample": s, "sample_n": i + 1} for i,s in enumerate(samples)]
  samples = sorted(samples, key = lambda x: x["sample"])
  
  for i in range(0, len(samples), 2):
    print(generateSample(samples[i], samples[i+1]))
      

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description = "Parser to convert samplesheet to config.yaml samples format")
  parser.add_argument("-s", "--samplesheet", required = True, help = "Samplesheet csv file path")

  args = parser.parse_args()

  if not os.path.exists(args.samplesheet):
    raise Exception("Samplesheet not found")

main(args)
