#!/usr/bin/env python3
import argparse
import csv
import os.path

'''
Read the samplesheet and return a list of samples

Usage:

with open("samplesheet.csv", "r") as csvfile:
  reader = csv.reader(csvfile)
  samples = readSamplesheet(reader)

Output:

[{"sample_name": "sample1_ADT", "sample_number": 1}, {"sample_name": "sample1_HTO", "sample_number": 1}, ...]
'''
def readSamplesheet(reader: csv.reader):
  samples = []
  skip = False
  foundSamples = False

  # Get the sample names and numbers, found after the [Data] section
  for row in reader:
    if row[0] == "[Data]":
      skip = True
      foundSamples = True
    elif skip:
      skip = False
    elif foundSamples:
      sample_name, sample_number = row[0], row[1]
      samples.append({"sample_name": sample_name, "sample_number": sample_number})

  return sorted(samples, key = lambda x: x["sample_name"])

'''
Identify the samples in the samplesheet

Usage:
identifySamples([{"sample_name": "sample1_ADT", "sample_number": 1}, {"sample_name": "sample1_HTO", "sample_number": 1}, ...])

Output:

{"adt": [{"sample_name": "sample1_ADT", "sample_number": 1}, ...], "hto": [{"sample_name": "sample1_HTO", "sample_number": 1}, ...], ...}
'''
def identifySamples(samples: list):
  splitSamples = {}
  for sample in samples:
    sample_type = sample["sample_name"].split("_")[-1]
    if sample_type not in splitSamples:
      splitSamples[sample_type] = []
    splitSamples[sample_type].append(sample)

  return splitSamples

'''
Group the samples in the samplesheet

Usage:

groupSamples({"adt": [{"sample_name": "sample1_ADT", "sample_number": 1}, ...], "hto": [{"sample_name": "sample1_HTO", "sample_number": 1}, ...], ...})

Output:

[{{"adt": "sample_name": "sample1_ADT", "sample_number": 1}, "hto": {"sample_name": "sample1_HTO", "sample_number": 1}}, ...]
'''
def groupSamples(samples: dict):
  groupedSamples = []
  for sample in samples:
    # I think this relies on the samples being sorted by sample name
    for i in range(len(samples[sample])):
      if i >= len(groupedSamples):
        groupedSamples.append({})
      groupedSamples[i][sample] = samples[sample][i]
  return groupedSamples

'''
Generate a partial config file for each sample in the samplesheet

Usage:

generateSample(adt={"sample_name": "sample1_ADT", "sample_number": 1}, hto={"sample_name": "sample1_HTO", "sample_number": 1})

Output:

- name: "sample1"
  adt_fastqs: ["sample1_ADT_S1_R1_001.fastq.gz", "sample1_ADT_S1_R3_001.fastq.gz"]
  hto_fastqs: ["sample1_HTO_S1_R1_001.fastq.gz", "sample1_HTO_S1_R3_001.fastq.gz"]
'''
def generateSample(**samples: dict):
  # Get name by truncating the postfix of any of the samples
  generalSampleName = "_".join(samples[list(samples.keys())[0]]["sample_name"].split("_")[0:-1])

  sampleConf = f'''\
    - name: "{generalSampleName}"'''
  
  for sample in samples:
    sampleName = samples[sample]["sample_name"]
    sampleNumber = samples[sample]["sample_number"]
    sampleConf += f'''\n      {sample.lower()}_fastqs: ["{sampleName}_S{sampleNumber}_R1_001.fastq.gz", "{sampleName}_S{sampleNumber}_R3_001.fastq.gz"]'''
  
  return sampleConf

'''
Main function
'''
def main(args: argparse.Namespace):
  with open(args.samplesheet, "r") as csvfile:
    reader = csv.reader(csvfile, delimiter = ",")
    samples = readSamplesheet(reader)
  
  identifiedSamples = identifySamples(samples)
  groupedSamples = groupSamples(identifiedSamples)

  generatedConfigs = []
  for sample in groupedSamples:
    generatedConfigs.append(generateSample(**sample))
  
  # Output all the generated configs
  print("\n".join(generatedConfigs))
      

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description = "Parser to convert samplesheet to config.yaml samples format")
  parser.add_argument("-s", "--samplesheet", required = True, help = "Samplesheet csv file path")

  args = parser.parse_args()

  if not os.path.exists(args.samplesheet):
    raise Exception("Samplesheet not found")

  main(args)
