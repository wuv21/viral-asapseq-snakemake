import argparse
import csv
from collections import OrderedDict
import logging

def version():
    print('0.0.2')

'''
Takes a csv-formatted file of feature barcode names and sequences. If your CSV file contains a header, set header=True.

Usage:
tags = get_tags('feature_barcodes.csv', header=True)
'''    
def get_tags(filename: str, header: bool=False, logger: logging.Logger=None):
    with open(filename, mode='r') as csv_file:
        csv_reader = csv.reader(csv_file)
        tags = {}
        if header: 
            next(csv_reader) #skip header unless the no_header option is specified
            if logger is not None:
                logger.info('CSV includes header row. Skipping.')
        else:
            if logger is not None:
                logger.info('CSV does not include header row.')
        for row in csv_reader:
            tags[row[0].strip()] = row[1].strip()
    return tags

'''
From a dictionary of feature barcode names and sequences, return a dictionary of feature barcode names and sequences with all possible single nucleotide mismatches.

Usage:
mismatch_tags = make_mismatch_map(tags)
'''
def make_mismatch_map(FeatureDict: dict, logger: logging.Logger=None):
    odict = OrderedDict()
    feature_barcode_length = len(FeatureDict[[name for name in FeatureDict][0]])
    if logger is not None:
        logger.info(f'Feature Barcode Length: {feature_barcode_length}')
        logger.info(f'Read {len(FeatureDict)} Feature Barcodes')
    for name in FeatureDict:
        seq = str(FeatureDict[name])
        if logger is not None:
            logger.info(f'Name: {name}')
            logger.info(f'Seq: {seq}')
        odict[f'{name}-*-*'] = seq[:feature_barcode_length]
        barcode=list(str(seq)[:feature_barcode_length])
        for pos in range(feature_barcode_length):
            letter = seq[pos]
            muts = ['T', 'G', 'A', 'C']
            muts.remove(letter)
            for mut in range(3):
                odictval = ''.join(barcode)
                odictval[mut] = muts[mut]
                odict[f'{name}-{pos}-{mut + 1}'] = odictval
    return odict

'''
From a dictionary of feature barcode names and sequences with all possible single nucleotide mismatches, write a t2g file and a fasta file.
'''
def write_mismatch_map(tag_map: dict, mismatch_t2g_path: str, mismatch_fasta_path: str, logger: logging.Logger=None):
    tagmap_file = open(mismatch_t2g_path, 'w+')
    tagmap_fasta = open(mismatch_fasta_path, 'w+')
    for i in list(tag_map.keys()):
        if i[-4:]=='-*-*':
            #print(i[:-4]+'\t'+i[:-4]+'\t'+i[:-4])
            tagmap_file.write(i[:-4]+'\t'+i[:-4]+'\t'+i[:-4]+'\n')
            tagmap_fasta.write('>' + i[:-4] + '\n' +tag_map[i] + '\n')
        else:
            #print(i+'\t'+'-'.join(i.split('-')[:-2])+'\t'+'-'.join(i.split('-')[:-2]))
            tagmap_file.write(i+'\t'+'-'.join(i.split('-')[:-2])+'\t'+'-'.join(i.split('-')[:-2])+'\n')
            tagmap_fasta.write('>' + i + '\n' +tag_map[i] + '\n')
    tagmap_file.close()
    tagmap_fasta.close()
    
'''
wrapper function for make_mismatch_map and write_mismatch_map
'''
def kite_mismatch_maps(FeatureDict: dict, mismatch_t2g_path: str, mismatch_fasta_path: str, logger: logging.Logger=None):
    write_mismatch_map(make_mismatch_map(FeatureDict), mismatch_t2g_path, mismatch_fasta_path, logger=logger)
    if logger is not None:
        logger.info('The t2g and fasta files are now ready')

def main(args: argparse.Namespace, logger: logging.Logger=None):
    tags = get_tags(args.FeatureRefCSV, args.header, logger=logger)
    kite_mismatch_maps(tags, args.t2g, args.fa, logger=logger)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('FeatureRefCSV', help='path to csv file with whitelist Feature Barcodes')
    parser.add_argument('--t2g', help='path to output t2g file, default ./FeaturesMismatch.t2g', default='./FeaturesMismatch.t2g', type=str)
    parser.add_argument('--fa', help='path to output fasta file, default ./FeaturesMismatch.fa', default='./FeaturesMismatch.fa', type=str)
    parser.add_argument('--header', help='include this flag if FeatureRefCSV has a header, omit if not', action='store_true')
    parser.add_argument('-v', '--verbose', help='print processed Feature Barcodes', action='store_true')

    args = parser.parse_args()

    logger = logging.getLogger(__name__)

    if args.verbose:
        logger.setLevel(logging.WARNING)
    else:
        logger.setLevel(logging.INFO)
    
    main(args, logger=logger)
    