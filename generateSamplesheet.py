#!/usr/bin/env python

import sys, argparse
from util import merge_items, read_list_of_items

def generateSamplesheet(samplename_file, tumour_file, tumourid_file, normal_file, normalid_file, gender_file, bb_file, variants_file, output_file):
    """
    Takes various single column files and generates a samplesheet. The i'th row of each of these files will be
    joined together (as if the Unix command line paste was called).
    """
    # Read in the various files
    samplenames = read_list_of_items(samplename_file)
    tumours = read_list_of_items(tumour_file)
    tumour_ids = read_list_of_items(tumourid_file)
    normals = read_list_of_items(normal_file)
    normal_ids = read_list_of_items(normalid_file)
    genders = read_list_of_items(gender_file)
    # BB is optional (could be ran after creation of this project. Set default placeholder if this is the case
    if not bb_file is None:
        bb = read_list_of_items(bb_file)
    else:
        bb = ['NA']*len(samplenames)
    variants = read_list_of_items(variants_file)
    
    # Write the output, joining line i from all vectors together
    outf = open(output_file, 'w')
    outf.write(merge_items(["#sample","tumour_id","tumour","normal_id","normal","bb_dir","gender","variants"], sep="\t")+"\n")
    for i in range(0, len(samplenames)):
        outf.write(merge_items([samplenames[i],
                                tumour_ids[i],
                                tumours[i],
                                normal_ids[i],
                                normals[i],
                                bb[i],
                                genders[i],
                                variants[i]], sep="\t")+"\n")
        
    outf.close()

def main(argv):
    parser = argparse.ArgumentParser(prog='generateSamplesheet',
                             formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-s", required=True, type=str, help="File containing list of samplenames directories")
    parser.add_argument("--bt", required=True, type=str, help="File containing list of tumour BAM files")
    parser.add_argument("--bn", required=True, type=str, help="File containing list of normal BAM files")
    parser.add_argument("--idt", required=True, type=str, help="File containing list of tumour IDs")
    parser.add_argument("--idn", required=True, type=str, help="File containing list of normal IDs")
    parser.add_argument("-v", required=True, type=str, help="File containing list of VCF files")
    parser.add_argument("-x", required=True, type=str, help="File containing list of genders for each samplename")
    parser.add_argument("-o", required=True, type=str, help="Output file")
    parser.add_argument("-b", required=True, type=str, help="Full path to file containing list of Battenberg directories")    
    
    args = parser.parse_args()
    generateSamplesheet(args.s, args.bt, args.idt, args.bn, args.idn, args.x, args.b, args.v, args.o)
    
if __name__ == '__main__':
    main(sys.argv[0:])