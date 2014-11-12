import sys, argparse
from path import path
import numpy as np
from util import run_command, merge_items
from generateClonalityPipeline_util import read_sample_infile, match_sample_to_file

# python ~/repo/generateClonalityPipeline/generateAlleleFrequencyLoci.py -i /lustre/scratch110/sanger/sd11/epitax/samplelist.txt -v /lustre/scratch110/sanger/sd11/epitax/caveman_v3_variants/filtered_vcf/ -o /lustre/scratch110/sanger/sd11/epitax/haplotype/mutation_loci

def generateAlleleFrequencyLoci(infile, vcf_dir, out_dir, col_to_match):
    '''
    Generates loci files that serve as input for the GetAlleleFrequency scripts
    
    TODO: This now uses the vcf ref/alt columns. Perhaps the ref value should be fetched from a reference genome
    '''
    ss = read_sample_infile(infile)
    
    if(col_to_match == 'sample'):
        list_of_ids = ss.getSamplenames()
    elif(col_to_match == 'tumour'):
        list_of_ids = ss.getAllTumoursList()

    vcf_files = np.array(path(vcf_dir).listdir('*.vcf*'))
    mapping, unmapped = match_sample_to_file(list_of_ids, vcf_files)
    
    for sample,files in mapping.iteritems():
        if (len(files) == 1):
            print(sample)
            generateAlleleFrequencyLociSample(sample, files[0], out_dir)
        else:
            print("Found more than one VCF file for sample "+sample)
            
    for sample in unmapped:
        print("Found no VCF file for sample "+sample)
            
def generateAlleleFrequencyLociSample(sample, vcf_file, out_dir):
    m, e = run_command(merge_items(['zcat',vcf_file,'| grep -v ^\# | cut -f1,2,4,5 >',path.joinpath(out_dir,sample+'.loci')]))

def main(argv):
    parser = argparse.ArgumentParser(prog='GenerateAlleleFrequencyLoci',
                             formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", required=True, type=str, help="Full path to sample sheet")
    parser.add_argument("-v", required=True, type=str, help="Full path to directory where VCF files are stored")
    parser.add_argument("-o", required=True, type=str, help="Full path to directory where the loci should be saved")
    parser.add_argument("-c", required=True, type=str, help="Column in the sample sheet to use to match the VCF files")
    args = parser.parse_args()
    
    generateAlleleFrequencyLoci(args.i, args.v, args.o, args.c)

if __name__ == '__main__':
    main(sys.argv[0:])
