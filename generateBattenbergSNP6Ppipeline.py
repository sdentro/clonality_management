#!/usr/bin/env python

import os, stat, sys, argparse
from path import path


LOG_DIR_NAME = "logs"
MIN_COUNT = 1
PROB_LOCI_FILE = "NA"
USE_LOCI_FILE = "NA"
HETEROZYGOUS_FILTER = "NA"
FILL_IN_SNPS = 0
USE_TUMOUR_SNPS=0
USE_HETEROZYGOUS_SNPS_ONLY=0
HOMOZYGOUS_CAVEMAN_CALLS_FILE="NA"
PLATFORM_GAMMA=0.55
PHASING_GAMMA=1
SEGMENTATION_GAMMA=5
CLONALITY_DIST_METRIC=0
ASCAT_DIST_METRIC=1
MIN_PLOIDY=1.6
MAX_PLOIDY=4.8
MIN_RHO=0.1
MIN_GOODNESS_OF_FIT=0.63
BALANCED_THRESHOLD=0.51

PIPE_DIR="/nfs/users/nfs_s/sd11/software/pipelines/battenberg_snp6_v1.0/"
IMPUTEINFOFILE="/lustre/scratch110/sanger/sd11/Documents/GenomeFiles/battenberg_impute/impute_info.txt"
IMPUTE_EXE="impute2"
SNPPOS="/lustre/scratch112/sanger/cancer_external/DBGap/TCGA_phs000178.v8.p7/analysis/ASCAT/pvl/PRAD/SNPpos.txt"
GC_SNP6="/lustre/scratch112/sanger/cancer_external/DBGap/TCGA_phs000178.v8.p7/analysis/ASCAT/pvl/PRAD/GC_SNP6.txt"
ANNO_FILE="/lustre/scratch110/sanger/sd11/Documents/GenomeFiles/battenberg_snp6/GenomeWideSNP_6.na32.annot.subset.csv"
SNP6_REF_INFO_FILE="/lustre/scratch110/sanger/sd11/Documents/GenomeFiles/battenberg_snp6/snp6_ref_info_file.txt"
BIRDSEED_REPORT_FILE="birdseed.report.txt"
APT_PROBESET_GENOTYPE_EXE="~pvl/PennCNV/apt-1.12.0-20091012-i386-intel-linux/bin/apt-probeset-genotype"
APT_PROBESET_SUMMARIZE_EXE="~pvl/PennCNV/apt-1.12.0-20091012-i386-intel-linux/bin/apt-probeset-summarize"
NORM_GENO_CLUST_EXE="~pvl/PennCNV/gw6/bin/normalize_affy_geno_cluster.pl"

RUN_SCRIPT = "RunCommands2014farm3_SNP6.sh"

SAMPLENAME_COL = 0
NORMAL_COL = 1
TUMOUR_COL = 2
GENDER_COL = 3

def parseInputFile(filename):
    '''
    Function that parses a simple samplesheet
    
    Expects input of the following format
    samplename, tumourcel, normalcel, gender
    '''
    samples = dict()
    for line in open(filename, 'r'):
        
        if not line.startswith("#"):
            words = line.strip().split("\t")
            if words[SAMPLENAME_COL] in samples.keys():
                print("Samplename "+words[SAMPLENAME_COL]+" found more than once. Make sure samplenames are unique.")
                sys.exit(1)
            else:
                samples[words[SAMPLENAME_COL]] = {"normal": words[NORMAL_COL], "tumour": words[TUMOUR_COL], 'gender': words[GENDER_COL]}
            
    return samples

def generateBattenbergSNP6Pipeline(inputfile, pipe_dir, run_dir):
    '''
    Sets up a pipeline directory and parameters file for each sample mentioned in the inputfile
    '''
    samples = parseInputFile(inputfile)
    master_run_script = open(path.joinpath(run_dir, "RunCommands.sh"), 'w')
    
    for samplename,values in samples.iteritems():
        # create dir for samplename and logs
        path.joinpath(run_dir, samplename).mkdir()
        path.joinpath(run_dir, samplename, "logs").mkdir()
        
        # stick params file in there
        paramsFile = path.joinpath(run_dir, samplename, "params"+samplename+".txt")
        generateParams(paramsFile, samplename, values['normal'], values['tumour'], values['gender'], pipe_dir, run_dir)
        
        # Append run command
        master_run_script.write(path.joinpath(pipe_dir, RUN_SCRIPT)+" "+paramsFile+"\n")
        
    master_run_script.close()
        
def generateParams(outfile, samplename, normal_file, tumour_file, gender, pipe_dir, run_dir):
    '''
    Creates the parameters file that BB_snp6 uses
    '''
    fout = open(outfile, 'w')
    fout.write("RUN_DIR="+path.joinpath(run_dir, samplename)+"\n")
    fout.write("LOG_DIR="+path.joinpath(run_dir, samplename, LOG_DIR_NAME)+"\n")
    fout.write("PIPELINE_DIR="+pipe_dir+"\n")
    fout.write("TUMOURNAME="+samplename+"\n")
    fout.write("NORMALCEL="+normal_file+"\n")
    fout.write("TUMOURCEL="+tumour_file+"\n")

    if gender == "male" or gender == "Male":
        is_male = "TRUE"
    elif gender == "Female" or gender == "female":
        is_male = "FALSE"
    else:
        print("Supplied gender for sample "+samplename+" not male or female")
        sys.exit(1)
        
    fout.write("IS_MALE="+is_male+"\n")
    
    fout.write("HETEROZYGOUS_FILTER="+HETEROZYGOUS_FILTER+"\n")
    fout.write("MIN_COUNT="+str(MIN_COUNT)+"\n")
    fout.write("FILL_IN_SNPS="+str(FILL_IN_SNPS)+"\n")
    fout.write("USE_TUMOUR_SNPS="+str(USE_TUMOUR_SNPS)+"\n")
    fout.write("USE_HETEROZYGOUS_SNPS_ONLY="+str(USE_HETEROZYGOUS_SNPS_ONLY)+"\n")
    fout.write("HOMOZYGOUS_CAVEMAN_CALLS_FILE="+HOMOZYGOUS_CAVEMAN_CALLS_FILE+"\n")
    fout.write("PLATFORM_GAMMA="+str(PLATFORM_GAMMA)+"\n")
    fout.write("PHASING_GAMMA="+str(PHASING_GAMMA)+"\n")
    fout.write("SEGMENTATION_GAMMA="+str(SEGMENTATION_GAMMA)+"\n")
    fout.write("CLONALITY_DIST_METRIC="+str(CLONALITY_DIST_METRIC)+"\n")
    fout.write("ASCAT_DIST_METRIC="+str(ASCAT_DIST_METRIC)+"\n")
    fout.write("MIN_PLOIDY="+str(MIN_PLOIDY)+"\n")
    fout.write("MAX_PLOIDY="+str(MAX_PLOIDY)+"\n")
    fout.write("MIN_RHO="+str(MIN_RHO)+"\n")
    fout.write("MIN_GOODNESS_OF_FIT="+str(MIN_GOODNESS_OF_FIT)+"\n")
    fout.write("BALANCED_THRESHOLD="+str(BALANCED_THRESHOLD)+"\n")
    fout.write("PROB_LOCI_FILE="+PROB_LOCI_FILE+"\n")
    fout.write("USE_LOCI_FILE="+USE_LOCI_FILE+"\n")
    
    fout.write("IMPUTEINFOFILE="+IMPUTEINFOFILE+"\n")
    fout.write("IMPUTE_EXE="+IMPUTE_EXE+"\n")    
    fout.write("SNPPOS="+SNPPOS+"\n")
    fout.write("GC_SNP6="+GC_SNP6+"\n")
    fout.write("ANNO_FILE="+ANNO_FILE+"\n")
    fout.write("SNP6_REF_INFO_FILE="+SNP6_REF_INFO_FILE+"\n")
    fout.write("BIRDSEED_REPORT_FILE="+BIRDSEED_REPORT_FILE+"\n")
    fout.write("APT_PROBESET_GENOTYPE_EXE="+APT_PROBESET_GENOTYPE_EXE+"\n")
    fout.write("APT_PROBESET_SUMMARIZE_EXE="+APT_PROBESET_SUMMARIZE_EXE+"\n")
    fout.write("NORM_GENO_CLUST_EXE="+NORM_GENO_CLUST_EXE+"\n")
    
    fout.close()
    

def main(argv):
    parser = argparse.ArgumentParser(prog='Battenberg SNP6',
                             formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", required=True, type=str, help="Input file containing: column of sample names, column of normals, column of tumours, column of genders. Header must start with #.")
    parser.add_argument("-p", required=False, type=str, help="Directory where the pipeline is installed")
    parser.add_argument("-r", required=True, type=str, help="Directory where pipelines will be ran")
    
    parser.set_defaults(p=PIPE_DIR)
    args = parser.parse_args()
    
    generateBattenbergSNP6Pipeline(args.i, args.p, args.r)

if __name__ == '__main__':
    main(sys.argv[0:])
