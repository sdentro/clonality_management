#!/usr/bin/env python

import os, stat, sys, argparse
from path import path

LOG_DIR_NAME = "logs"
PROBLEMLOCI="/lustre/scratch110/sanger/sd11/Documents/GenomeFiles/battenberg_probloci/probloci.txt"
IMPUTEINFOFILE="/lustre/scratch110/sanger/sd11/Documents/GenomeFiles/battenberg_impute/impute_info.txt"
IMPUTE_EXE="impute2"
G1000_PREFIX="/lustre/scratch110/sanger/sd11/Documents/GenomeFiles/battenberg_1000genomesloci2012/1000genomesAlleles2012_chr"
G1000_PREFIX_AC="/lustre/scratch110/sanger/sd11/Documents/GenomeFiles/alleleCount_1000genomesloci2012/1000genomesloci2012_chr"
PLATFORM_GAMMA=1
PHASING_GAMMA=1
SEGMENTATION_GAMMA=10
CLONALITY_DIST_METRIC=0
ASCAT_DIST_METRIC=1
MIN_PLOIDY=1.6
MAX_PLOIDY=4.8
MIN_RHO=0.1
MIN_GOODNESS_OF_FIT=0.63
BALANCED_THRESHOLD=0.51
MIN_NORMAL_DEPTH=10

PIPE_DIR="/nfs/users/nfs_s/sd11/software/pipelines/battenberg_v1.0/"

RUN_SCRIPT = "RunCommands.sh"

SAMPLENAME_COL = 0
NORMAL_ID_COL = 1
NORMAL_BAM_COL = 2
TUMOUR_ID_COL = 3
TUMOUR_BAM_COL = 4
GENDER_COL = 5

def parseInputFile(filename):
    '''
    Function that parses a simple samplesheet
    
    Expects input of the following format
    samplename, normalbam, tumourbam, gender
    '''
    samples = dict()
    for line in open(filename, 'r'):
        
        if not line.startswith("#"):
            words = line.strip().split("\t")
            if words[SAMPLENAME_COL] in samples.keys():
                print("Samplename "+words[SAMPLENAME_COL]+" found more than once. Make sure samplenames are unique.")
                sys.exit(1)
            else:
                samples[words[SAMPLENAME_COL]] = {"normal_id": words[NORMAL_ID_COL], "normal": words[NORMAL_BAM_COL], "tumour_id": words[TUMOUR_ID_COL], "tumour": words[TUMOUR_BAM_COL], 'gender': words[GENDER_COL]}
            
    return samples

def generateBattenbergSNP6Pipeline(inputfile, pipe_dir, run_dir):
    '''
    Sets up a pipeline directory and parameters file for each sample mentioned in the inputfile
    '''
    samples = parseInputFile(inputfile)
    master_run_script = open(path.joinpath(run_dir, "RunCommands.sh"), 'w')
    
    for samplename,values in samples.iteritems():
        # create dir for tumour and logs
        path.joinpath(run_dir, values['tumour_id']).mkdir()
        path.joinpath(run_dir, values['tumour_id'], "logs").mkdir()
        
        # stick params file in there
        paramsFile = path.joinpath(run_dir, samplename, "params"+samplename+".txt")
        generateParams(paramsFile, samplename, values['normal_id'], values['normal'], values['tumour_id'], values['tumour'], values['gender'], pipe_dir, run_dir)
        
        # Append run command
        master_run_script.write(path.joinpath(pipe_dir, RUN_SCRIPT)+" "+paramsFile+"\n")
        
    master_run_script.close()
        
def generateParams(outfile, samplename, normal_id, normal_file, tumour_id, tumour_file, gender, pipe_dir, run_dir):
    '''
    Creates the parameters file that BB_snp6 uses
    '''
    fout = open(outfile, 'w')
    fout.write("RUN_DIR="+path.joinpath(run_dir, samplename)+"\n")
    fout.write("LOG_DIR="+path.joinpath(run_dir, samplename, LOG_DIR_NAME)+"\n")
    fout.write("PIPELINE_DIR="+pipe_dir+"\n")
    fout.write("TUMOURNAME="+tumour_id+"\n")
    fout.write("NORMALNAME="+normal_id+"\n")
    fout.write("TUMOURBAM="+tumour_file+"\n")
    fout.write("NORMALBAM="+normal_file+"\n")

    if gender == "male" or gender == "Male":
        is_male = "TRUE"
    elif gender == "Female" or gender == "female":
        is_male = "FALSE"
    else:
        print("Supplied gender for sample "+samplename+" not male or female")
        sys.exit(1)
        
    fout.write("IS_MALE="+is_male+"\n")
    
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
    
    fout.write("PROBLEMLOCI="+PROBLEMLOCI+"\n")
    fout.write("IMPUTEINFOFILE="+IMPUTEINFOFILE+"\n")
    fout.write("IMPUTE_EXE="+IMPUTE_EXE+"\n")
    fout.write("G1000_PREFIX="+G1000_PREFIX+"\n")
    fout.write("G1000_PREFIX_AC="+G1000_PREFIX_AC+"\n")

    fout.close()
    

def main(argv):
    parser = argparse.ArgumentParser(prog='Battenberg WGS',
                             formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", required=True, type=str, help="Input file containing: column of sample names, column of normal ids, column of normal bams, column of tumour ids, column of tumour bams, column of genders. Header must start with #.")
    parser.add_argument("-p", required=False, type=str, help="Directory where the pipeline is installed")
    parser.add_argument("-r", required=True, type=str, help="Directory where pipelines will be ran")
    
    parser.set_defaults(p=PIPE_DIR)
    args = parser.parse_args()
    
    generateBattenbergSNP6Pipeline(args.i, args.p, args.r)

if __name__ == '__main__':
    main(sys.argv[0:])
