#!/usr/bin/env python

'''
Script to generate a pure R BB pipeline for a table of samples
'''

import os, stat, sys, argparse
from path import path

#BB_PURE_R_SNP6 = "/nfs/users/nfs_s/sd11/repo/Battenberg/inst/example/battenberg_snp6.R"
BB_PURE_R_SNP6 = "/lustre/scratch112/sanger/cancer_external/DBGap/TCGA_phs000178.v8.p7/TCGA_pancan/pipelines/battenberg/battenberg_snp6_tcga_pancan.R"
BB_PURE_R_WGS = ""

QUEUE = "basement"

NUM_THREADS = 1
MEMORY = 20000

def generateBattenbergPureR(sample, normal, tumour, run_dir_sample, log_dir, pipe_script):
    '''
    Setup a pure R BB pipeline for the sample with its associated normal, tumour, run_dir and log_dir
    '''
    # Write a runscript
    runscript = path.joinpath(run_dir_sample, "submit.sh")
    outf = open(runscript, 'w')
    outf.write("TUMOURNAME="+sample+"\n")
    outf.write("NORMALCEL="+normal+"\n")
    outf.write("TUMOURCEL="+tumour+"\n")
    outf.write("RUN_DIR="+run_dir_sample+"\n")
    outf.write("bsub -n "+str(NUM_THREADS) + \
	       " -q "+QUEUE + \
               " -R\"select[mem>"+str(MEMORY)+"] rusage[mem="+str(MEMORY)+"] span[hosts=1]\" -M"+str(MEMORY)+" " + \
               " -o "+path.joinpath(log_dir, sample+".%J.out") + " -J"+sample + \
               " \"R CMD BATCH '--no-restore-data --no-save --args "+sample+" "+normal+" "+tumour+" "+run_dir_sample+" "+str(NUM_THREADS)+"' " + \
               pipe_script+" "+path.joinpath(log_dir, "battenberg_snp6."+sample+".Rout")+"\"")
    outf.close()

    # Make the runscript executable
    st = os.stat(runscript)
    os.chmod(runscript, st.st_mode | stat.S_IEXEC)
    
    return runscript

def main(argv):
    parser = argparse.ArgumentParser(prog='Generate a pure R Battenberg pipeline',
                             formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-s", required=True, type=str, help="Sample sheet that contains a line per sample")
    parser.add_argument("-r", required=True, type=str, help="Directory where the pipelines will run")
    parser.add_argument("--snp6", action="store_true", help="Set up a SNP6 pipeline")
    parser.set_defaults(snp6=True)
    args = parser.parse_args()
    
    run_dir = args.r
    log_dir = path.joinpath(run_dir, "logs")
    if not log_dir.exists(): log_dir.mkdir()
    
    runcommands = []    
    for line in open(args.s, 'r'):
        # Skip headers
        if line.startswith("#"): continue
        
        # Unpack the sample information
        words = line.strip().split("\t")
        sample = words[0]
        tumour = words[1]
        normal = words[2]
        
        # Create the run dir for this sample
        run_dir_sample = path.joinpath(run_dir, sample)
        if not run_dir_sample.exists(): run_dir_sample.mkdir()
        
        # Set the pipeline to be set up
        if args.snp6:
            pipe_script = BB_PURE_R_SNP6
        else:
            pipe_script = BB_PURE_R_WGS
        
        # Create the pipeline and store the path to the submit script
        runscript = generateBattenbergPureR(sample, normal, tumour, run_dir_sample, log_dir, pipe_script)
        runcommands.append(runscript)
        
    # Create a master file for easy submission
    runscript = path.joinpath(run_dir, "submit_all.sh")
    outf = open(runscript, 'w')
    for cmd in runcommands:
        outf.write(cmd+"\n")
    outf.close()
        
    # Make executable
    st = os.stat(runscript)
    os.chmod(runscript, st.st_mode | stat.S_IEXEC)
        

if __name__ == '__main__':
    main(sys.argv[0:])
