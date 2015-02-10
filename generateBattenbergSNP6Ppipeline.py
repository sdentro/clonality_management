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

IMPUTEINFOFILE="/lustre/scratch110/sanger/sd11/Documents/GenomeFiles/battenberg_impute/impute_info.txt"
IMPUTE_EXE="/nfs/users/nfs_s/sd11/repo/battenberg_snp6/impute_v2.2.2_x86_64_static/impute2"
SNPPOS="/lustre/scratch112/sanger/cancer_external/DBGap/TCGA_phs000178.v8.p7/analysis/ASCAT/pvl/PRAD/SNPpos.txt"
GC_SNP6="/lustre/scratch112/sanger/cancer_external/DBGap/TCGA_phs000178.v8.p7/analysis/ASCAT/pvl/PRAD/GC_SNP6.txt"
KNOWN_SNPS_AUTOSOMES_PREFIX="/lustre/scratch110/sanger/dw9/haplotype_pipeline/impute/ALL_1000G_phase1integrated_feb2012_impute/ALL_1000G_phase1integrated_feb2012_chr"
KNOWN_SNPS_X_PAR1="/lustre/scratch110/sanger/dw9/haplotype_pipeline/impute/ALL_1000G_phase1integrated_feb2012_impute/ALL_1000G_phase1integrated_feb2012_chrX_PAR1_impute.legend"
KNOWN_SNPS_X_PAR2="/lustre/scratch110/sanger/dw9/haplotype_pipeline/impute/ALL_1000G_phase1integrated_feb2012_impute/ALL_1000G_phase1integrated_feb2012_chrX_PAR2_impute.legend"
KNOWN_SNPS_X_NONPAR="/lustre/scratch110/sanger/dw9/haplotype_pipeline/impute/ALL_1000G_phase1integrated_feb2012_impute/ALL_1000G_phase1integrated_feb2012_chrX_nonPAR_impute.legend"
ANNO_FILE="/lustre/scratch110/sanger/sd11/Documents/GenomeFiles/battenberg_snp6/GenomeWideSNP_6.na32.annot.subset.csv"

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
    fout.write("RUN_DIR="+run_dir+"\n")
    fout.write("LOG_DIR="+path.joinpath(run_dir, samplename, LOG_DIR_NAME)+"\n")
    fout.write("PIPELINE_DIR="+pipe_dir+"\n")
    fout.write("SAMPLENAME="+samplename+"\n")
    fout.write("NORMAL_CEL="+normal_file+"\n")
    fout.write("TUMOUR_CEL="+tumour_file+"\n")

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
    fout.write("KNOWN_SNPS_AUTOSOMES_PREFIX="+KNOWN_SNPS_AUTOSOMES_PREFIX+"\n")
    fout.write("KNOWN_SNPS_X_PAR1="+KNOWN_SNPS_X_PAR1+"\n")
    fout.write("KNOWN_SNPS_X_PAR2="+KNOWN_SNPS_X_PAR2+"\n")
    fout.write("KNOWN_SNPS_X_NONPAR="+KNOWN_SNPS_X_NONPAR+"\n")
    
    fout.close()
    

def main(argv):
    parser = argparse.ArgumentParser(prog='Battenberg SNP6',
                             formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", required=True, type=str, help="Input file containing list of samples")
    parser.add_argument("-p", required=True, type=str, help="Directory where the pipeline is installed")
    parser.add_argument("-r", required=True, type=str, help="Directory where pipelines will be ran")
    args = parser.parse_args()
    
    generateBattenbergSNP6Pipeline(args.i, args.p, args.r)

if __name__ == '__main__':
    main(sys.argv[0:])