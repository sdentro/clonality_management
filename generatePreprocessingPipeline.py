import os, stat, sys, argparse
from path import path
from util import merge_items

samplename = "PD7404a"
bam_file = "/lustre/scratch110/sanger/sd11/epitax/bam/PD7404a.bam"
bai_file = "/lustre/scratch110/sanger/sd11/epitax/bam/PD7404a.bam.bai"
vcf_file = "/lustre/scratch110/sanger/sd11/epitax/variants/filtered_vcf/PD7404a.filt.vcf.gz"
run_dir = "/nfs/users/nfs_s/sd11/repo/dirichlet_preprocessing/test/PD7404a/"

PIPE_DIR = "/nfs/users/nfs_s/sd11/repo/dirichlet_preprocessing"

CHROMS_FILE = "/nfs/users/nfs_s/sd11/repo/dirichlet_preprocessing/chroms_human.txt"
# Parse the CHROMS file into the values we need to create jobs
chroms = []
sex_chroms = []
for line in open(CHROMS_FILE, 'r'):
    words = line.strip().split("\t")
    if int(words[2]) == 0:
        chroms.append(words[0])
    else:
        sex_chroms.append(words[0])

no_aut_chroms = len(chroms)
no_chroms = len(chroms) + len(sex_chroms)


def generateBsubCmd(jobname, logdir, cmd, queue="normal", mem=1, depends=None, isArray=False, threads=None):
    '''
    Transforms the cmd into a bsub command with the supplied parameters.
    '''
    bcmd = merge_items(["bsub","-q", queue, "-J \""+jobname+"\""])
    
    if isArray:
        bcmd = merge_items([bcmd, "-o", path.joinpath(logdir, jobname)+".%J.%I.out", "-e", path.joinpath(logdir, jobname+".%J.%I.err")])
    else:
        bcmd = merge_items([bcmd, "-o", path.joinpath(logdir, jobname)+".%J.out", "-e", path.joinpath(logdir, jobname+".%J.err")])

    mem = str(mem)+"000"
    bcmd = merge_items([bcmd, "-M", mem, "-R", "'span[hosts=1] select[mem>" + mem + "] rusage[mem=" + mem + "]'"])

    if depends is not None:
        depends_str = map(lambda x: "done("+x+")", depends)
        depends_str = "&&".join(depends_str)    
        bcmd = merge_items([bcmd, "-w\""+depends_str+"\""])
        
                
    if threads is not None:
        bcmd = merge_items([bcmd, "-n", str(threads)])

    bcmd = merge_items([bcmd, "'"+cmd+"'"])

    return(bcmd)

def writeSimpleShellScript(rundir, scriptname, cmds):
    '''
    Creates a simple script with the commands specified in cmds contained within.
    This script works with jobarrays. 
    Note: It returns the status of the last run command.
    '''
    #scriptfile = path.joinpath(rundir, 'GetAlleleFrequenciesFromBAMByChromosome'+samplename+'.sh')
    scriptfile = path.joinpath(rundir, scriptname)
    samplecommands = open(scriptfile,'w')
    samplecommands.write('#$LSB_JOBINDEX\n')
    for item in cmds:
        samplecommands.write(item+"\n")

    samplecommands.write('exit $?\n')
    samplecommands.close()
    st = os.stat(scriptfile)
    os.chmod(scriptfile, st.st_mode | stat.S_IEXEC)
    
    return(scriptfile)

def createGenerateAFLociCmd(samplename, vcf_file, pipe_dir, run_dir):
    return(merge_items(["python", path.joinpath(pipe_dir, "dirichlet_preprocessing.py -c generateAFLoci"),
                        "-s", samplename,
                        "-o", samplename+".loci",
                        "-v", vcf_file,
                        "-r", run_dir,
                        "-p", pipe_dir]))
    
def createSplitLociCmd(samplename, loci_file, prefix, postfix, chroms_file, pipe_dir, run_dir):
    return(merge_items(["python", path.joinpath(pipe_dir, "dirichlet_preprocessing.py -c splitLociFile"),
                        "-s", samplename,
                        "--loci", loci_file,
                        "--prefix", prefix,
                        "--postfix", postfix,
                        "--chroms", chroms_file,
                        "-r", run_dir,
                        "-p", pipe_dir]))
    
def createGetAlleleFrequencyCmd(samplename, loci_file_prefix, bam_file, out_file_prefix, pipe_dir, run_dir):
    return(merge_items(["python", path.joinpath(pipe_dir, "dirichlet_preprocessing.py -c getAlleleFrequency"),
                        "-s", samplename,
                        "--bam", bam_file,
                        "--loci", loci_file_prefix+"${LSB_JOBINDEX}.txt",
                        "-o", out_file_prefix+"${LSB_JOBINDEX}.txt",
                        "-r", run_dir,
                        "-p", pipe_dir]))
    
def createConcatSplitFilesCmd(samplename, infile_list, outfile, haveHeader, run_dir, pipe_dir):
    cmd = ["python", path.joinpath(pipe_dir, "dirichlet_preprocessing.py -c concatSplitFiles"),
            "-s", samplename,
            "--files", merge_items(infile_list, sep=","),
            "-o", outfile,
            "-r", run_dir,
            "-p", pipe_dir]
    
    if haveHeader:
        cmd.append("--haveHeader")
    
    return(merge_items(cmd))
    
def createMutMutPhasingCmd(samplename, loci_file_prefix, out_file_prefix, bam_file, bai_file, max_distance, bb_dir, run_dir, pipe_dir):
    return(merge_items(["python", path.joinpath(pipe_dir, "dirichlet_preprocessing.py -c mutMutPhasing"),
                        "-s", samplename,
                        "--loci", loci_file_prefix+"${LSB_JOBINDEX}.txt",
                        "-o", out_file_prefix+"${LSB_JOBINDEX}.txt",
                        "--bam", bam_file,
                        "--bai", bai_file,
                        "--max_distance", str(max_distance),
                        "-b", bb_dir,
                        "-r", run_dir,
                        "-p", pipe_dir]))
    
def createMutCnPhasingCmd(samplename, loci_file_prefix, baf_file, hap_info_prefix, hap_info_suffix, outfile_prefix, bam_file, bai_file, max_distance, bb_dir, run_dir, pipe_dir):
    return(merge_items(["python", path.joinpath(pipe_dir, "dirichlet_preprocessing.py -c mutCNPhasing"),
                        "-s", samplename,
                        "--loci", loci_file_prefix+"${LSB_JOBINDEX}.txt",
                        "--phased_baf", baf_file,
                        "--hap_info", hap_info_prefix+"${LSB_JOBINDEX}"+hap_info_suffix,
                        "-o", outfile_prefix+"${LSB_JOBINDEX}.txt",
                        "--bam", bam_file,
                        "--bai", bai_file,
                        "-b", bb_dir,
                        "-r", run_dir,
                        "-p", pipe_dir]))
    
def createDpInputCmd(samplename, loci_file, allele_freq_file, subclone_file, rho_psi_file, mut_mut_phase_file, mut_cn_phase_file, gender, bb_dir, run_dir, pipe_dir):
    return(merge_items(["python", path.joinpath(pipe_dir, "dirichlet_preprocessing.py -c dpInput"),
                        "-s", samplename,
                        "--loci", loci_file,
                        "--all_freq", allele_freq_file,
                        "--subclones", subclone_file,
                        "--rhopsi", rho_psi_file,
                        "--mut_mut", mut_mut_phase_file,
                        "--mut_cn", mut_cn_phase_file,
                        "-x", gender,
                        "-o", samplename+"_allDirichletProcessInfo.txt",
                        "-b", bb_dir,
                        "-r", run_dir,
                        "-p", pipe_dir]))
    
def createDpIn2VcfCmd(vcf_file, dpIn_file, outfile, pipe_dir):
    return(merge_items(["python", path.joinpath(pipe_dir, "dpIn2vcf.py"),
                        "-v", vcf_file,
                        "-i", dpIn_file,
                        "-o", outfile]))
    
def dp_preprocessing_pipeline(samplename, vcf_file, bam_file, bai_file, baf_file, hap_info_prefix, hap_info_suffix, subclone_file, rho_psi_file, chroms_file, max_distance, gender, bb_dir, log_dir, pipe_dir, run_dir):
    '''
    Creates a list of commands that together form the preprocessing pipeline. It consists of 3 separate threads (a,b,c)
    that come together in the last step. 
        1) Get list of loci of interest from vcf file
        2) Split the loci per chromosome
        3a1) Obtain allele counts for each of the split loci files in parallel
        3a2) Concatenate the allele counts
        3b1) Perform mutation to mutation phasing for those pairs of mutations less then max_distance apart, per chromosome in parallel
        3b2) Concatenate the mutation to mutation phasing files
        3c1) Perform mutation to copynumber phasing for pairs less then max_distance apart, per chromosome in parallel
        3c2) Concatenate the mutation to copynumber phasing files
        4) Create the Dirichlet input file using all the above information
    '''
    runscript = path.joinpath(run_dir, "RunCommands_"+samplename+".sh")
    outf = open(runscript, 'w')
    
    # Generate the loci file from vcf
    cmd = createGenerateAFLociCmd(samplename, vcf_file, pipe_dir, run_dir)
    outf.write(generateBsubCmd("loci_"+samplename, log_dir, cmd, queue="normal", mem=1, depends=None, isArray=False) + "\n")
    
    # Split the loci file per chromosome
    cmd = createSplitLociCmd(samplename, samplename+".loci", samplename+"_loci_chr", ".txt", chroms_file, pipe_dir, run_dir)
    outf.write(generateBsubCmd("splitLoci_"+samplename, log_dir, cmd, queue="normal", mem=1, depends=["loci_"+samplename], isArray=False) + "\n")
    
    # Get the allele frequencies in parallel per chromosome
    cmd = createGetAlleleFrequencyCmd(samplename, samplename+"_loci_chr", bam_file, samplename+"_alleleFrequency_chr", pipe_dir, run_dir)
    writeSimpleShellScript(run_dir, "RunGetAlleleFrequency_"+samplename+".sh", [cmd])
    cmd = path.joinpath(run_dir, "RunGetAlleleFrequency_"+samplename+".sh")
    outf.write(generateBsubCmd("allCount_"+samplename+_arrayJobNameExt(no_chroms), log_dir, cmd, queue="normal", mem=1, depends=["splitLoci_"+samplename], isArray=True) + "\n")
    
    # Merge the counts together into a single file
    infile_list = [item[0]+str(item[1])+item[2] for item in zip([samplename+"_alleleFrequency_chr"]*no_chroms, range(1,no_chroms+1), [".txt"]*no_chroms)]
    cmd = createConcatSplitFilesCmd(samplename, infile_list, samplename+"_alleleFrequency.txt", True, run_dir, pipe_dir)
    outf.write(generateBsubCmd("concCounts_"+samplename, log_dir, cmd, queue="normal", mem=1, depends=["allCount_"+samplename], isArray=False) + "\n")
    
    '''
    ########################################################### Mut Mut Phasing ###########################################################
    '''
    cmd = createMutMutPhasingCmd(samplename, samplename+"_loci_chr", samplename+"_phasedmuts_chr", bam_file, bai_file, max_distance, bb_dir, run_dir, pipe_dir)
    writeSimpleShellScript(run_dir, "RunMutMutPhasing_"+samplename+".sh", [cmd])
    cmd = path.joinpath(run_dir, "RunMutMutPhasing_"+samplename+".sh")
    outf.write(generateBsubCmd("mmp_"+samplename+_arrayJobNameExt(no_chroms), log_dir, cmd, queue="normal", mem=2, depends=["splitLoci_"+samplename], isArray=True) + "\n")
    
    infile_list = [item[0]+str(item[1])+item[2] for item in zip([samplename+"_phasedmuts_chr"]*no_chroms, range(1,no_chroms+1), [".txt"]*no_chroms)]
    cmd = createConcatSplitFilesCmd(samplename, infile_list, samplename+"_phasedmuts.txt", True, run_dir, pipe_dir)
    outf.write(generateBsubCmd("concMMP_"+samplename, log_dir, cmd, queue="normal", mem=1, depends=["mmp_"+samplename], isArray=False) + "\n")
    
    '''
    ########################################################### Mut CN Phasing ###########################################################
    '''
    cmd = createMutCnPhasingCmd(samplename, samplename+"_loci_chr", baf_file, hap_info_prefix, hap_info_suffix, samplename+"_phased_mutcn_chr", bam_file, bai_file, max_distance, bb_dir, run_dir, pipe_dir)
    writeSimpleShellScript(run_dir, "RunMutCnPhasing_"+samplename+".sh", [cmd])
    cmd = path.joinpath(run_dir, "RunMutCnPhasing_"+samplename+".sh")
    # Note: We run this bit only for the autosomal chromosomes. The Y chrom can never be phased, while X is not as simple to do.
    outf.write(generateBsubCmd("mcp_"+samplename+_arrayJobNameExt(no_aut_chroms), log_dir, cmd, queue="normal", mem=2, depends=["splitLoci_"+samplename], isArray=True) + "\n")
    
    # Note: We run this bit only for the autosomal chromosomes. The Y chrom can never be phased, while X is not as simple to do.
    infile_list = [item[0]+str(item[1])+item[2] for item in zip([samplename+"_phased_mutcn_chr"]*no_chroms, range(1,no_aut_chroms+1), [".txt"]*no_chroms)]
    cmd = createConcatSplitFilesCmd(samplename, infile_list, samplename+"_phasedmutCN.txt", True, run_dir, pipe_dir)
    outf.write(generateBsubCmd("concMCP_"+samplename, log_dir, cmd, queue="normal", mem=1, depends=["mcp_"+samplename], isArray=False) + "\n")   
    
    '''
    ########################################################### Generate DP input ###########################################################
    '''
    cmd = createDpInputCmd(samplename, samplename+".loci", samplename+"_alleleFrequency.txt", subclone_file, rho_psi_file, samplename+"_phasedmuts.txt", samplename+"_phasedmutCN.txt", gender, bb_dir, run_dir, pipe_dir)
    outf.write(generateBsubCmd("dpIn_"+samplename, log_dir, cmd, queue="normal", mem=2, depends=["allCount_"+samplename, "concMMP_"+samplename, "concMCP_"+samplename], isArray=False) + "\n")
    
    '''
    ########################################################### DP input to VCF ###########################################################
    '''
    cmd = createDpIn2VcfCmd(vcf_file, path.joinpath(run_dir,samplename+"_allDirichletProcessInfo.txt"), path.joinpath(run_dir, samplename+".dpIn.vcf"), pipe_dir)
    outf.write(generateBsubCmd("dpIn2Vcf_"+samplename, log_dir, cmd, queue="normal", mem=2, depends=["dpIn_"+samplename], isArray=False) + "\n")
    
    outf.close()
    
    # Make executable
    st = os.stat(runscript)
    os.chmod(runscript, st.st_mode | stat.S_IEXEC)
    
    return(runscript)
    
def _readSampleSheet(infile):
    list_of_info = []
    for line in open(infile, "r"):
        words = line.strip().split("\t")
        list_of_info.extend([words])
    return list_of_info

def _arrayJobNameExt(no_chroms):
    return("[1-"+str(no_chroms)+"]")


def main(argv):
    parser = argparse.ArgumentParser(prog='Dirichlet_preprocessing pipeline',
                             formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-s", required=True, type=str, help="Sample sheet that contains a line per sample")
    parser.add_argument("-r", required=True, type=str, help="Directory where the pipeline will be run")
    
    # Parameters
    parser.add_argument("--min_baq", type=int, help="Minimum BAQ for a base to be included")
    parser.add_argument("--min_maq", type=int, help="Minimum MAQ for a base to be included")
    parser.add_argument("--max_distance", type=int, help="Maximum distance for a pair of mutations to be considered for phasing. Use when either mut_mut or mut_cn phasing")
    
    parser.set_defaults(min_baq=10, min_maq=10, max_distance=700, debug=False)
    
    args = parser.parse_args()
    
    # read in a samplesheet
    samples = _readSampleSheet(args.s)
    
    runscripts_sample = []
    for i in range(0,len(samples)):
        samplename = samples[i][0]
        print(samplename)
        vcf_file = samples[i][1]
        bam_file = samples[i][2]
        bai_file = samples[i][3]
        bb_dir = samples[i][4]
        gender = samples[i][5]
        baf_file = samples[i][6]
        subclone_file = samples[i][7]
        rho_psi_file = samples[i][8]
        hap_info_prefix = samples[i][9]
        hap_info_suffix = samples[i][10]
    
        run_dir = path.joinpath(args.r, samplename)
        log_dir = path.joinpath(run_dir, "logs")
        
        if not run_dir.exists():
            run_dir.mkdir()
            log_dir.mkdir()
        
        runscript = dp_preprocessing_pipeline(samplename=samplename, 
                                  vcf_file=vcf_file, 
                                  bam_file=bam_file, 
                                  bai_file=bai_file, 
                                  baf_file=baf_file, 
                                  hap_info_prefix=hap_info_prefix,
                                  hap_info_suffix=hap_info_suffix,
                                  subclone_file=subclone_file, 
                                  rho_psi_file=rho_psi_file, 
                                  chroms_file=CHROMS_FILE, 
                                  max_distance=args.max_distance, 
                                  gender=gender, 
                                  bb_dir=bb_dir, 
                                  log_dir=log_dir, 
                                  pipe_dir=PIPE_DIR, 
                                  run_dir=run_dir)
        runscripts_sample.append(runscript)
        
    # Create a master script that contains pointers to all sample specific runscripts
    scriptname = path.joinpath(args.r, "..","RunCommands.sh")
    runscript = open(scriptname, 'w')
    for item in runscripts_sample:
        print(item)
        runscript.write(item+"\n")
    runscript.close()
    
    # Make executable
    st = os.stat(scriptname)
    os.chmod(scriptname, st.st_mode | stat.S_IEXEC)
    
    
if __name__ == '__main__':
    main(sys.argv[0:])
        