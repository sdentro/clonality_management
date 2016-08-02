#!/usr/bin/env python

import sys, argparse, os, stat
from path import path
# from clonalityPipelineConfig import IMPUTEFILESDIR, G1000LOCIDIR
from generateClonalityPipeline_util import read_sample_infile, read_basic_sample_infile, generateBsubCmd #, writeSimpleShellScript
from util import merge_items
from battenberg_util import bb_pipeline_config

'''
#################################################################################################################
CGP-IT BB setup 
#################################################################################################################
'''
BBSCRIPT = "/software/perl-5.16.3/bin/perl ~/repo/cgpBattenberg/perl/bin/battenberg.pl"
GENOME_INDEX = "/lustre/scratch116/casm/cgp/pancancer/reference/genome.fa.fai"
IGNORE_FILE = "/lustre/scratch116/casm/cgp/pancancer/reference/battenberg_full/ignore_contigs.txt"
PROTOCOL = "WGS"

'''
#################################################################################################################
Pipeline BB setup - WGS
#################################################################################################################
'''
PLATFORM_GAMMA_WGS=1
MIN_COUNT=10 # TODO: should this be removed from SNP6?

'''
#################################################################################################################
Pipeline BB setup - SNP6
#################################################################################################################
'''
USE_LOCI_FILE="NA"
HETEROZYGOUS_FILTER="NA"
FILL_IN_SNPS=0
USE_TUMOUR_SNPS=0
USE_HETEROZYGOUS_SNPS_ONLY=0
HOMOZYGOUS_CAVEMAN_CALLS_FILE="NA"
PLATFORM_GAMMA_SNP6=0.55

SNPPOS="/lustre/scratch112/sanger/cancer_external/DBGap/TCGA_phs000178.v8.p7/analysis/ASCAT/pvl/PRAD/SNPpos.txt"
GC_SNP6="/lustre/scratch112/sanger/cancer_external/DBGap/TCGA_phs000178.v8.p7/analysis/ASCAT/pvl/PRAD/GC_SNP6.txt"
ANNO_FILE="/lustre/scratch110/sanger/sd11/Documents/GenomeFiles/battenberg_snp6/GenomeWideSNP_6.na32.annot.subset.csv"
SNP6_REF_INFO_FILE="/lustre/scratch110/sanger/sd11/Documents/GenomeFiles/battenberg_snp6/snp6_ref_info_file.txt"
BIRDSEED_REPORT_FILE="birdseed.report.txt"
APT_PROBESET_GENOTYPE_EXE="~pvl/PennCNV/apt-1.12.0-20091012-i386-intel-linux/bin/apt-probeset-genotype"
APT_PROBESET_SUMMARIZE_EXE="~pvl/PennCNV/apt-1.12.0-20091012-i386-intel-linux/bin/apt-probeset-summarize"
NORM_GENO_CLUST_EXE="~pvl/PennCNV/gw6/bin/normalize_affy_geno_cluster.pl"

'''
#################################################################################################################
General options
#################################################################################################################
'''
PHASING_GAMMA=1
SEGMENTATION_GAMMA=10
CLONALITY_DIST_METRIC=0
ASCAT_DIST_METRIC=1
MIN_PLOIDY=1.6
MAX_PLOIDY=4.8
MIN_RHO=0.13
MAX_RHO=1.02
MIN_GOODNESS_OF_FIT=0.63
BALANCED_THRESHOLD=0.51
SEED=123
USE_SV_BREAKPOINTS_FILE=False
MAX_CN_STATE=250

IMPUTE_EXE='impute2' #/lustre/scratch110/sanger/sd11/epitax/battenberg/PD7404a/impute_v2.2.2_x86_64_static/impute2'
IMPUTEINFOFILE='/lustre/scratch116/casm/cgp/pancancer/reference/battenberg_full/impute/impute_info.txt'
# This is a little confusing, but the allele counter needs the loci while the later steps need the alleles
G1000_ALLELES_PREFIX="/lustre/scratch116/casm/cgp/pancancer/reference/battenberg_full/1000genomesloci/1000genomesloci2012_chr"
G1000_PREFIX="/lustre/scratch116/casm/cgp/pancancer/reference/battenberg_full/1000genomesloci/1000genomesAlleles2012_chr"
PROBLEMLOCIFILE='/lustre/scratch116/casm/cgp/pancancer/reference/battenberg_full/probloci.txt'
GCCORRECTPREFIX = "/lustre/scratch110/sanger/sd11/Documents/GenomeFiles/battenberg_wgs_gc_correction_1000g_v3/1000_genomes_GC_corr_chr_"

'''
#################################################################################################################
Versions
#################################################################################################################
'''
PIPE_BASE_DIR="/nfs/users/nfs_s/sd11/software/pipelines/"

PIPE_DIR_WGS_1_0=PIPE_BASE_DIR+'battenberg_v1.0'
PIPE_DIR_WGS_1_0_EXE=path.joinpath(PIPE_DIR_WGS_1_0, "RunCommands.sh")
PIPE_DIR_WGS_1_0_RERUN_EXE=path.joinpath(PIPE_DIR_WGS_1_0, "RunCommandsRerunFitCopynumber.sh")
PIPE_DIR_WGS_1_0_RERUN_MANUAL_EXE=path.joinpath(PIPE_DIR_WGS_1_0, "RunCommandsRerunFitCopynumberManual.sh")

PIPE_DIR_WGS_2_0_0=PIPE_BASE_DIR+"battenberg_v2.0.0/"
PIPE_DIR_WGS_2_0_0_EXE=path.joinpath(PIPE_DIR_WGS_2_0_0, "RunCommands.sh")
PIPE_DIR_WGS_2_0_0_RERUN_EXE=path.joinpath(PIPE_DIR_WGS_2_0_0, "RunCommandsRerunFitCopynumber.sh")
PIPE_DIR_WGS_2_0_0_RERUN_MANUAL_EXE=path.joinpath(PIPE_DIR_WGS_2_0_0, "RunCommandsRerunFitCopynumberManual.sh")

PIPE_DIR_WGS_DEV="/nfs/users/nfs_s/sd11/repo/battenberg"
PIPE_DIR_WGS_DEV_EXE=path.joinpath(PIPE_DIR_WGS_DEV, "RunCommands.sh")
PIPE_DIR_WGS_DEV_RERUN_EXE=path.joinpath(PIPE_DIR_WGS_DEV, "RunCommandsRerunFitCopynumber.sh")
PIPE_DIR_WGS_DEV_RERUN_MANUAL_EXE=path.joinpath(PIPE_DIR_WGS_DEV, "RunCommandsRerunFitCopynumberManual.sh")

PIPE_DIR_SNP6_1_0=PIPE_BASE_DIR+'battenberg_snp6_v1.0'
PIPE_DIR_SNP6_1_0_EXE=path.joinpath(PIPE_DIR_SNP6_1_0, "RunCommands2014farm3_SNP6.sh")
PIPE_DIR_SNP6_1_0_RERUN_EXE=path.joinpath(PIPE_DIR_SNP6_1_0, "RunCommandsRerunFitCopynumber.sh")
PIPE_DIR_SNP6_1_0_RERUN_MANUAL_EXE=path.joinpath(PIPE_DIR_SNP6_1_0, "RunCommandsRerunFitCopynumberManual.sh")

PIPE_DIR_SNP6_2_0_0=PIPE_BASE_DIR+"battenberg_v2.0.0/"
PIPE_DIR_SNP6_2_0_0_EXE=path.joinpath(PIPE_DIR_SNP6_2_0_0, "Battenberg_SNP6_LSF.sh")
PIPE_DIR_SNP6_2_0_0_RERUN_EXE=path.joinpath(PIPE_DIR_SNP6_2_0_0, "RunCommandsRerunFitCopynumber.sh")
PIPE_DIR_SNP6_2_0_0_RERUN_MANUAL_EXE=path.joinpath(PIPE_DIR_SNP6_2_0_0, "RunCommandsRerunFitCopynumberManual.sh")

PIPE_DIR_SNP6_DEV="/nfs/users/nfs_s/sd11/repo/battenberg"
PIPE_DIR_SNP6_DEV_EXE=path.joinpath(PIPE_DIR_WGS_DEV, "Battenberg_SNP6_LSF.sh")
PIPE_DIR_SNP6_DEV_RERUN_EXE=path.joinpath(PIPE_DIR_WGS_DEV, "RunCommandsRerunFitCopynumber.sh")
PIPE_DIR_SNP6_DEV_RERUN_MANUAL_EXE=path.joinpath(PIPE_DIR_WGS_DEV, "RunCommandsRerunFitCopynumberManual.sh")

'''
#################################################################################################################
CGP IT BB version
#################################################################################################################
'''
def createAlleleCountCmd(bb_conf):
	cmd = [BBSCRIPT, bb_conf.getStandardOptions_cgpBB(), bb_conf.getThreadsOption_cgpBB(),
		"-p", "allelecount"]
	return(merge_items(cmd))

def createBafLogCmd(bb_conf):
	cmd = [BBSCRIPT, bb_conf.getStandardOptions_cgpBB(), bb_conf.getThreadsOption_cgpBB(),
		"-p", "baflog"]
	return(merge_items(cmd))

def createImputeFromBafCmd(bb_conf):
	cmd = [BBSCRIPT, bb_conf.getStandardOptions_cgpBB(), bb_conf.getThreadsOption_cgpBB(),
		"-p", "imputefromaf"]
	return(merge_items(cmd))	

def createImputeCmd(bb_conf):
	cmd = [BBSCRIPT, bb_conf.getStandardOptions_cgpBB(), bb_conf.getThreadsOption_cgpBB(),
		"-p", "impute"]
	return(merge_items(cmd))	

def createCombineImputeCmd(bb_conf):
	cmd = [BBSCRIPT, bb_conf.getStandardOptions_cgpBB(), bb_conf.getThreadsOption_cgpBB(),
		"-p", "combineimpute"]
	return(merge_items(cmd))	

def createHaplotypeBafsCmd(bb_conf):
	cmd = [BBSCRIPT, bb_conf.getStandardOptions_cgpBB(), bb_conf.getThreadsOption_cgpBB(),
		"-p", "haplotypebafs"]
	return(merge_items(cmd))	

def createCleanupPostBafCmd(bb_conf):
	cmd = [BBSCRIPT, bb_conf.getStandardOptions_cgpBB(), bb_conf.getThreadsOption_cgpBB(),
		"-p", "cleanuppostbaf"]
	return(merge_items(cmd))

def createPlotHaplotypesCmd(bb_conf):
	cmd = [BBSCRIPT, bb_conf.getStandardOptions_cgpBB(), bb_conf.getThreadsOption_cgpBB(),
		"-p", "plothaplotypes"]
	return(merge_items(cmd))

def createCombineBafsCmd(bb_conf):
	cmd = [BBSCRIPT, bb_conf.getStandardOptions_cgpBB(),
		"-p", "combinebafs"]
	return(merge_items(cmd))

def createSegmentPhasedCmd(bb_conf):
	cmd = [BBSCRIPT, bb_conf.getStandardOptions_cgpBB(),
		"-p", "segmentphased"]
	return(merge_items(cmd))

def createFitcnCmd(bb_conf):
	cmd = [BBSCRIPT, bb_conf.getStandardOptions_cgpBB(),
		"-p", "fitcn"]
	return(merge_items(cmd))

def createSubclonesCmd(bb_conf):
	cmd = [BBSCRIPT, bb_conf.getStandardOptions_cgpBB(),
		"-p", "subclones"]
	return(merge_items(cmd))

def createFinaliseCmd(bb_conf):
	cmd = [BBSCRIPT, bb_conf.getStandardOptions_cgpBB(),
		"-p", "finalise"]
	return(merge_items(cmd))

def generateBattenbergPipeline_CGPIT(bb_conf):
	samplename = bb_conf.get_samplename()
	log_dir = bb_conf.get_log_dir()
	threads = bb_conf.get_threads()
	
	runscript = path.joinpath(bb_conf.get_run_dir(), "RunCommands_"+samplename+".sh")
	outf = open(runscript, 'w')
	
	cmd = createAlleleCountCmd(bb_conf)
	outf.write(generateBsubCmd("loci_"+samplename, log_dir, cmd, queue="normal", mem=4, depends=None, isArray=False, threads=threads) + "\n")
	
	cmd = createBafLogCmd(bb_conf)
	outf.write(generateBsubCmd("baflog_"+samplename, log_dir, cmd, queue="normal", mem=28, depends=["loci_"+samplename], isArray=False, threads=threads) + "\n")
	
	cmd = createImputeFromBafCmd(bb_conf)
	outf.write(generateBsubCmd("imputebaf_"+samplename, log_dir, cmd, queue="normal", mem=7, depends=["loci_"+samplename], isArray=False, threads=threads) + "\n")

	cmd = createImputeCmd(bb_conf)
	outf.write(generateBsubCmd("impute_"+samplename, log_dir, cmd, queue="long", mem=25, depends=["imputebaf_"+samplename], isArray=False, threads=threads) + "\n")

	cmd = createCombineImputeCmd(bb_conf)
	outf.write(generateBsubCmd("combineimpute_"+samplename, log_dir, cmd, queue="long", mem=25, depends=["impute_"+samplename], isArray=False, threads=threads) + "\n")

	cmd = createHaplotypeBafsCmd(bb_conf)
	outf.write(generateBsubCmd("haplotypebafs_"+samplename, log_dir, cmd, queue="normal", mem=12, depends=["combineimpute_"+samplename], isArray=False, threads=threads) + "\n")
	
	cmd = createCleanupPostBafCmd(bb_conf)
	outf.write(generateBsubCmd("cleanuppostbaf_"+samplename, log_dir, cmd, queue="normal", mem=12, depends=["haplotypebafs_"+samplename], isArray=False, threads=threads) + "\n")
	
	cmd = createPlotHaplotypesCmd(bb_conf)
	outf.write(generateBsubCmd("plothaplotypes_"+samplename, log_dir, cmd, queue="normal", mem=1, depends=["cleanuppostbaf_"+samplename], isArray=False, threads=threads) + "\n")

	cmd = createCombineBafsCmd(bb_conf)
	outf.write(generateBsubCmd("combinebafs_"+samplename, log_dir, cmd, queue="normal", mem=4, depends=["plothaplotypes_"+samplename], isArray=False, threads=None) + "\n")

	cmd = createSegmentPhasedCmd(bb_conf)
	outf.write(generateBsubCmd("segmentphased_"+samplename, log_dir, cmd, queue="normal", mem=4, depends=["combinebafs_"+samplename], isArray=False, threads=None) + "\n")
	
	cmd = createFitcnCmd(bb_conf)
	outf.write(generateBsubCmd("fitcn_"+samplename, log_dir, cmd, queue="normal", mem=16, depends=["segmentphased_"+samplename, "baflog_"+samplename], isArray=False, threads=None) + "\n")
	
	cmd = createSubclonesCmd(bb_conf)
	outf.write(generateBsubCmd("subclones_"+samplename, log_dir, cmd, queue="normal", mem=16, depends=["fitcn_"+samplename], isArray=False, threads=None) + "\n")
	
	cmd = createFinaliseCmd(bb_conf)
	outf.write(generateBsubCmd("finalise_"+samplename, log_dir, cmd, queue="normal", mem=2, depends=["subclones_"+samplename], isArray=False, threads=None) + "\n")
	
	outf.close()
	
	# Make executable
	st = os.stat(runscript)
	os.chmod(runscript, st.st_mode | stat.S_IEXEC)
	
	return(runscript)

'''
#################################################################################################################
Regular LSF pipeline versions
#################################################################################################################
'''
def generateBattenbergPipeline_WGS(bb_config, pipe_exe, pipe_rerun_exe, pipe_rerun_manual_exe, no_allelecount):
	# Write the params file
	if no_allelecount:
		allecount_param="1"
	else:
		allecount_param="0"
	
	config_file = path.joinpath(bb_config.get_run_dir(), "params"+bb_config.get_samplename()+".txt")
	bb_config.generateParamsFile_WGS(config_file)
	
	# Return run commands
	return (pipe_exe+' '+config_file+' '+allecount_param), (pipe_rerun_exe+' '+config_file+' '+allecount_param), (pipe_rerun_manual_exe+' '+config_file+' '+allecount_param)
	
def generateBattenbergPipeline_SNP6(bb_config, pipe_exe, pipe_rerun_exe, pipe_rerun_manual_exe):
	# Write the params file
	config_file = path.joinpath(bb_config.get_run_dir(), "params"+bb_config.get_samplename()+".txt")
	bb_config.generateParamsFile_SNP6(config_file)
	
	# Return run commands
	return (pipe_exe+' '+config_file), (pipe_rerun_exe+' '+config_file), (pipe_rerun_manual_exe+' '+config_file)
		
'''
#################################################################################################################
Main that strings it all together
#################################################################################################################
'''		
def main(argv):
	parser = argparse.ArgumentParser(prog='GenerateBattenbergPipeline',
							 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("-i", type=str, required=True, help="Input file containing at least four columns: samplename, tumour_file, normal_file, gender. If this is a samplesheet, also specify --ss.")
	parser.add_argument("-r", type=str, required=True, help="Full path to a directory where the pipelines will be run")
	
	parser.add_argument("--ss", action="store_true", help="Input file is a samplesheet")
	
	# SNP6 and WGS options
	parser.add_argument("--type", choices=["wgs", "snp6"], type=str, help='Type of pipeline to be set up')
	parser.add_argument("-v", "--version", type=str, help='Version of this pipeline to set up')
	
	# General options
	parser.add_argument("--imputeinfofile", help='Path to impute info file')
	parser.add_argument("--impute_exe", help='Path to impute exe')
	parser.add_argument("--g1000_prefix_alleles", help='Prefix to 1000 Genomes alleles reference files')
	parser.add_argument("--g1000_prefix_loci", help='Prefix to 1000 Genomes loci reference files')
	parser.add_argument("--gc_correction_prefix", help='Prefix to GC correction reference files')
	parser.add_argument("--phasing_gamma", help="Phasing gamma parameter")
	parser.add_argument("--segmentation_gamma", help="Phasing segmentation parameter")
	parser.add_argument("--clonality_dist_metric", help="Type of distance metric used when fitting the clonal copy number profile")
	parser.add_argument("--ascat_dist_metric", help="Type of distance metric used when fitting the ASCAT profile, before BB specific fitting")
	parser.add_argument("--min_ploidy", help="Minimum ploidy to be considered")
	parser.add_argument("--max_ploidy", help="Maximum ploidy to be considered")
	parser.add_argument("--min_rho", help="Minimum rho (cellularity) to be considered")
	parser.add_argument("--max_rho", help="Maximum rho (cellularity) to be considered")
	parser.add_argument("--min_goodness_of_fit", help="Minimum goodness of fit to be allowed")
	parser.add_argument("--balanced_threshold", help="TODO") # TODO: explain this option 
	parser.add_argument("--problemlocifile", help="File containing problem loci")
	parser.add_argument("--seed", help='Seed to use')
	parser.add_argument("--use_sv_breakpoints_file", action="store_true", help='Supply when BB should expect a file with breakpoints in its working directory')
	parser.add_argument("--max_cn_state", help='Maximum copy number state to allow before masking is performed. Anything above this state is masked out')
	parser.add_argument("--no_allelecount", action="store_true", help="Supply when not to run the allelecounter, for now only works with WGS and not with cgpBB")
	
	# SNP6 options	
	parser.add_argument("--min_count", type=int, help="Minimum read count") # TODO: shouldn't this be added to the WGS version?? => Added to WGS, remove from SNP6 now?
	parser.add_argument("--fill_in_snps", type=int, choices=[0,1], help="Fill in SNPs, should be either 0 or 1") # TODO: still used?
	parser.add_argument("--heterozygous_filter", type=str, help="Heterozygous filter to be used") # TODO: still used?
	parser.add_argument("--use_tumour_snps", type=int, choices=[0,1], help="Use tumour SNPs, either 0 or 1") # TODO: still used?
	parser.add_argument("--use_het_snps_only", type=int, choices=[0,1], help="Use heterozygous SNPs only, either 0 or 1") # TODO: still used?
	parser.add_argument("--hom_caveman_snps", type=str, help="File containing Caveman called homozygous SNPs") # TODO: still used?
	parser.add_argument("--use_loci_file", type=str, help="File containing loci locations") # TODO: still used?
	parser.add_argument("--snppos", type=str, help="File SNP6 SNP locations")
	parser.add_argument("--gc_snp6", type=str, help="File SNP6 SNP GC content information")
	parser.add_argument("--anno", type=str, help="File with SNP annotation information")
	parser.add_argument("--snp6_ref_info_file", type=str, help="Full path to file with SNP6 reference info") # TODO: better explanation
	parser.add_argument("--birdseed_report_file", type=str, help="Name of birdseed output file")
	parser.add_argument("--apt_probeset_geno_exe", type=str, help="Full path to apt_probeset_genotype exe from AFFY tools")
	parser.add_argument("--apt_probeset_summ_exe", type=str, help="Full path to apt_probeset_summarise exe from AFFY tools")
	parser.add_argument("--norm_geno_clust_exe", type=str, help="Full path to normalise_genotype_clusters exe from PennCNV")
	
	# cgpBB options
	parser.add_argument("-t", type=int, help="Number of threads to use")
	parser.add_argument("--genome_index", type=str, help="Full path to a reference genome index FAI file")
	parser.add_argument("--protocol", type=str, choices=["WGS"], help="Sequencing protocol used")
	parser.add_argument("--ignore_file", type=str, help="File with chromosomes to ignore")
	
	parser.set_defaults(ss=False, t=1, version="2.0.0", type="wgs", \
					imputeinfofile=IMPUTEINFOFILE, impute_exe=IMPUTE_EXE, g1000_prefix_alleles=G1000_PREFIX, phasing_gamma=PHASING_GAMMA, \
					segmentation_gamma=SEGMENTATION_GAMMA, clonality_dist_metric=CLONALITY_DIST_METRIC, ascat_dist_metric=ASCAT_DIST_METRIC, \
					min_ploidy=MIN_PLOIDY, max_ploidy=MAX_PLOIDY, min_rho=MIN_RHO, max_rho=MAX_RHO, min_goodness_of_fit=MIN_GOODNESS_OF_FIT, \
					balanced_threshold=BALANCED_THRESHOLD, problemlocifile=PROBLEMLOCIFILE, min_count=MIN_COUNT, fill_in_snps=FILL_IN_SNPS, \
					heterozygous_filter=HETEROZYGOUS_FILTER, use_tumour_snps=USE_TUMOUR_SNPS, use_het_snps_only=USE_HETEROZYGOUS_SNPS_ONLY, \
					hom_caveman_snps=HOMOZYGOUS_CAVEMAN_CALLS_FILE, use_loci_file=USE_LOCI_FILE, snppos=SNPPOS, gc_snp6=GC_SNP6, anno=ANNO_FILE, \
					snp6_ref_info_file=SNP6_REF_INFO_FILE, birdseed_report_file=BIRDSEED_REPORT_FILE, apt_probeset_geno_exe=APT_PROBESET_GENOTYPE_EXE, \
					apt_probeset_summ_exe=APT_PROBESET_SUMMARIZE_EXE, norm_geno_clust_exe=NORM_GENO_CLUST_EXE, genome_index=GENOME_INDEX, \
					protocol=PROTOCOL, ignore_file=IGNORE_FILE, g1000_prefix_loci=G1000_ALLELES_PREFIX, seed=SEED, max_cn_state=MAX_CN_STATE, \
					use_sv_breakpoints_file=USE_SV_BREAKPOINTS_FILE, no_allelecount=False, gc_correction_prefix=GCCORRECTPREFIX)
	args = parser.parse_args()
	
	'''
	#######################################################################
	Setting Gamma platform parameter according to BB type
	#######################################################################
	'''
	if args.type=="wgs" or args.type=="WGS":
		platform_gamma = PLATFORM_GAMMA_WGS
	elif args.type=="snp6" or args.type=="SNP6":
		platform_gamma = PLATFORM_GAMMA_SNP6
	else:
		print("Unknown BB type supplied, do not know how to set the platform gamma parameter")
		sys.exit(1)
		
	'''
	#######################################################################
	Checking for BB type and version compatibility
	#######################################################################
	'''
	if (args.type=="wgs" or args.type=="WGS"):
		if (args.version=="1.0"):
			pipe_dir = PIPE_DIR_WGS_1_0
			pipe_exe = PIPE_DIR_WGS_1_0_EXE
			pipe_rerun_exe = PIPE_DIR_WGS_1_0_RERUN_EXE
			pipe_rerun_manual_exe = PIPE_DIR_WGS_1_0_RERUN_MANUAL_EXE
		elif (args.version=="1.1"):
			print("not implemented")
		elif (args.version=="2.0.0"):
			pipe_dir = PIPE_DIR_WGS_2_0_0
			pipe_exe = PIPE_DIR_WGS_2_0_0_EXE
			pipe_rerun_exe = PIPE_DIR_WGS_2_0_0_RERUN_EXE
			pipe_rerun_manual_exe = PIPE_DIR_WGS_2_0_0_RERUN_MANUAL_EXE
		elif (args.version=="dev"):
			pipe_dir = PIPE_DIR_WGS_DEV
			pipe_exe = PIPE_DIR_WGS_DEV_EXE
			pipe_rerun_exe = PIPE_DIR_WGS_DEV_RERUN_EXE
			pipe_rerun_manual_exe = PIPE_DIR_WGS_DEV_RERUN_MANUAL_EXE
		elif (args.version=="cgp"):
			# Set bogus values, not needed for this version
			pipe_dir = ""
			pipe_exe = ""
			pipe_rerun_exe = ""
			pipe_rerun_manual_exe = ""
		else:
			print("Unsupported BB WGS version supplied")
			sys.exit(1)
	elif (args.type=="snp6" or args.type=="SNP6"):
		if (args.version=="1.0"):
			pipe_dir = PIPE_DIR_SNP6_1_0
			pipe_exe = PIPE_DIR_SNP6_1_0_EXE
			pipe_rerun_exe = PIPE_DIR_SNP6_1_0_RERUN_EXE
			pipe_rerun_manual_exe = PIPE_DIR_SNP6_1_0_RERUN_MANUAL_EXE
		elif (args.version=="2.0.0"):
			pipe_dir = PIPE_DIR_SNP6_2_0_0
			pipe_exe = PIPE_DIR_SNP6_2_0_0_EXE
			pipe_rerun_exe = PIPE_DIR_SNP6_2_0_0_RERUN_EXE
			pipe_rerun_manual_exe = PIPE_DIR_SNP6_2_0_0_RERUN_MANUAL_EXE
		elif (args.version=="dev"):
			pipe_dir = PIPE_DIR_SNP6_DEV
			pipe_exe = PIPE_DIR_SNP6_DEV_EXE
			pipe_rerun_exe = PIPE_DIR_SNP6_DEV_RERUN_EXE
			pipe_rerun_manual_exe = PIPE_DIR_SNP6_DEV_RERUN_MANUAL_EXE
		else:
			print("Unsupported BB SNP6 version supplied")
			sys.exit(1)
			
	runscripts = []
	'''
	#######################################################################
	Read in input data
	#######################################################################
	'''
	if (args.ss):
		# Read in samplesheet
		ss = read_sample_infile(args.i)
	else:
		ss = read_basic_sample_infile(args.i, args.r)
	
	# For every entry:
	for sample in ss.getSamplenames():
		print(sample)
		
		tn_pair = ss.getTumour2NormalPairingBam(sample)
		for tb,nb in tn_pair:
			tumourid = ss.getIdByTumourBam(tb)
			normalid = ss.getNormals(sample)[0] # TODO: what happens when multiple normals?
			run_dir = path.joinpath(args.r, tumourid)
			log_dir = path.joinpath(args.r, tumourid, "logs")
			if not run_dir.exists(): run_dir.makedirs()
			if not log_dir.exists(): log_dir.makedirs()
			
			'''
			#######################################################################
			Create the BB config master object
			#######################################################################
			'''
			bb_conf = bb_pipeline_config(pipe_type=args.type, pipe_version=args.version, pipe_dir=pipe_dir, \
									tumour_file=tb, normal_file=nb, run_dir=run_dir, log_dir=log_dir, samplename=sample, gender=ss.getSex(sample), \
									tumour_id=tumourid, normal_id=normalid, impute_info=args.imputeinfofile, impute_exe=args.impute_exe, \
									g1000_loci_dir=args.g1000_prefix_loci, platform_gamma=platform_gamma, phasing_gamma=args.phasing_gamma, \
									segmentation_gamma=args.segmentation_gamma, clonality_dist_metric=args.clonality_dist_metric, min_count=args.min_count, \
									ascat_dist_metric=args.ascat_dist_metric, min_ploidy=args.min_ploidy, max_ploidy=args.max_ploidy, max_rho=args.max_rho, \
									min_rho=args.min_rho, min_goodness_of_fit=args.min_goodness_of_fit, balanced_threshold=args.balanced_threshold, \
									prob_loci_file=args.problemlocifile, heterozygous_filter=args.heterozygous_filter, use_tumour_snps=args.use_tumour_snps, \
									use_het_snps_only=args.use_het_snps_only, hom_caveman_file=args.hom_caveman_snps, use_loci_file=args.use_loci_file, \
									snppos_file=args.snppos, gc_snp6_file=args.gc_snp6, snp6_anno_file=args.anno, snp6_ref_info_file=args.snp6_ref_info_file, \
									birdseed_report_file=args.birdseed_report_file, apt_probeset_geno_exe=args.apt_probeset_geno_exe, \
									apt_probeset_summ_exe=args.apt_probeset_summ_exe, norm_geno_clust_exe=args.norm_geno_clust_exe, \
									fill_in_snps=args.fill_in_snps, threads=args.t, genome_index=args.genome_index, protocol=args.protocol, \
									ignore_file=args.ignore_file, g1000_alleles_dir=args.g1000_prefix_alleles, seed=args.seed, max_cn_state=args.max_cn_state, \
									use_sv_breakpoints_file=args.use_sv_breakpoints_file, gc_correction_prefix=args.gc_correction_prefix)
			
			'''
			#######################################################################
			Set up the pipelines
			#######################################################################
			'''
			if (args.type == "wgs"):
				if (args.version=="cgp"):
					runscript = generateBattenbergPipeline_CGPIT(bb_conf)
				else:
					runscript, _, _ = generateBattenbergPipeline_WGS(bb_conf, pipe_exe, pipe_rerun_exe, pipe_rerun_manual_exe, args.no_allelecount)
			elif (args.type == "snp6"):
				runscript, _, _ = generateBattenbergPipeline_SNP6(bb_conf, pipe_exe, pipe_rerun_exe, pipe_rerun_manual_exe)

			runscripts.append(runscript)
	
	# Create a master script
	scriptname = path.joinpath(args.r, "RunCommands.sh")
	runscript = open(scriptname, 'w')
	for item in runscripts:
		runscript.write(item+"\n")
	runscript.close()
	
	# Make executable
	st = os.stat(scriptname)
	os.chmod(scriptname, st.st_mode | stat.S_IEXEC)

if __name__ == '__main__':
	main(sys.argv[0:])
