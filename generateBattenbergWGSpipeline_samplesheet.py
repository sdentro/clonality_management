import sys, argparse, os, stat
from path import path
# from clonalityPipelineConfig import IMPUTEFILESDIR, G1000LOCIDIR
from generateClonalityPipeline_util import read_sample_infile, generateBsubCmd, writeSimpleShellScript
from util import merge_items


#################################################################################################################
# CGP-IT BB setup 
#################################################################################################################

BBSCRIPT = "perl ~/repo/cgpBattenberg/perl/bin/battenberg.pl"
GENOME_INDEX = "/lustre/scratch110/sanger/sd11/Documents/GenomeFiles/refs_icgc_pancan/genome.fa.fai"
IMPUTE_INFO = "/lustre/scratch110/sanger/sd11/Documents/GenomeFiles/battenberg_impute/impute_info.txt"
G1000_LOCI_DIR = "/lustre/scratch110/sanger/sd11/Documents/GenomeFiles/battenberg_1000genomesloci2012/"
PROB_LOCI_FILE = "/lustre/scratch110/sanger/sd11/Documents/GenomeFiles/battenberg_probloci/probloci.txt"
IGNORE_FILE = "/lustre/scratch110/sanger/sd11/Documents/GenomeFiles/battenberg_ignore/ignore.txt"
PROTOCOL = "WGS"

class bb_pipeline_config(object):
	
	def __init__(self, tumour_bam, normal_bam, run_dir, genome_index=None, impute_info=None, g1000_loci_dir=None, prob_loci_file=None, ignore_file=None, protocol=None, threads=None):
		self._tumour_bam = tumour_bam
		self._normal_bam = normal_bam
		self._run_dir = run_dir
		self._genome_index = GENOME_INDEX if genome_index is None else genome_index
		self._impute_info = IMPUTE_INFO if impute_info is None else impute_info
		self._g1000_loci_dir = G1000_LOCI_DIR if g1000_loci_dir is None else g1000_loci_dir
		self._prob_loci = PROB_LOCI_FILE if prob_loci_file is None else prob_loci_file
		self._ignore_file = IGNORE_FILE if ignore_file is None else ignore_file
		self._protocol = PROTOCOL if protocol is None else protocol
		self._threads = threads if threads is not None else None 
		
	def getStandardOptions(self):
		options = ["-o", self._run_dir,
				"-r", self._genome_index,
				"-e", self._impute_info,
				"-u", self._g1000_loci_dir,
				"-c", self._prob_loci,
				"-ig", self._ignore_file,
				"-pr", self._protocol,
				"-nb", self._normal_bam,
				"-tb", self._tumour_bam]
		return(merge_items(options))
	
	def getThreadsOption(self):
		if self._threads is None:
			return("")
		else:
			return(merge_items(["-t", str(self._threads)]))

def createAlleleCountCmd(bb_conf):
	cmd = [BBSCRIPT, bb_conf.getStandardOptions(), bb_conf.getThreadsOption(),
		"-p", "allelecount"]
	return(merge_items(cmd))

def createBafLogCmd(bb_conf):
	cmd = [BBSCRIPT, bb_conf.getStandardOptions(), bb_conf.getThreadsOption(),
		"-p", "baflog"]
	return(merge_items(cmd))

def createImputeFromBafCmd(bb_conf):
	cmd = [BBSCRIPT, bb_conf.getStandardOptions(), bb_conf.getThreadsOption(),
		"-p", "imputefromaf"]
	return(merge_items(cmd))	

def createImputeCmd(bb_conf):
	cmd = [BBSCRIPT, bb_conf.getStandardOptions(), bb_conf.getThreadsOption(),
		"-p", "impute"]
	return(merge_items(cmd))	

def createCombineImputeCmd(bb_conf):
	cmd = [BBSCRIPT, bb_conf.getStandardOptions(), bb_conf.getThreadsOption(),
		"-p", "combineimpute"]
	return(merge_items(cmd))	

def createHaplotypeBafsCmd(bb_conf):
	cmd = [BBSCRIPT, bb_conf.getStandardOptions(), bb_conf.getThreadsOption(),
		"-p", "haplotypebafs"]
	return(merge_items(cmd))	

def createCleanupPostBafCmd(bb_conf):
	cmd = [BBSCRIPT, bb_conf.getStandardOptions(), bb_conf.getThreadsOption(),
		"-p", "cleanuppostbaf"]
	return(merge_items(cmd))

def createPlotHaplotypesCmd(bb_conf):
	cmd = [BBSCRIPT, bb_conf.getStandardOptions(), bb_conf.getThreadsOption(),
		"-p", "plothaplotypes"]
	return(merge_items(cmd))

def createCombineBafsCmd(bb_conf):
	cmd = [BBSCRIPT, bb_conf.getStandardOptions(),
		"-p", "combinebafs"]
	return(merge_items(cmd))

def createSegmentPhasedCmd(bb_conf):
	cmd = [BBSCRIPT, bb_conf.getStandardOptions(),
		"-p", "segmentphased"]
	return(merge_items(cmd))

def createFitcnCmd(bb_conf):
	cmd = [BBSCRIPT, bb_conf.getStandardOptions(),
		"-p", "fitcn"]
	return(merge_items(cmd))

def createSubclonesCmd(bb_conf):
	cmd = [BBSCRIPT, bb_conf.getStandardOptions(),
		"-p", "subclones"]
	return(merge_items(cmd))

def createFinaliseCmd(bb_conf):
	cmd = [BBSCRIPT, bb_conf.getStandardOptions(),
		"-p", "finalise"]
	return(merge_items(cmd))

def generateBattenbergPipelineCGPIT(run_dir, log_dir, samplename, tumour_bam, normal_bam, threads):
	bb_conf = bb_pipeline_config(tumour_bam, normal_bam, run_dir, threads=threads)
	
	runscript = path.joinpath(run_dir, "RunCommands_"+samplename+".sh")
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

#################################################################################################################
# Original setup 
#################################################################################################################

# Set some default options
PIPE_DIR='/nfs/users/nfs_s/sd11/software/pipelines/battenberg_v1.0'
IMPUTEINFOFILE='/lustre/scratch110/sanger/sd11/Documents/GenomeFiles/battenberg_impute_v3/impute_info.txt'
IMPUTE_EXE='impute2' #/lustre/scratch110/sanger/sd11/epitax/battenberg/PD7404a/impute_v2.2.2_x86_64_static/impute2'
PROBLEMLOCIFILE='/lustre/scratch110/sanger/sd11/Documents/GenomeFiles/battenberg_probloci/probloci.txt'
G1000_PREFIX="/lustre/scratch110/sanger/sd11/Documents/GenomeFiles/battenberg_1000genomesloci2012_v3/1000genomesAlleles2012_chr"
IS_MALE=False
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

def generateBattenbergConfig(tumourname, normalname, run_dir, pipeline_dir, log_dir, g1000_prefix=G1000_PREFIX, is_male=IS_MALE, platform_gamma=PLATFORM_GAMMA, phasing_gamma=PHASING_GAMMA, segmentation_gamma=SEGMENTATION_GAMMA, clonality_dist_metric=CLONALITY_DIST_METRIC, ascat_dist_metric=ASCAT_DIST_METRIC, min_ploidy=MIN_PLOIDY, max_ploidy=MAX_PLOIDY, min_rho=MIN_RHO, min_goodness_of_fit=MIN_GOODNESS_OF_FIT, balanced_threshold=BALANCED_THRESHOLD, imputeinfofile=IMPUTEINFOFILE, impute_exe=IMPUTE_EXE, problemlocifile=PROBLEMLOCIFILE):
	config_file = path.joinpath(run_dir, 'params'+tumourname+'.txt')
	f = open(config_file, 'w')
	f.write('RUN_DIR='+run_dir+'\n')
	f.write('PIPELINE_DIR='+pipeline_dir+'\n')
	f.write('LOG_DIR='+log_dir+'\n')
	f.write('TUMOURNAME='+tumourname+'\n')
	f.write('NORMALNAME='+normalname+'\n')
	f.write('PROBLEMLOCI='+problemlocifile+'\n')
	f.write('IMPUTEINFOFILE='+imputeinfofile+'\n')
	f.write('IMPUTE_EXE='+impute_exe+'\n')
	f.write('G1000_PREFIX='+g1000_prefix+'\n')
	if is_male: # == 'male' or is_male == 'Male':
		f.write('IS_MALE=TRUE\n')
	else:
		f.write('IS_MALE=FALSE\n')
	f.write('PLATFORM_GAMMA='+str(platform_gamma)+'\n')
	f.write('PHASING_GAMMA='+str(phasing_gamma)+'\n')
	f.write('SEGMENTATION_GAMMA='+str(segmentation_gamma)+'\n')
	f.write('CLONALITY_DIST_METRIC='+str(clonality_dist_metric)+'\n')
	f.write('ASCAT_DIST_METRIC='+str(ascat_dist_metric)+'\n')
	f.write('MIN_PLOIDY='+str(min_ploidy)+'\n')
	f.write('MAX_PLOIDY='+str(max_ploidy)+'\n')
	f.write('MIN_RHO='+str(min_rho)+'\n')
	f.write('MIN_GOODNESS_OF_FIT='+str(min_goodness_of_fit)+'\n')
	f.write('BALANCED_THRESHOLD='+str(balanced_threshold)+'\n')
	f.close()
	return config_file


def generateBattenbergPipeline(tumourname, normalname, run_dir, pipeline_dir, log_dir, g1000_prefix=G1000_PREFIX, is_male=IS_MALE, platform_gamma=PLATFORM_GAMMA, phasing_gamma=PHASING_GAMMA, segmentation_gamma=SEGMENTATION_GAMMA, clonality_dist_metric=CLONALITY_DIST_METRIC, ascat_dist_metric=ASCAT_DIST_METRIC, min_ploidy=MIN_PLOIDY, max_ploidy=MAX_PLOIDY, min_rho=MIN_RHO, min_goodness_of_fit=MIN_GOODNESS_OF_FIT, balanced_threshold=BALANCED_THRESHOLD, imputeinfofile=IMPUTEINFOFILE, impute_exe=IMPUTE_EXE, problemlocifile=PROBLEMLOCIFILE):

	# Write the params file
	config_file = generateBattenbergConfig(tumourname, normalname, run_dir, pipeline_dir, log_dir, g1000_prefix, is_male, platform_gamma, phasing_gamma, segmentation_gamma, clonality_dist_metric, ascat_dist_metric, min_ploidy, max_ploidy, min_rho, min_goodness_of_fit, balanced_threshold, imputeinfofile, impute_exe, problemlocifile)
	# Return run commands
	return (path.joinpath(pipeline_dir, 'RunCommands.sh')+' '+config_file), (path.joinpath(pipeline_dir, 'RunCommandsRerunFitCopynumber.sh')+' '+config_file)
	
def generateBattenbergPipelines(infile, run_dir, pipeline_dir, log_dir, g1000_prefix=G1000_PREFIX, platform_gamma=PLATFORM_GAMMA, phasing_gamma=PHASING_GAMMA, segmentation_gamma=SEGMENTATION_GAMMA, clonality_dist_metric=CLONALITY_DIST_METRIC, ascat_dist_metric=ASCAT_DIST_METRIC, min_ploidy=MIN_PLOIDY, max_ploidy=MAX_PLOIDY, min_rho=MIN_RHO, min_goodness_of_fit=MIN_GOODNESS_OF_FIT, balanced_threshold=BALANCED_THRESHOLD, imputeinfofile=IMPUTEINFOFILE, impute_exe=IMPUTE_EXE, problemlocifile=PROBLEMLOCIFILE):
	ss = read_sample_infile(infile)

	# Create a directory for each normal vs tumour taking only the first normal mentioned in the samplesheet.
	for samplename in ss.getSamplenames():
		print(samplename)
		for tumour in ss.getTumours(samplename):
			print("\t"+tumour)
			generateBattenbergPipeline(tumourname=tumour, 
									normalname=ss.getNormals(samplename)[0], 
									run_dir=run_dir, 
									pipeline_dir=pipeline_dir, 
									log_dir=log_dir,
									g1000_prefix=g1000_prefix,
									is_male=ss.isMale(samplename), 
									platform_gamma=platform_gamma, 
									phasing_gamma=phasing_gamma, 
									segmentation_gamma=segmentation_gamma, 
									clonality_dist_metric=clonality_dist_metric, 
									ascat_dist_metric=ascat_dist_metric, 
									min_ploidy=min_ploidy, 
									max_ploidy=max_ploidy, 
									min_rho=min_rho, 
									min_goodness_of_fit=min_goodness_of_fit, 
									balanced_threshold=balanced_threshold, 
									imputeinfofile=imputeinfofile, 
									impute_exe=impute_exe, 
									problemlocifile=problemlocifile)

def main(argv):
	parser = argparse.ArgumentParser(prog='GenerateBattenbergPipeline',
							 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("--ss", type=str, required=True, help="Full path to a samplesheet")
	parser.add_argument("-r", type=str, required=True, help="Full path to a directory where the pipelines will be ran")
	
	parser.add_argument("--cgpit", action='store_true')
	parser.add_argument("-t", type=int, help="Number of threads to use")
	
	parser.set_defaults(cgpit=False, t=1)
	args = parser.parse_args()

	runscripts = []
	# Read in samplesheet
	ss = read_sample_infile(args.ss)
	
	# For every entry:
	for sample in ss.getSamplenames():
		print(sample)
		
		tn_pair = ss.getTumour2NormalPairingBam(sample)
		for tb,nb in tn_pair:
			tumourid = ss.getIdByTumourBam(tb)
	
			run_dir = path.joinpath(args.r, tumourid)
			log_dir = path.joinpath(args.r, tumourid, "logs")
			if not run_dir.exists(): run_dir.makedirs()
			if not log_dir.exists(): log_dir.makedirs()
		
			if (args.cgpit):
				runscript = generateBattenbergPipelineCGPIT(run_dir, log_dir, tumourid, tb, nb, args.t)
			else:
				runscript, _ = generateBattenbergPipeline(tumourname=tumourid, 
									normalname=ss.getNormals(sample)[0], 
									run_dir=run_dir, 
									pipeline_dir=PIPE_DIR, 
									log_dir=log_dir,
									is_male=ss.isMale(sample),
									g1000_prefix=G1000_PREFIX, 
									platform_gamma=PLATFORM_GAMMA, 
									phasing_gamma=PHASING_GAMMA, 
									segmentation_gamma=SEGMENTATION_GAMMA, 
									clonality_dist_metric=CLONALITY_DIST_METRIC, 
									ascat_dist_metric=ASCAT_DIST_METRIC, 
									min_ploidy=MIN_PLOIDY, 
									max_ploidy=MAX_PLOIDY, 
									min_rho=MIN_RHO, 
									min_goodness_of_fit=MIN_GOODNESS_OF_FIT, 
									balanced_threshold=BALANCED_THRESHOLD, 
									imputeinfofile=IMPUTEINFOFILE, 
									impute_exe=IMPUTE_EXE, 
									problemlocifile=PROBLEMLOCIFILE)
				

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
