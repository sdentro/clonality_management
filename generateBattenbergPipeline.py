import sys, argparse, os, stat
from path import path
# from clonalityPipelineConfig import IMPUTEFILESDIR, G1000LOCIDIR
from generateClonalityPipeline_util import read_sample_infile, generateBsubCmd, writeSimpleShellScript
from util import merge_items

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

def generateBattenbergPipeline(run_dir, log_dir, samplename, tumour_bam, normal_bam, threads):
	bb_conf = bb_pipeline_config(tumour_bam, normal_bam, run_dir, threads=threads)
	
	runscript = path.joinpath(run_dir, "RunCommands_"+samplename+".sh")
	outf = open(runscript, 'w')
	
	cmd = createAlleleCountCmd(bb_conf)
	outf.write(generateBsubCmd("loci_"+samplename, log_dir, cmd, queue="normal", mem=2, depends=None, isArray=False, threads=threads) + "\n")
	
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

def main(argv):
	parser = argparse.ArgumentParser(prog='GenerateBattenbergPipeline',
							 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("--ss", type=str, required=True, help="Full path to a samplesheet")
	parser.add_argument("-r", type=str, required=True, help="Full path to a directory where the pipelines will be ran")
	parser.add_argument("-t", type=int, required=True, help="Number of threads to use")
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
		
			runscript = generateBattenbergPipeline(run_dir, log_dir, tumourid, tb, nb, args.t)
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
