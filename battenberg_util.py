import sys
from util import merge_items

# WGS cgpIT and SNP6
class bb_pipeline_config(object):
    
    def __init__(self, pipe_type, pipe_version, \
                 tumour_file, normal_file, run_dir, log_dir, samplename, tumour_id=None, normal_id=None, \
                 gender=None, genome_index=None, impute_info=None, g1000_loci_dir=None, g1000_alleles_dir=None, \
                 prob_loci_file=None, ignore_file=None, protocol=None, threads=None, pipe_dir=None, \
                 platform_gamma=None, phasing_gamma=None, segmentation_gamma=None, clonality_dist_metric=None, \
                 ascat_dist_metric=None, min_ploidy=None, max_ploidy=None, min_rho=None, min_goodness_of_fit=None, \
                 balanced_threshold=None, impute_exe=None, min_count=None, heterozygous_filter=None, \
                 fill_in_snps=None, use_tumour_snps=None, use_het_snps_only=None, hom_caveman_file=None, \
                 use_loci_file=None, snppos_file=None, gc_snp6_file=None, snp6_anno_file=None, \
                 snp6_ref_info_file=None, birdseed_report_file=None, apt_probeset_geno_exe=None, \
                 apt_probeset_summ_exe=None, norm_geno_clust_exe=None):
        
        # General
        self._pipe_type = pipe_type
        self._pipe_version = pipe_version
        
        self._tumour_file = tumour_file
        self._normal_file = normal_file
        self._tumour_id = tumour_id
        self._normal_id = normal_id
        self._gender = gender
        self._samplename = samplename
        
        self._run_dir = run_dir
        self._log_dir = log_dir
        self._pipe_dir = pipe_dir
        
        self._impute_exe = impute_exe
        self._impute_info = impute_info
        self._g1000_loci_dir = g1000_loci_dir
	self._g1000_alleles_dir = g1000_alleles_dir
        self._prob_loci = prob_loci_file
        
        self._platform_gamma = platform_gamma
        self._phasing_gamma = phasing_gamma
        self._segmentation_gamma = segmentation_gamma
        self._clonality_dist_metric = clonality_dist_metric
        self._ascat_dist_metric = ascat_dist_metric
        
        self._min_ploidy = min_ploidy
        self._max_ploidy = max_ploidy
        self._min_rho = min_rho
        self._min_goodness_of_fit = min_goodness_of_fit
        self._balanced_threshold = balanced_threshold
        
        self._min_count = min_count
        
        # SNP6
        self._heterozygous_filter = heterozygous_filter
        self._fill_in_snps = fill_in_snps
        self._use_tumour_snps = use_tumour_snps
        self._use_het_snps_only = use_het_snps_only
        self._hom_caveman_file = hom_caveman_file
        self._use_loci_file = use_loci_file
        self._snppos_file = snppos_file
        self._gc_snp6_file = gc_snp6_file
        self._snp6_anno_file = snp6_anno_file
        self._snp6_ref_info_file = snp6_ref_info_file
        self._birdseed_report_file = birdseed_report_file
        self._apt_probeset_geno_exe = apt_probeset_geno_exe
        self._apt_probeset_summ_exe = apt_probeset_summ_exe
        self._norm_geno_clust_exe = norm_geno_clust_exe
        
        # cgpBB
        self._genome_index = genome_index
        self._ignore_file = ignore_file
        self._protocol = protocol
        self._threads = threads 
    
    '''
    ####################################################
    Simple getters
    ####################################################
    '''
    def get_samplename(self):
        return self._samplename    
    
    def get_run_dir(self):
        return self._run_dir
    
    def get_log_dir(self):
        return self._log_dir
    
    def get_pipe_dir(self):
        return self._pipe_dir
    
    def get_threads(self):
        return self._threads
        
    '''
    ####################################################
    GCP IT convenience functions
    ####################################################
    '''
    def getStandardOptions_cgpBB(self):
        '''
        Returns the default options that always need to be supplied when calling cgpBB
        '''
        options = ["-o", self._run_dir,
                "-r", self._genome_index,
                "-e", self._impute_info,
                "-u", self._g1000_loci_dir,
                "-c", self._prob_loci,
                "-ig", self._ignore_file,
                "-pr", self._protocol,
                "-nb", self._normal_file,
                "-tb", self._tumour_file]
        return(merge_items(options))
    
    
    def getThreadsOption_cgpBB(self):
        '''
        Returns the number of threads parameter as command line option for cgpBB
        '''
        if self._threads is None:
            return("")
        else:
            return(merge_items(["-t", str(self._threads)]))
        
    '''
    ####################################################
    Params files generators
    ####################################################
    '''
    def _generateSharedParams_sample(self, fout):
        fout.write("RUN_DIR="+self._run_dir+"\n")
        fout.write("LOG_DIR="+self._log_dir+"\n")
        fout.write("PIPELINE_DIR="+self._pipe_dir+"\n")
        fout.write("TUMOURNAME="+self._tumour_id+"\n")
        if self._gender == "male" or self._gender == "Male":
            is_male = "TRUE"
        elif self._gender == "Female" or self._gender == "female":
            is_male = "FALSE"
        else:
            print("Supplied gender for sample "+self._samplename+" not male or female")
            sys.exit(1)
        fout.write("IS_MALE="+is_male+"\n")
            
    def _generateSharedParams_general(self, fout):
        fout.write("PLATFORM_GAMMA="+str(self._platform_gamma)+"\n")
        fout.write("PHASING_GAMMA="+str(self._phasing_gamma)+"\n")
        fout.write("SEGMENTATION_GAMMA="+str(self._segmentation_gamma)+"\n")
        fout.write("CLONALITY_DIST_METRIC="+str(self._clonality_dist_metric)+"\n")
        fout.write("ASCAT_DIST_METRIC="+str(self._ascat_dist_metric)+"\n")
        fout.write("MIN_PLOIDY="+str(self._min_ploidy)+"\n")
        fout.write("MAX_PLOIDY="+str(self._max_ploidy)+"\n")
        fout.write("MIN_RHO="+str(self._min_rho)+"\n")
        fout.write("MIN_GOODNESS_OF_FIT="+str(self._min_goodness_of_fit)+"\n")
        fout.write("BALANCED_THRESHOLD="+str(self._balanced_threshold)+"\n")
        fout.write("IMPUTEINFOFILE="+self._impute_info+"\n")
        fout.write("IMPUTE_EXE="+self._impute_exe+"\n")
        
        
    def generateParamsFile_WGS(self, outfile):
        '''
        Creates the parameters file that BB_WGS uses and writes them to the file supplied as outfile
        '''
        fout = open(outfile, 'w')
        self._generateSharedParams_sample(fout)
        fout.write("NORMALNAME="+self._normal_id+"\n")
        fout.write("TUMOURBAM="+self._tumour_file+"\n")
        fout.write("NORMALBAM="+self._normal_file+"\n")
    
        self._generateSharedParams_general(fout)
        
        fout.write("PROBLEMLOCI="+self._prob_loci+"\n")
        fout.write("G1000_PREFIX="+self._g1000_alleles_dir+"\n")
        fout.write("G1000_PREFIX_AC="+self._g1000_loci_dir+"\n") # TODO: This option should be removed??
	fout.write("MIN_NORMAL_DEPTH="+str(self._min_count)+"\n")
    
        fout.close()
        
    def generateParamsFile_SNP6(self, outfile):
        '''
        Creates the parameters file that BB_snp6 uses
        '''
        fout = open(outfile, 'w')
        self._generateSharedParams_sample(fout)
        fout.write("NORMALCEL="+self._normal_file+"\n")
        fout.write("TUMOURCEL="+self._tumour_file+"\n")
    
        self._generateSharedParams_general(fout)
        fout.write("HETEROZYGOUS_FILTER="+self._heterozygous_filter+"\n")
        fout.write("MIN_COUNT="+str(self._min_count)+"\n")
        fout.write("FILL_IN_SNPS="+str(self._fill_in_snps)+"\n") # TODO: used?
        fout.write("USE_TUMOUR_SNPS="+str(self._use_tumour_snps)+"\n") # TODO: used?
        fout.write("USE_HETEROZYGOUS_SNPS_ONLY="+str(self._use_het_snps_only)+"\n") # TODO: used?
        fout.write("HOMOZYGOUS_CAVEMAN_CALLS_FILE="+self._hom_caveman_file+"\n") # TODO: used?
        fout.write("PROB_LOCI_FILE="+self._prob_loci+"\n") # TODO: merge with PROB_LOCI ?
        fout.write("USE_LOCI_FILE="+self._use_loci_file+"\n") # TODO: used?
        
        fout.write("SNPPOS="+self._snppos_file+"\n")
        fout.write("GC_SNP6="+self._gc_snp6_file+"\n")
        fout.write("ANNO_FILE="+self._snp6_anno_file+"\n")
        fout.write("SNP6_REF_INFO_FILE="+self._snp6_ref_info_file+"\n")
        fout.write("BIRDSEED_REPORT_FILE="+self._birdseed_report_file+"\n")
        fout.write("APT_PROBESET_GENOTYPE_EXE="+self._apt_probeset_geno_exe+"\n")
        fout.write("APT_PROBESET_SUMMARIZE_EXE="+self._apt_probeset_summ_exe+"\n")
        fout.write("NORM_GENO_CLUST_EXE="+self._norm_geno_clust_exe+"\n")
        
        fout.close()
