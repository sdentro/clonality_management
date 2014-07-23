import sys, argparse
from path import path
from clonalityPipelineConfig import IMPUTEFILESDIR, G1000LOCIDIR
from generateClonalityPipeline_util import read_item_list, match_tumour_normal
'''
RUN_DIR=/lustre/scratch110/sanger/sd11/epitax/battenberg/PD7404a
PIPELINE_DIR=/lustre/scratch110/sanger/sd11/epitax/battenberg/PD7404a
TUMOURNAME=PD7404a
NORMALNAME=PD7404b
'''

# Set some default options
IMPUTEINFOFILE='impute_info.txt' #/lustre/scratch110/sanger/sd11/epitax/battenberg/PD7404a/impute_info.txt'i
IMPUTE_EXE='impute_v2.2.2_x86_64_static/impute2' #/lustre/scratch110/sanger/sd11/epitax/battenberg/PD7404a/impute_v2.2.2_x86_64_static/impute2'
PROBLEMLOCIFILE='probloci.txt'
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

def generateBattenbergConfig(tumourname, normalname, run_dir, pipeline_dir, log_dir, is_male=IS_MALE, platform_gamma=PLATFORM_GAMMA, phasing_gamma=PHASING_GAMMA, segmentation_gamma=SEGMENTATION_GAMMA, clonality_dist_metric=CLONALITY_DIST_METRIC, ascat_dist_metric=ASCAT_DIST_METRIC, min_ploidy=MIN_PLOIDY, max_ploidy=MAX_PLOIDY, min_rho=MIN_RHO, min_goodness_of_fit=MIN_GOODNESS_OF_FIT, balanced_threshold=BALANCED_THRESHOLD, imputeinfofile=IMPUTEINFOFILE, impute_exe=IMPUTE_EXE, problemlocifile=PROBLEMLOCIFILE):
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
	if is_male == 'male' or is_male == 'Male':
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


def generateBattenbergPipeline(tumourname, normalname, run_dir, pipeline_dir, log_dir, is_male=IS_MALE, platform_gamma=PLATFORM_GAMMA, phasing_gamma=PHASING_GAMMA, segmentation_gamma=SEGMENTATION_GAMMA, clonality_dist_metric=CLONALITY_DIST_METRIC, ascat_dist_metric=ASCAT_DIST_METRIC, min_ploidy=MIN_PLOIDY, max_ploidy=MAX_PLOIDY, min_rho=MIN_RHO, min_goodness_of_fit=MIN_GOODNESS_OF_FIT, balanced_threshold=BALANCED_THRESHOLD, imputeinfofile=IMPUTEINFOFILE, impute_exe=IMPUTE_EXE, problemlocifile=PROBLEMLOCIFILE):
	# Create sample dir
	run_dir_sample = path.joinpath(run_dir,tumourname)
	run_dir_sample.makedirs()
	# Make a symlink to the impute and 1000 genomes files that are needed
	path.joinpath(IMPUTEFILESDIR).symlink(path.joinpath(run_dir_sample, 'impute'))
	path.joinpath(G1000LOCIDIR).symlink(path.joinpath(run_dir_sample, '1000genomesloci'))
	
	# Create logs directory
	if log_dir is None:
		log_dir = path.joinpath(run_dir_sample, 'logs')
		log_dir.makedirs()

	# write the params file
	config_file = generateBattenbergConfig(tumourname, normalname, run_dir_sample, pipeline_dir, log_dir, is_male, platform_gamma, phasing_gamma, segmentation_gamma, clonality_dist_metric, ascat_dist_metric, min_ploidy, max_ploidy, min_rho, min_goodness_of_fit, balanced_threshold, imputeinfofile, impute_exe, problemlocifile)
	
	# write the RunCommands.sh with the only required parameter
	f = open(path.joinpath(run_dir_sample,'RunCommands.sh'),'w')
	f.write(path.joinpath(pipeline_dir, 'RunCommands.sh')+' '+config_file+'\n')
	f.close()
	

def generateBattenbergPipelines(infile, tumour_id, normal_id, run_dir, pipeline_dir, log_dir, is_male=IS_MALE, platform_gamma=PLATFORM_GAMMA, phasing_gamma=PHASING_GAMMA, segmentation_gamma=SEGMENTATION_GAMMA, clonality_dist_metric=CLONALITY_DIST_METRIC, ascat_dist_metric=ASCAT_DIST_METRIC, min_ploidy=MIN_PLOIDY, max_ploidy=MAX_PLOIDY, min_rho=MIN_RHO, min_goodness_of_fit=MIN_GOODNESS_OF_FIT, balanced_threshold=BALANCED_THRESHOLD, imputeinfofile=IMPUTEINFOFILE, impute_exe=IMPUTE_EXE, problemlocifile=PROBLEMLOCIFILE):
	samplenames = read_item_list(infile)
	is_male = read_item_list(is_male)
	# Create lookup table for each sample
	sample_is_male = dict(zip(samplenames, is_male))
	# Match tumours and normals together
	matched = match_tumour_normal(samplenames, tumour_id, normal_id)
	for sample in matched:
		generateBattenbergPipeline(sample[0], sample[1], run_dir, pipeline_dir, log_dir, sample_is_male[sample[0]], platform_gamma, phasing_gamma, segmentation_gamma, clonality_dist_metric, ascat_dist_metric, min_ploidy, max_ploidy, min_rho, min_goodness_of_fit, balanced_threshold, imputeinfofile, impute_exe, problemlocifile)


def main(argv):
	parser = argparse.ArgumentParser(prog='GenerateBattenbergPipeline',
							 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("-i", required=True, type=str, help="Full path to input file containing a list of samplenames")
	parser.add_argument("-t", required=True, type=str, help="Identifier of the tumour sample")
	parser.add_argument("-n", required=True, type=str, help="Identifier of the normal sample")
	parser.add_argument("-r", required=True, type=str, help="Full path to directory where the pipeline is going to be run. This directory cannot exist yet.")
	parser.add_argument("-p", required=True, type=str, help="Full path to directory where the Battenberg pipeline is stored.")
	parser.add_argument("--ismale", required=True, type=str, help='Path to file containing whether each sample is male or female.')
	
	parser.add_argument("--platform_gamma", type=float, help='The platform gamma')
	parser.add_argument("--phasing_gamma", type=float, help='The phasing gamma')
	parser.add_argument("--segmentation_gamma", type=float, help='The segmentation gamma')
	parser.add_argument("--clonality_dist_metric", type=float, help='The clonality dist metric')
	parser.add_argument("--ascat_dist_metric", type=float, help='The ASCAT dist metric')
	parser.add_argument("--min_ploidy", type=float, help='Minimum allowed ploidy')
	parser.add_argument("--max_ploidy", type=float, help='Maximum allowed ploidy')
	parser.add_argument("--min_rho", type=float, help='Minimum allowed rho')
	parser.add_argument("--min_goodness_of_fit", type=float, help='Minimum allowed goodness of fit')
	parser.add_argument("--balanced_threshold", type=float, help='Balanced threshold')
	parser.add_argument("--log_dir", type=str, help='Location where logfiles should be saved')
	parser.set_defaults(ismale=IS_MALE, platform_gamma=PLATFORM_GAMMA, phasing_gamma=PHASING_GAMMA, segmentation_gamma=SEGMENTATION_GAMMA, clonality_dist_metric=CLONALITY_DIST_METRIC, ascat_dist_metric=ASCAT_DIST_METRIC, min_ploidy=MIN_PLOIDY, max_ploidy=MAX_PLOIDY, min_rho=MIN_RHO, min_goodness_of_fit=MIN_GOODNESS_OF_FIT, balanced_threshold=BALANCED_THRESHOLD)

	args = parser.parse_args()
	# Default will be used, which is battenberg_dir/samplename/logs
	log_dir = None if args.log_dir is None else args.log_dir

	generateBattenbergPipelines(args.i, args.t, args.n, args.r, args.p, log_dir, args.ismale, args.platform_gamma, args.phasing_gamma, args.segmentation_gamma, args.clonality_dist_metric, args.ascat_dist_metric, args.min_ploidy, args.max_ploidy, args.min_rho, args.min_goodness_of_fit, args.balanced_threshold, path.joinpath(args.p, IMPUTEINFOFILE), path.joinpath(args.p, IMPUTE_EXE), path.joinpath(args.p, PROBLEMLOCIFILE))
	
	print("")
	print("Don't forget to copy the GetAlleleFrequencies output into the sample battenberg directories")
	print("")
	

if __name__ == '__main__':
	main(sys.argv[0:])