import sys, argparse
from path import path
from clonalityPipelineConfig import IMPUTEFILESDIR, G1000LOCIDIR
from generateClonalityPipeline_util import read_sample_infile
'''
RUN_DIR=/lustre/scratch110/sanger/sd11/epitax/battenberg/PD7404a
PIPELINE_DIR=/lustre/scratch110/sanger/sd11/epitax/battenberg/PD7404a
TUMOURNAME=PD7404a
NORMALNAME=PD7404b
'''

# Set some default options
IMPUTEINFOFILE='impute_info.txt' #/lustre/scratch110/sanger/sd11/epitax/battenberg/PD7404a/impute_info.txt'
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


def generateBattenbergPipeline(tumourname, normalname, run_dir, pipeline_dir, log_dir, rewrite_params=False, is_male=IS_MALE, platform_gamma=PLATFORM_GAMMA, phasing_gamma=PHASING_GAMMA, segmentation_gamma=SEGMENTATION_GAMMA, clonality_dist_metric=CLONALITY_DIST_METRIC, ascat_dist_metric=ASCAT_DIST_METRIC, min_ploidy=MIN_PLOIDY, max_ploidy=MAX_PLOIDY, min_rho=MIN_RHO, min_goodness_of_fit=MIN_GOODNESS_OF_FIT, balanced_threshold=BALANCED_THRESHOLD, imputeinfofile=IMPUTEINFOFILE, impute_exe=IMPUTE_EXE, problemlocifile=PROBLEMLOCIFILE):
	# Create sample dir
	run_dir_sample = path.joinpath(run_dir,normalname+"_vs_"+tumourname)
	
	if log_dir is None:
		log_dir = path.joinpath(run_dir_sample, 'logs')
	
	if not rewrite_params: # we're setting up a new pipeline, thus create dirs, symlinks and the lot
		run_dir_sample.makedirs()
		# Make a symlink to the impute and 1000 genomes files that are needed
		path.joinpath(IMPUTEFILESDIR).symlink(path.joinpath(run_dir_sample, 'impute'))
		path.joinpath(G1000LOCIDIR).symlink(path.joinpath(run_dir_sample, '1000genomesloci'))
		
		# Create logs directory
		if log_dir is None:
			log_dir.makedirs()

	# write the params file
	config_file = generateBattenbergConfig(tumourname, normalname, run_dir_sample, pipeline_dir, log_dir, is_male, platform_gamma, phasing_gamma, segmentation_gamma, clonality_dist_metric, ascat_dist_metric, min_ploidy, max_ploidy, min_rho, min_goodness_of_fit, balanced_threshold, imputeinfofile, impute_exe, problemlocifile)
	
	if not rewrite_params: # We're creating new BB pipelines, thus we need to write the RunCommands scripts
		# write the RunCommands.sh with the only required parameter
		f = open(path.joinpath(run_dir_sample,'RunCommands.sh'),'w')
		f.write(path.joinpath(pipeline_dir, 'RunCommands.sh')+' '+config_file+'\n')
		f.close()
		
		# Create a script to rerun from fitcopynumber onwards, after a change in parameters
		f = open(path.joinpath(run_dir_sample,'RunCommandsRerunFitCopynumber.sh'),'w')
		f.write(path.joinpath(pipeline_dir, 'RunCommandsRerunFitCopynumber.sh')+' '+config_file+'\n')
		f.close()
	
def generateBattenbergPipelines(infile, run_dir, pipeline_dir, log_dir, rewrite_params=False, platform_gamma=PLATFORM_GAMMA, phasing_gamma=PHASING_GAMMA, segmentation_gamma=SEGMENTATION_GAMMA, clonality_dist_metric=CLONALITY_DIST_METRIC, ascat_dist_metric=ASCAT_DIST_METRIC, min_ploidy=MIN_PLOIDY, max_ploidy=MAX_PLOIDY, min_rho=MIN_RHO, min_goodness_of_fit=MIN_GOODNESS_OF_FIT, balanced_threshold=BALANCED_THRESHOLD, imputeinfofile=IMPUTEINFOFILE, impute_exe=IMPUTE_EXE, problemlocifile=PROBLEMLOCIFILE):
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
									rewrite_params=rewrite_params, 
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
	parser.add_argument("-i", required=True, type=str, help="Full path to samplesheet. When multiple normals are defined only the first normal will be used.")
	parser.add_argument("-r", required=True, type=str, help="Full path to directory where the pipeline is going to be run. This directory cannot exist yet.")
	parser.add_argument("-p", required=True, type=str, help="Full path to directory where the Battenberg pipeline is stored.")

	parser.add_argument("--rewrite_params", action="store_true", help="Supply this boolean when the paramaters input files should be rewritten only into existing pipelines")
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
	parser.set_defaults(rewrite_params=False, platform_gamma=PLATFORM_GAMMA, phasing_gamma=PHASING_GAMMA, segmentation_gamma=SEGMENTATION_GAMMA, clonality_dist_metric=CLONALITY_DIST_METRIC, ascat_dist_metric=ASCAT_DIST_METRIC, min_ploidy=MIN_PLOIDY, max_ploidy=MAX_PLOIDY, min_rho=MIN_RHO, min_goodness_of_fit=MIN_GOODNESS_OF_FIT, balanced_threshold=BALANCED_THRESHOLD)

	args = parser.parse_args()
	# Default will be used, which is battenberg_dir/samplename/logs
	log_dir = None if args.log_dir is None else args.log_dir

	generateBattenbergPipelines(infile=args.i, 
							run_dir=args.r, 
							pipeline_dir=args.p, 
							log_dir=log_dir, 
							rewrite_params=args.rewrite_params, 
							platform_gamma=args.platform_gamma, 
							phasing_gamma=args.phasing_gamma, 
							segmentation_gamma=args.segmentation_gamma, 
							clonality_dist_metric=args.clonality_dist_metric, 
							ascat_dist_metric=args.ascat_dist_metric, 
							min_ploidy=args.min_ploidy, 
							max_ploidy=args.max_ploidy, 
							min_rho=args.min_rho, 
							min_goodness_of_fit=args.min_goodness_of_fit, 
							balanced_threshold=args.balanced_threshold, 
							imputeinfofile=path.joinpath(args.p, IMPUTEINFOFILE), 
							impute_exe=path.joinpath(args.p, IMPUTE_EXE), 
							problemlocifile=path.joinpath(args.p, PROBLEMLOCIFILE))
	
	print("")
	print("Don't forget to copy the GetAlleleFrequencies output into the sample battenberg directories")
	print("")
	

if __name__ == '__main__':
	main(sys.argv[0:])
