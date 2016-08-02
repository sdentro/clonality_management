#!/usr/bin/env python

import argparse, sys
import numpy as np
from path import path
from generateClonalityPipeline_util import read_sample_infile
from clonalityPipelineConfig import RHO_AND_PSI_REGEX

def generateDPDataFile(proj_name, infile, dp_in_dir, run_dir):
	ss = read_sample_infile(infile)
	
	# Collect the purity estimate for each tumour from BB output
	tumour2purity = getTumourPurity(ss)
	
	# Create an inventory of all available dp input files
	dp_in_files = np.array(path(dp_in_dir).listdir("*allDirichletProcessInfo.txt"))
	cn_dp_in_files = np.array(path(dp_in_dir).listdir("*cnDirichletInput.txt"))
	print(cn_dp_in_files)
	
	outfile = open(path.joinpath(run_dir, proj_name+".txt"), 'w')	
	outfile.write("sample\tsubsample\tdatafile\tcellularity\tsex\tcndatafile\n")
	for sample in ss.getSamplenames():
		lines = []
		for tumour in ss.getTumours(sample):
			dp_in_file = dp_in_files[np.array([tumour in item for item in dp_in_files])]
			if not len(dp_in_file) == 1:
				print(dp_in_file)
				print("Found different than expected dp input matches for "+tumour)
				continue
			
			# Adding in CNA input files
			if len(cn_dp_in_files) == 0:
				cn_dp_in_file = np.array(["NA"])
			else:
				cn_dp_in_file = cn_dp_in_files[np.array([tumour in item for item in cn_dp_in_files])]
			if len(cn_dp_in_file) > 1:
				print("Found more than one CNA dp input matches for "+tumour)
				continue
			elif len(cn_dp_in_file) == 0:
				cn_dp_in_file = np.array(["NA"])

			if ss.getSex(tumour) is None:
				sex = "male"
			else:
				sex = ss.getSex(tumour)

			if tumour in tumour2purity.keys():
				lines = lines + [sample+"\t"+tumour+"\t"+path(dp_in_file[0]).basename()+"\t"+tumour2purity[tumour]+"\t"+sex+"\t"+path(cn_dp_in_file[0]).basename()+"\n"]
				#outfile.write(sample+"\t"+tumour+"\t"+path(dp_in_file[0]).basename()+"\t"+tumour2purity[tumour]+"\t"+sex+"\t"+path(cn_dp_in_file[0]).basename()+"\n")
			else:
				print("Did not find purity estimate for "+tumour)

		if len(lines) == len(ss.getTumours(sample)):
			for line in lines: outfile.write(line)
		else:
			print("Incomplete multi sample case "+sample)
	outfile.close()
		
def getTumourPurity(ss):
	'''
	Fetches the rho_and_psi file for each tumour_id, for each sample
	'''
	sample2purity = dict()
	
	# For every sample
	for sample in ss.getSamplenames():
		# For every tumour_id available for that sample
		for tumour_id in ss.getTumours(sample):
			bb_dir = ss.getBbDirByTumourId(sample, tumour_id)
			if path(bb_dir).exists():
				listing = np.array(path(bb_dir).listdir(tumour_id+RHO_AND_PSI_REGEX))
				
				if (len(listing) == 1):
					rho_and_psi = listing[np.array([tumour_id in filename for filename in listing])] # Should yield only a single file
				else:
					print("Found different than expected number of matches for "+tumour_id)
					continue
				
				if not (len(rho_and_psi) == 1):
					print("Found different than expected number of matches for "+tumour_id)
					continue
				if path(rho_and_psi[0]).exists():
					sample2purity[tumour_id] = getPurityPerSample(rho_and_psi[0])
				else:
					print("Did not find rho_and_psi file for "+tumour_id)
					continue
			else:
				print("Did not find BB directory for "+tumour_id)
				continue
	
	return sample2purity
	
def getPurityPerSample(rho_and_psi_file):
	f = open(rho_and_psi_file,'r')
	f.readline()
	f.readline()
	frac_genome_line = f.readline().split("\t")
	f.close()
	if(frac_genome_line[0] != 'FRAC_GENOME'):
		print("Unexpected rho_and_psi format in file "+rho_and_psi_file)
		return("")
	else:
		return(frac_genome_line[1])

def main(argv):
	parser = argparse.ArgumentParser(prog='GenerateDPDataFile',
							 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("-i", required=True, type=str, help="Full path to a sample sheet")
	parser.add_argument("-d", required=True, type=str, help="Full path to where Dirichlet input files are stored")
	parser.add_argument("-r", required=True, type=str, help="Full path to directory where the pipeline is going to be run.")
	parser.add_argument("-p", required=True, type=str, help="Name of the project.")
#	 parser.add_argument("--no_iters", type=int, help="Number of iterations the 1D DP has to run.")
#	 parser.add_argument("--no_iters_burn_in", type=int, help="Number of iterations used for burn in.")
	
	args = parser.parse_args()
	generateDPDataFile(args.p, args.i, args.d, args.r)
	
if __name__ == '__main__':
	main(sys.argv[0:])
