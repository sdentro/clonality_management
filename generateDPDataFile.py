
import argparse, sys
import numpy as np
from path import path
from generateClonalityPipeline_util import read_sample_infile
from clonalityPipelineConfig import RHO_AND_PSI_REGEX

def generateDPDataFile(proj_name, infile, bb_dir, dp_in_dir, run_dir):
    ss = read_sample_infile(infile)
    
    # Collect the purity estimate for each tumour from BB output
    tumour2purity = getTumourPurity(ss, bb_dir, run_dir)
    
    # Create an inventory of all available dp input files
    dp_in_files = np.array(path.listdir(dp_in_dir, "*.txt"))
    
    outfile = open(path.joinpath(run_dir, proj_name+".txt"))    
    outfile.write("sample\tsubsample\tdatafile\tcellularity")
    for sample in ss.getSamplenames():
        for tumour in ss.getTumours(sample):
            dp_in_file = dp_in_files[np.array([tumour in item for item in dp_in_files])]
            outfile.write(sample+"\t"+tumour+"\t"+dp_in_file+"\t"+tumour2purity[tumour])
    outfile.close()
        

def getTumourPurity(ss, bb_dir, run_dir):
    ''' 
    Create an inventory of tumour and purity and return it as dict
    '''
    sample2purity = dict()
    # TODO: Perhaps here determine whether rho_and_psi files are in subdirs?
    bb_dirs = np.array(path.joinpath(bb_dir).listdir())
    for sample in ss.getSamplenames():
        
        normal = ss.getNormals(sample)[0]
        
        for tumour in ss.getTumours(sample): # bb is run as first normal against all tumours
            sample_dir = normal+"_vs_"+tumour
        
            if len(bb_dirs[np.array([sample_dir in item for item in bb_dirs])]) > 1:
                print("Found more than one Battenberg dir for sample "+sample_dir)
                sys.exit(1)
                
            indir = bb_dirs[np.array([sample_dir in item for item in bb_dirs])][0] # Should yield only a single dir
            rho_and_psi = path(indir).listdir(RHO_AND_PSI_REGEX)[0] # Should yield only a single file
            sample2purity[tumour] = getPurityPerSample(rho_and_psi)
        
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
    parser.add_argument("-b", required=True, type=str, help="Full path to where Battenberg output is stored")
    parser.add_argument("-d", required=True, type=str, help="Full path to where Dirichlet input files are stored")
    parser.add_argument("-r", required=True, type=str, help="Full path to directory where the pipeline is going to be run.")
    parser.add_argument("-p", required=True, type=str, help="Name of the project.")
#     parser.add_argument("--no_iters", type=int, help="Number of iterations the 1D DP has to run.")
#     parser.add_argument("--no_iters_burn_in", type=int, help="Number of iterations used for burn in.")
    
    args = parser.parse_args()
    generateDPDataFile(args.i, args.b, args.r, args.d)
    
if __name__ == '__main__':
    main(sys.argv[0:])