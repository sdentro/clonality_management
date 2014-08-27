import sys, argparse
from path import path
import numpy as np
from generateClonalityPipeline_util import read_item_list
from clonalityPipelineConfig import RHO_AND_PSI_REGEX, DIRICHLET_1D_RUNSCRIPT


def generateDirichlet1D(infile, bb_dir, run_dir, dirichlet_indir, no_iters, no_iters_burn_in):
    samplenames = read_item_list(infile)
    
    # make list of purities per sample
    purity_file = generateSample2PurityFile(samplenames, bb_dir, run_dir)  
    
    # create run script using list of samples, list of purities and run number per sample
    counter = 1
    fout = open(path.joinpath(run_dir, "RunCommands.sh"), 'w')
    for sample in samplenames:
        fout.write("bsub -q long -M 5000 -R 'span[hosts=1] select[mem>5000] rusage[mem=5000]' -o "+
                   path.joinpath(run_dir, 'logs', sample+".%J.out")+" -e "+
                   path.joinpath(run_dir, 'logs', sample+".%J.err ")+
                   "-J "+sample+"_1D "+
                   "Rscript "+DIRICHLET_1D_RUNSCRIPT+" "+
                   str(counter)+" "+str(no_iters)+" "+str(no_iters_burn_in)+" "+dirichlet_indir+" "+purity_file+"\n")
	counter+=1
    fout.close()

    
def generateSample2PurityFile(samplenames, bb_dir, run_dir):
    # Create an inventory of samplename and purity
    sample2purity = dict()
    # TODO: Perhaps here determine whether rho_and_psi files are in subdirs?
    bb_dirs = np.array(path.joinpath(bb_dir).listdir())
    for sample in samplenames:
        if len(bb_dirs[np.array([sample in item for item in bb_dirs])]) > 1:
            print("Found more than one Battenberg dir for sample "+sample)
        indir = bb_dirs[np.array([sample in item for item in bb_dirs])][0] # Should yield only a single dir
        rho_and_psi = path(indir).listdir(RHO_AND_PSI_REGEX)[0] # Should yield only a single file
        sample2purity[sample] = getPurityPerSample(rho_and_psi)
        
    # Write the inventory to a file
    fout_name = path.joinpath(run_dir, "samplepurities.txt")
    fout = open(fout_name, 'w')
    fout.write("sample\tpurity\n")
    for k,v in sample2purity.iteritems():
        fout.write(k+"\t"+v+"\n")
    fout.close()
    
    return fout_name
    
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
    parser = argparse.ArgumentParser(prog='GenerateDirichlet1D',
                             formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", required=True, type=str, help="Full path to input file containing a list of samplenames")
    parser.add_argument("-b", required=True, type=str, help="Full path to where Battenberg output is stored")
    parser.add_argument("-d", required=True, type=str, help="Full path to where Dirichlet input is stored")
    parser.add_argument("-r", required=True, type=str, help="Full path to directory where the pipeline is going to be run.")
    parser.add_argument("--no_iters", type=int, help="Number of iterations the 1D DP has to run.")
    parser.add_argument("--no_iters_burn_in", type=int, help="Number of iterations used for burn in.")
    
    args = parser.parse_args()
    generateDirichlet1D(args.i, args.b, args.r, args.d, args.no_iters, args.no_iters_burn_in)

if __name__ == '__main__':
    main(sys.argv[0:])
