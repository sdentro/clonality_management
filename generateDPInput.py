import sys, argparse
from path import path
from util import merge_items

from clonalityPipelineConfig import GETDIRICHLETPROCESINFO_SCRIPT
from generateClonalityPipeline_util import read_sample_infile

#  R CMD BATCH '--no-restore-data --args PD7404a ../battenberg/PD7404a/PD7404a_rho_and_psi.txt ../haplotype/mutation_loci/PD7404a.loci ../battenberg/PD7404a/PD7404a_subclones.txt ../haplotype/mutations/PD7404a_alleleFrequencies.txt female' ~/repo/dirichlet/GetDirichletProcessInfo.R PD7404a.Rout

def generateDPInputs(infile, rundir, outdir):
#     samplenames = read_item_list(infile)
#     ismale = read_item_list(ismale)
    
    ss = read_sample_infile(infile)
    
    f = open(path.joinpath(outdir, "RunCommands.sh"),'w')
    for sample in ss.getSamplenames(): 

        normal = ss.getNormals(sample)[0]
        sex = ss.getSex(sample)
        
        for tumour in ss.getTumours(sample): # We've ran BB for each tumour specifically and we have variants for each tumour. Create DPInput for all
            cmd = generateDPInput(normal, tumour, rundir, sex, outdir)
            f.write(cmd+"\n")
    
    f.close()

def generateDPInput(normal, tumour, rundir, sex, outdir):
    cmd = ["bsub",
           "-J\"DPinput_"+tumour+"\"",
           "-o "+path.joinpath(outdir, 'logs',tumour+".%J.out"),
           "R CMD BATCH '--no-restore-data --args",
           tumour,
           
           this should be properly fixed
           
           #path.joinpath(rundir,'battenberg/'+normal+"_vs"+tumour,normal+"_vs"+tumour+'_rho_and_psi.txt'),
           path.joinpath(rundir,'battenberg/',tumour+'_rho_and_psi.txt'),
           path.joinpath(rundir,'haplotype/mutation_loci',tumour+'.loci'),
           #path.joinpath(rundir,'battenberg/'+normal+"_vs"+tumour,normal+"_vs"+tumour+'_subclones.txt'),
           path.joinpath(rundir,'battenberg/',tumour+'_subclones.txt'),
           path.joinpath(rundir,'haplotype/mutation_loci',tumour+'_alleleFrequencies.txt'),
           sex,
           outdir,
           "'",
           GETDIRICHLETPROCESINFO_SCRIPT,
           path.joinpath(outdir,tumour+".Rout")]
    
    return(merge_items(cmd))

def main(argv):
    parser = argparse.ArgumentParser(prog='generateDPInput',
                             formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", required=True, type=str, help="A sample sheet")
    parser.add_argument("-r", required=True, type=str, help="Base dir where pipelines were run")
    parser.add_argument("-o", required=True, type=str, help="Directory where to save the output")
    parser.add_argument("--bb_subdir", action='store_true', help="Give this flag when BB output is stored in a subdir.")
    parser.add_argument("--n_vs_t", )
    
    args = parser.parse_args()
    
    generateDPInputs(args.i, args.r, args.o)

if __name__ == '__main__':
    main(sys.argv[0:])