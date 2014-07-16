import sys, argparse
from path import path
from util import merge_items, run_command

from generateClonalityPipeline_util import read_item_list

#  R CMD BATCH '--no-restore-data --args PD7404a ../battenberg/PD7404a/PD7404a_rho_and_psi.txt ../haplotype/mutation_loci/PD7404a.loci ../battenberg/PD7404a/PD7404a_subclones.txt ../haplotype/mutations/PD7404a_alleleFrequencies.txt female' ~/repo/dirichlet/GetDirichletProcessInfo.R PD7404a.Rout
GETDIRICHLETPROCESINFO_SCRIPT = "~/repo/dirichlet/GetDirichletProcessInfo.R"

def generateDPInputs(infile, rundir, ismale, outdir):
    samplenames = read_item_list(infile)
    ismale = read_item_list(ismale)
    f = open(path.joinpath(outdir, "RunCommands.sh"),'w')
    for item in zip(samplenames, ismale):
        cmd = generateDPInput(item[0], rundir, item[1], outdir)
        f.write(cmd+"\n")
    f.close()

def generateDPInput(samplename, rundir, ismale, outdir):
    cmd = ["bsub",
           "-J\"DPinput_"+samplename+"\"",
           "-o "+path.joinpath(outdir, 'logs',samplename+".%J.out"),
           "R CMD BATCH '--no-restore-data --args",
           samplename,
           path.joinpath(rundir,'battenberg/'+samplename,samplename+'_rho_and_psi.txt'),
           path.joinpath(rundir,'haplotype/mutation_loci',samplename+'.loci'),
           path.joinpath(rundir,'battenberg/'+samplename,samplename+'_subclones.txt'),
           path.joinpath(rundir,'haplotype/mutations',samplename+'_alleleFrequencies.txt'),
           ismale,
           outdir,
           "'",
           GETDIRICHLETPROCESINFO_SCRIPT,
           path.joinpath(outdir,samplename+".Rout")]
    
    return(merge_items(cmd))

def main(argv):
    parser = argparse.ArgumentParser(prog='generateDPInput',
                             formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-s", required=True, type=str, help="File containing samplename")
    parser.add_argument("-r", required=True, type=str, help="Base dir where pipelines were run")
    parser.add_argument("-o", required=True, type=str, help="Directory where to save the output")
    parser.add_argument("--ismale", required=True, type=str, help='File containing whether each samplename is male or female.')
    
    args = parser.parse_args()
    
    generateDPInputs(args.s, args.r, args.ismale, args.o)

if __name__ == '__main__':
    main(sys.argv[0:])