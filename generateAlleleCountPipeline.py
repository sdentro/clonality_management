import sys, argparse, os, stat
from util import merge_items
from path import path
from generateClonalityPipeline_util import read_sample_infile, writeSimpleShellScript

ALLELECOUNTER = "alleleCounter"

def generateAlleleCountCommand(bam, loci, outfile):
    # Generate: alleleCounter -b bam/PD7422a.bam -l /lustre/scratch110/sanger/sd11/Documents/GenomeFiles/battenberg_1000genomesloci2012_v3/1000genomesloci2012_chr22.txt -o battenberg/PD7422a/^C7422a_alleleFrequencies_chr22.txt
    return(merge_items([ALLELECOUNTER,
                        "-b", bam,
                        "-l", loci,
                        "-o", outfile]))

def generateAlleleCountPipeline(run_dir, samplename, bam, loci, outfile):
    cmds = []
    cmds.append(generateAlleleCountCommand(bam, loci, outfile))
    return writeSimpleShellScript(run_dir, "alleleCount_"+samplename+".sh", cmds)

def generateAlleleCountPipelines(run_dir, samplesheet, loci):
    runscripts = []
    
    for sample in samplesheet.getSamplenames():
        print(sample)
        tn_pairs = samplesheet.getTumour2NormalPairingBam(sample)
        for tb,nb in tn_pairs:
            run_dir_sample = path.joinpath(run_dir, sample)
            
            runscript = generateAlleleCountPipeline(run_dir_sample, 
                                                    samplesheet.getIdByTumourBam(tb), 
                                                    tb, 
                                                    loci,
                                                    path.joinpath(run_dir_sample, samplesheet.getIdByTumourBam(tb)+"_alleleFrequencies_chr${LSB_JOBINDEX}.txt"))
            runscripts.append("bsub -q long -J alleleCount_"+samplesheet.getIdByTumourBam(tb)+"[1-23] -o "+path.joinpath(run_dir_sample, "logs", "alleleCount_"+samplesheet.getIdByTumourBam(tb)+".%J.out")+" "+runscript)
            
            runscript = generateAlleleCountPipeline(run_dir_sample, 
                                                    samplesheet.getIdByNormalBam(nb), 
                                                    nb, 
                                                    loci,
                                                    path.joinpath(run_dir_sample, samplesheet.getIdByNormalBam(nb)+"_alleleFrequencies_chr${LSB_JOBINDEX}.txt"))
            runscripts.append("bsub -q long -J alleleCount_"+samplesheet.getIdByNormalBam(nb)+"[1-23] -o "+path.joinpath(run_dir_sample, "logs", "alleleCount_"+samplesheet.getIdByNormalBam(nb)+".%J.out")+" "+runscript)
    
    
    # Create a master script
    scriptname = path.joinpath(run_dir, "RunAlleleCountCommands.sh")
    runscript = open(scriptname, 'w')
    for item in runscripts:
        runscript.write(item+"\n")
    runscript.close()
    
    # Make executable
    st = os.stat(scriptname)
    os.chmod(scriptname, st.st_mode | stat.S_IEXEC)
    

def main(argv):
    parser = argparse.ArgumentParser(prog='GenerateAlleleCountPipeline',
                             formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--ss", type=str, required=True, help="Full path to a samplesheet")
    parser.add_argument("-r", type=str, required=True, help="Full path to a directory where the pipelines will be ran")
    parser.add_argument("-l", type=str, help="Loci prefix")
    parser.add_argument("--l_postfix", type=str, help="Loci file postfix")
    parser.set_defaults(l_postfix=".txt", l="/lustre/scratch110/sanger/sd11/Documents/GenomeFiles/battenberg_1000genomesloci2012_v3/1000genomesloci2012_chr")
    
    args = parser.parse_args()

    ss = read_sample_infile(args.ss)
    generateAlleleCountPipelines(args.r, ss, args.l+"${LSB_JOBINDEX}"+args.l_postfix)

if __name__ == '__main__':
    main(sys.argv[0:])