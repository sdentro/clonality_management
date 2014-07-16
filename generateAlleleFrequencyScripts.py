import sys, argparse, stat, os
from path import path

from clonalityPipelineConfig import PERL, ALLELECOUNTER, G1000LOCI
from generateClonalityPipeline_util import read_item_list, match_sample_to_file

def generateAlleleFrequencyScripts(infile, rundir, bamdir, locidir, logdir, isArrayJob):
    '''
    Creates a run script for each file specified in infile and lists bsub commands for each of those scripts in a single submit script.
    '''
    list_of_samplenames = read_item_list(infile)
    # Create a mapping from sample to loci file, depending on whether custom loci were supplied or not
    if locidir is not None:
        list_of_loci_files = path(locidir).listdir()
        mapping,_ = match_sample_to_file(list_of_samplenames, list_of_loci_files)
    else:
        mapping = dict()
        for sample in list_of_samplenames:
            mapping[sample] = G1000LOCI

    # Generate an allele frequency script for every combination
    commands_file = open(path.joinpath(rundir, "GetAlleleFrequenciesRunCommands.txt"), 'w')
    for sample,loci in mapping.iteritems():
        # TODO: throw error when mutliple loci files found!
        cmd = generateAlleleFrequencyScript(sample, rundir, bamdir, loci[0], logdir, isArrayJob)
        commands_file.write(cmd)
    commands_file.close()

def generateAlleleFrequencyScript(samplename, rundir, bamdir, locifile, logdir, isArrayJob):
    '''
    Creates a shell script for the sample. Returns the bsub command
    '''
    scriptfile = path.joinpath(rundir, 'GetAlleleFrequenciesFromBAMByChromosome'+samplename+'.sh')
    samplecommands = open(scriptfile,'w')
    samplecommands.write('#$LSB_JOBINDEX\n')
    if isArrayJob:
        samplecommands.write(PERL+' '+ALLELECOUNTER+' -b '+path.joinpath(bamdir,samplename+'.bam')+ ' -o '+path.joinpath(rundir,samplename+'_alleleFrequencies_chr$LSB_JOBINDEX.txt')+ ' -l '+locifile +' -m 20\n')
    else:
        samplecommands.write(PERL+' '+ALLELECOUNTER+' -b '+path.joinpath(bamdir,samplename+'.bam')+ ' -o '+path.joinpath(rundir,samplename+'_alleleFrequencies.txt')+ ' -l '+locifile +' -m 20\n')
    samplecommands.write('exit $?\n')
    samplecommands.close()
    st = os.stat(scriptfile)
    os.chmod(scriptfile, st.st_mode | stat.S_IEXEC)
    if isArrayJob:
        return("bsub -J\"GetAlleleFrequencies"+samplename+"[1-23]\" -o GetAlleleFrequencies_"+samplename+".%J.%I.out "+ path.joinpath(rundir,"GetAlleleFrequenciesFromBAMByChromosome"+samplename+".sh")+"\n")
    else:
        return("bsub -J\"GetAlleleFrequencies"+samplename+"\" -o " + path.joinpath(logdir, "GetAlleleFrequencies_"+samplename+".%J.out") +" "+ path.joinpath(rundir,"GetAlleleFrequenciesFromBAMByChromosome"+samplename+".sh")+"\n")


def main(argv):
    parser = argparse.ArgumentParser(prog='GenerateAlleleFrequencyScripts',
                             formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", required=True, type=str, help="Full path to input file containing a list of samplenames")
    parser.add_argument("-r", required=True, type=str, help="Full path to directory where the pipelines will be run")
    parser.add_argument("-b", required=True, type=str, help="Full path to directory where BAM files are residing")
    parser.add_argument("-l", type=str, help="Full path to directory where loci files are residing")
    parser.add_argument("--log", type=str, help="Full path to directory where log files should go")
#     parser.add_argument("--isarray", dest='isarray', action='store_true', help='Supply this option when submitting an array job')
    args = parser.parse_args()

    infile = args.i
    rundir = args.r
    bamdir = args.b
    locidir = args.l if not args.l is None else None
    isarray = False if not args.l is None else True # This is not working yet. When supplying manual loci, it is not an array job.
    logdir = args.log if not args.log is None else ''
    print(locidir)
    print(isarray)

    generateAlleleFrequencyScripts(infile, rundir, bamdir, locidir, logdir, isarray)

if __name__ == '__main__':
    main(sys.argv[0:])
