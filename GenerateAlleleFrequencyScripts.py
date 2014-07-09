import sys, argparse
from path import path

# Pointers to external files
PERL = 'perl-5.16.3 -I /software/CGP/pancan/lib/perl5'
ALLELECOUNTER = '/software/CGP/pancan/bin/alleleCounter.pl'

def _read_sample_list(infile):
    '''
    Reads in a list of samplenames specified one per line
    TODO: Add specification for what is tumour and normal to match these together here
    '''
    f = open(infile, 'r')

    samplenames = list()

    for line in f:
        l = line.strip()
        samplenames.append(l)

    return(samplenames)

def generateAlleleFrequencyScripts(infile, rundir, bamdir):
        '''
        Creates a run script for each file specified in infile and lists bsub commands for each of those scripts in a single submit script.
        '''
	list_of_samplenames = _read_sample_list(infile)

        commands_file = open(path.joinpath(rundir, "GetAlleleFrequenciesRunCommands.txt"), 'w')
        for samplename in list_of_samplenames:
            print(samplename)
            samplecommands = open(path.joinpath(rundir, 'GetAlleleFrequenciesFromBAMByChromosome'+samplename+'.sh'),'w')
            samplecommands.write('#$LSB_JOBINDEX\n')
            samplecommands.write(PERL+' '+ALLELECOUNTER+' -b '+path.joinpath(bamdir,samplename+'.bam')+ ' -o '+path.joinpath(rundir,samplename+'_alleleFrequencies_quality10_chr$LSB_JOBINDEX.txt')+ ' -l /lustre/scratch110/sanger/dw9/haplotype_pipeline/1000genomesloci2012_chr$LSB_JOBINDEX.txt -m 20\n')
            samplecommands.write('exit $?\n')
            samplecommands.close()
            commands_file.write("bsub -J\"GetAlleleFrequencies"+samplename+"[1-23]\" -o GetAlleleFrequencies_"+samplename+".%J.%I.out "+ path.joinpath(rundir,"GetAlleleFrequenciesFromBAMByChromosome"+samplename+".sh")+"\n")

        commands_file.close()



def main(argv):
    parser = argparse.ArgumentParser(prog='GenerateAlleleFrequencyScripts',
                             formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", required=True, type=str, help="Full path to input file containing a list of samplenames")
    parser.add_argument("-r", required=True, type=str, help="Full path to directory where the pipelines will be run")
    parser.add_argument("-b", required=True, type=str, help="Full path to directory where BAM files are residing")
    args = parser.parse_args()

    infile = args.i
    rundir = args.r
    bamdir = args.b

    generateAlleleFrequencyScripts(infile, rundir, bamdir)

if __name__ == '__main__':
    main(sys.argv[0:])
