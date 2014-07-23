import sys, argparse
from path import path

def setupProject(basedir):
    print("Making output directories")
    path.joinpath(basedir,'bam').makedirs()
    path.joinpath(basedir,'battenberg').makedirs()
    path.joinpath(basedir,'variants').makedirs()
    path.joinpath(basedir,'haplotype','logs').makedirs()
    path.joinpath(basedir,'dirichlet_input','logs').makedirs()
    path.joinpath(basedir,'dirichlet_1d','logs').makedirs()
    path.joinpath(basedir,'dirichlet_nd','logs').makedirs()
    path.joinpath(basedir,'dirichlet_tree','logs').makedirs()

    print("")
    print("Now the following must be done:")
    print("- Create a file that lists all samplenames")
    print("- Create a file that lists all tumournames")
    print("- Create a file that lists which samples are \"male\" or \"female\"")
    print("- Put the BAM files or symlinks in the bam directory")
    print("")

def main(argv):
    parser = argparse.ArgumentParser(prog='generateDPInput',
                             formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-b", required=True, type=str, help="Full path to basedir")
    args = parser.parse_args()
    setupProject(args.b)
    
if __name__ == '__main__':
    main(sys.argv[0:])