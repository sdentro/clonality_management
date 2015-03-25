#!/usr/bin/env python

import sys, argparse
from path import path
from clonalityPipelineConfig import DIRICHLET_SUBDIR_NAME


def setupProject(basedir):
    print("Making output directories")
    path.joinpath(basedir,'samplesheet','input').makedirs()
    path.joinpath(basedir,'samplesheet','output').makedirs()
    path.joinpath(basedir,'bam').makedirs()
    path.joinpath(basedir,'battenberg').makedirs()
    path.joinpath(basedir,'variants').makedirs()
    path.joinpath(basedir,'dirichlet_preprocessing', 'input').makedirs()
    path.joinpath(basedir,'dirichlet_preprocessing', 'output').makedirs()
    path.joinpath(basedir,'dirichlet_input','logs').makedirs()
    path.joinpath(basedir,DIRICHLET_SUBDIR_NAME,'logs').makedirs()
    path.joinpath(basedir,'qc').makedirs()

    print("")
    print("Now the following must be done:")
    print("- Create a file that lists all samplenames")
    print("- Create a file that lists which samples are \"male\" or \"female\"")
    print("- Create a file that lists where the Battenberg input can be found for these samples")
    print("- Create a file that lists where the VCF files can be found")
    print("- Create a file that lists tumour bam files")
    print("- Create a file that lists normal bam files")
    print("- Put the BAM files or symlinks in the bam directory")
    print("- Put VCF files or symlinks in the variants directory")
    print("")

def main(argv):
    parser = argparse.ArgumentParser(prog='generateDPInput',
                             formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-b", required=True, type=str, help="Full path to basedir")
    args = parser.parse_args()
    setupProject(args.b)
    
if __name__ == '__main__':
    main(sys.argv[0:])
