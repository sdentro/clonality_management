import argparse, sys, os, stat
from path import path

from util import merge_items

SCRIPT = "/software/R-3.1.0/bin/Rscript ~/repo/dirichlet/dp_combined/qc.R" 


def generateQcrunScript(dp_master_file, dp_in_dir, qc_dir):
    """
    Creates the qc runscript for easy qc generation.
    """
    scriptfile = path.joinpath(qc_dir, "runQc.sh")
    outf = open(scriptfile, "w")
    outf.write(merge_items([SCRIPT, dp_master_file, dp_in_dir, qc_dir])+"\n")
    outf.write("convert *alleleFreq*png alleleFrequency.pdf\n")
    outf.write("convert *copyNumberAdjustment*png copyNumberAdjustment.pdf\n")
    outf.write("convert *depth*png depth.pdf\n")
    outf.write("convert *kappa*png kappa.pdf\n")
    outf.write("convert *mutation.copy.num*png mutation.copy.number.pdf\n")
    outf.write("convert *totalCopy*png totalCopyNumber.pdf\n")
    outf.write("convert *_fractionOfCells*png fractionOfCells.pdf\n")
    outf.write("convert *subclonalFractionPerChromosome*png subclonalFractionPerChromosome.pdf\n")
    outf.write("convert *large.subclonal.fraction.by.chrom*png large.subclonal.fraction.by.chrom.pdf\n")
    outf.write("convert *depth.vs.frac.mutCount.png depth.vs.frac.mutCount.pdf\n")
    outf.write("convert *_cellularityCorrectedAF.png cellularityCorrectedAF.pdf\n")
    outf.close()
    # Make executable
    st = os.stat(scriptfile)
    os.chmod(scriptfile, st.st_mode | stat.S_IEXEC)

def main(argv):
    parser = argparse.ArgumentParser(prog='generateQCrunScript',
                             formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", required=True, type=str, help="Full path to dirichlet master file")
    parser.add_argument("-d", required=True, type=str, help="Full path to directory where dirichlet input files are stored")
    parser.add_argument("-q", required=True, type=str, help="Full path to directory where qc will be run")
    args = parser.parse_args()
    generateQcrunScript(args.i, args.d, args.q)
    
if __name__ == '__main__':
    main(sys.argv[0:])
