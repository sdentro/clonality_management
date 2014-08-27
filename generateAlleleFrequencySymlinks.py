import sys, argparse
from path import path
from generateClonalityPipeline_util import read_sample_infile

def generateAlleleFrequencySymlinks(infile, in_dir, out_dir, out_sample_subdir, in_sample_subdir):
    """
    Creates symlinks of the alleleFrequency files in the out_dir pointing to the in_dir. The files
    can be stored in a sub directory. Use the sample_subdir flags to determine whether output and/or
    input files are in a subdir per sample
    
    """
    ss = read_sample_infile(infile)
    
    for samplename in ss.getSamplenames():
        for tumour in ss.getTumours(samplename):
            normal = ss.getNormals(samplename)[0]
        
            if in_sample_subdir:
                print("Unclear whether this is working properly")
                in_dir = path.joinpath(in_dir, tumour)
            if out_sample_subdir:
                outpath = path.joinpath(out_dir, normal+"_vs_"+tumour)    
    
            for sample in [tumour, normal]:
                alleleFrequencyFiles = path(in_dir).listdir(sample+'_allele*')
                if len(alleleFrequencyFiles) == 23:
                    for alleleFrequencyFile in alleleFrequencyFiles:
                        chrom = _get_chr_nr(alleleFrequencyFile)
                        alleleFrequencyFile.symlink(path.joinpath(outpath, sample+'_alleleFrequencies_chr'+chrom+'.txt'))
                else:
                    print("Found a different number of alleleFrequency files than expected for sample "+sample) 
                    print(alleleFrequencyFiles)

def _get_chr_nr(alleleFrequencyFile):
    """
    Assumes that the chromosome number is right before the extension of the file.
    """
    return alleleFrequencyFile.split("_")[-1].split(".")[0].strip('chr')

def main(argv):
    parser = argparse.ArgumentParser(prog='GenerateAlleleFrequencySymlinks',
                             formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", required=True, type=str, help="Full path to a sample sheet")
    parser.add_argument("-s", required=True, type=str, help="Full path to where AlleleFrequencies are stored")
    parser.add_argument("-o", required=True, type=str, help="Full path to where the symlinks should go")
    parser.add_argument("--out_sample_subdir", action='store_true', help="Supply this if symlinks should be created in a directory with the same name as the sample within the output directory")
    parser.add_argument("--in_sample_subdir", action='store_true', help="Supply this if symlinks should point to a directory with the same name as the sample within the input directory")
    parser.set_defaults(sample_subdir=False, sample_subdir_input=False)
    
    args = parser.parse_args()
    generateAlleleFrequencySymlinks(args.i, args.s, args.o, args.out_sample_subdir, args.in_sample_subdir)

if __name__ == '__main__':
    main(sys.argv[0:])