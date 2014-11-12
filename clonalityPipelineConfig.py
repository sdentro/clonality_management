
# generateBattenbergPipeline
IMPUTEFILESDIR='/lustre/scratch110/sanger/sd11/Documents/GenomeFiles/battenberg_impute'
G1000LOCIDIR='/lustre/scratch110/sanger/sd11/Documents/GenomeFiles/battenberg_1000genomesloci2012'

# generateAlleleFrequencyScripts
PERL = 'perl-5.16.3 -I /software/CGP/pancan/lib/perl5'
ALLELECOUNTER = '/software/CGP/pancan/bin/alleleCounter.pl'
ALLELECOUNTER_BAQ = 1
ALLELECOUNTER_MAQ = 1
G1000LOCI = '/lustre/scratch110/sanger/sd11/Documents/GenomeFiles/battenberg_1000genomesloci2012/1000genomesloci2012_chr$LSB_JOBINDEX.txt'

# generateDPInput
GETDIRICHLETPROCESINFO_SCRIPT = "~/repo/dirichlet/GetDirichletProcessInfo.R"

# generateDirichlet1D
RHO_AND_PSI_REGEX = '*rho_and_psi*'
DIRICHLET_1D_RUNSCRIPT = "~/repo/dirichlet/dp_1d/RunDirichlet_1D_pipeline.R"

# generateDPDataFile
DIRICHLET_SUBDIR_NAME = "dirichlet"


