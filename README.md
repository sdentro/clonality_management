# generateClonalityPipeline
A script that will create an LSF based analysis pipeline

### Make project directory
mkdir /lustre/scratch110/sanger/sd11/epitax

### Create project setup
python ~/repo/generateClonalityPipeline/setupProject.py -b /lustre/scratch110/sanger/sd11/dp_test

### Create the required files
 * Create a file that lists per line (tab separated):
 	* individual id
 	* "male" or "female"
 	* comma separated list of normals
 	* comma separated list of tumours
 * Put the BAM files or symlinks in the bam directory

### Generate allele frequency loci for mutation calls
python ~/repo/generateClonalityPipeline/generateAlleleFrequencyLoci.py -i /lustre/scratch110/sanger/sd11/dp_test/samplesheet.txt -v /lustre/scratch110/sanger/sd11/dp_test/variants -o /lustre/scratch110/sanger/sd11/dp_test/haplotype/mutation_loci

### Generate GetAlleleFrequency scripts for the mutation call loci
python ~/repo/generateClonalityPipeline/generateAlleleFrequencyScripts.py -i /lustre/scratch110/sanger/sd11/dp_test/samplesheet.txt -r /lustre/scratch110/sanger/sd11/dp_test/haplotype/mutation_loci -b /lustre/scratch110/sanger/sd11/dp_test/bam -l /lustre/scratch110/sanger/sd11/dp_test/haplotype/mutation_loci --log /lustre/scratch110/sanger/sd11/dp_test/haplotype/mutations/logs --nonormals

### 1000 genomes loci
python ~/repo/generateClonalityPipeline/generateAlleleFrequencyScripts.py -i /lustre/scratch110/sanger/sd11/dp_test/samplesheet.txt -r /lustre/scratch110/sanger/sd11/dp_test/haplotype/G1000 -b /lustre/scratch110/sanger/sd11/dp_test/bam --log /lustre/scratch110/sanger/sd11/dp_test/haplotype/G1000/logs

### Battenberg
> Don't forget to copy the GetAlleleFrequencies output for both tumour and normal into the sample Battenberg directory

python ~/repo/generateClonalityPipeline/generateBattenbergPipeline.py -i /lustre/scratch110/sanger/sd11/dp_test/samplesheet.txt -r /lustre/scratch110/sanger/sd11/dp_test/battenberg/ -p ~/repo/battenberg/

### Battenberg - update params

python ~/repo/generateClonalityPipeline/generateBattenbergPipeline.py -i /lustre/scratch110/sanger/sd11/dp_test/samplesheet.txt -r /lustre/scratch110/sanger/sd11/dp_test/battenberg/ -p ~/repo/battenberg/ --rewrite_params

### Create symlinks to AlleleFrequencies files
> This interface has changed!

python ~/repo/generateClonalityPipeline/generateAlleleFrequencySymlinks.py -i /lustre/scratch110/sanger/sd11/dp_test/samplesheet.txt -s /lustre/scratch110/sanger/sd11/dp_test/haplotype/G1000/ -o /lustre/scratch110/sanger/sd11/dp_test/battenberg/ --out_sample_subdir

==== Tested to here ====

### Generate dp input

python ~/repo/generateClonalityPipeline/generateDPInput.py -i /lustre/scratch110/sanger/sd11/dp_test/samplesheet.txt -r /lustre/scratch110/sanger/sd11/dp_test/ -o /lustre/scratch110/sanger/sd11/dp_test/dirichlet_input/
 
### Generate 1D Dirichlet Process pipeline
> Note: deprecated

python ~/repo/generateClonalityPipeline/generateDirichlet1D.py -i /lustre/scratch110/sanger/sd11/epitax/samplelist.tumours.txt -b /lustre/scratch110/sanger/sd11/epitax/battenberg/ -r /lustre/scratch110/sanger/sd11/epitax/dirichlet_1d/ -d /lustre/scratch110/sanger/sd11/epitax/dirichlet_input/ --no_iters 1300 --no_iters_burn_in 300


### TODO: script that synchronises all the dp_input files per sample



### Generate dp data file that lists all the samples

python ~/repo/generateClonalityPipeline/generateDPDataFile.py -p dp_test -i /lustre/scratch110/sanger/sd11/dp_test/samplesheet.txt -b /lustre/scratch110/sanger/sd11/dp_test/battenberg/ -d /lustre/scratch110/sanger/sd11/dp_test/dirichlet_input/ -r /lustre/scratch110/sanger/sd11/dp_test/dirichlet/

## ICGC
### Battenberg input for ICGC
/lustre/scratch112/sanger/cancer_external/DBGap/TCGA_phs000178.v8.p7/analysis/Battenberg_sd11/scripts/createBattenbergInputFile.sh

### Battenberg ismale file for ICGC
/lustre/scratch112/sanger/cancer_external/DBGap/TCGA_phs000178.v8.p7/analysis/Battenberg_sd11/scripts/determineMaleFemale.sh > /lustre/scratch112/sanger/cancer_external/DBGap/TCGA_phs000178.v8.p7/analysis/Battenberg_sd11/ismale.train1.txt

### Battenberg pipelines
python ~/repo/generateClonalityPipeline/generateBattenbergPipeline.py -i /lustre/scratch112/sanger/cancer_external/DBGap/TCGA_phs000178.v8.p7/analysis/Battenberg_sd11/battenberg/battenberg_input_train1.txt -r $PWD -p ~/repo/battenberg/ --ismale /lustre/scratch112/sanger/cancer_external/DBGap/TCGA_phs000178.v8.p7/analysis/Battenberg_sd11/ismale.train1.txt

### Symlink allelefrequency files
python ~/repo/generateClonalityPipeline/generateAlleleFrequencySymlinks.py -i /lustre/scratch112/sanger/cancer_external/DBGap/TCGA_phs000178.v8.p7/analysis/Battenberg_sd11/battenberg/battenberg_input_train1.txt -s /lustre/scratch112/sanger/cancer_external/DBGap/TCGA_phs000178.v8.p7/analysis/Battenberg_sd11/haplotype/G1000/ -o /lustre/scratch112/sanger/cancer_external/DBGap/TCGA_phs000178.v8.p7/analysis/Battenberg_sd11/battenberg/ --out_sample_subdir


TODO: Write parser for ICGC file and a modifier that morphs the ICGC layout into something that fits with this layout