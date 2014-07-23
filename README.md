# generateClonalityPipeline
A script that will create an LSF based analysis pipeline

### Make project directory
mkdir /lustre/scratch110/sanger/sd11/epitax

### Create project setup
python ~/repo/generateClonalityPipeline/setupProject.py -b /lustre/scratch110/sanger/sd11/epitax

### Create the required files
 * Create a file that lists all samplenames
 * Create a file that lists all tumournames
 * Create a file that lists which samples are "male" or "female"
 * Put the BAM files or symlinks in the bam directory

### Generate allele frequency loci for mutation calls
python ~/repo/generateClonalityPipeline/generateAlleleFrequencyLoci.py -i /lustre/scratch110/sanger/sd11/epitax/samplelist.txt -v /lustre/scratch110/sanger/sd11/epitax/variants/filtered_vcf -o /lustre/scratch110/sanger/sd11/epitax/haplotype/mutation_loci

### Generate GetAlleleFrequency scripts for the mutation call loci
python ~/repo/generateClonalityPipeline/generateAlleleFrequencyScripts.py -i /lustre/scratch110/sanger/sd11/epitax/samplelist.txt -r /lustre/scratch110/sanger/sd11/epitax/haplotype/mutations -b /lustre/scratch110/sanger/sd11/epitax/bam -l /lustre/scratch110/sanger/sd11/epitax/haplotype/mutation_loci --log /lustre/scratch110/sanger/sd11/epitax/haplotype/mutations/logs

### 1000 genomes loci
python ~/repo/generateClonalityPipeline/generateAlleleFrequencyScripts.py -i /lustre/scratch110/sanger/sd11/epitax/samplelist.txt -r /lustre/scratch110/sanger/sd11/epitax/haplotype/G1000 -b /lustre/scratch110/sanger/sd11/epitax/bam --log /lustre/scratch110/sanger/sd11/epitax/haplotype/G1000/logs

### Battenberg
> Don't forget to copy the GetAlleleFrequencies output for both tumour and normal into the sample Battenberg directory

python ~/repo/generateClonalityPipeline/generateBattenbergPipeline.py -i /lustre/scratch110/sanger/sd11/epitax/samplelist.small.txt -t a -n b -r /lustre/scratch110/sanger/sd11/epitax/battenberg/ -p ~/repo/battenberg/ --ismale /lustre/scratch110/sanger/sd11/epitax/ismale.small.txt

### Create symlinks to AlleleFrequencies files
python ~/repo/generateClonalityPipeline/generateAlleleFrequencySymlinks.py -i /lustre/scratch110/sanger/sd11/epitax/samplelist.genome.2del.txt -s /lustre/scratch110/sanger/sd11/epitax/haplotype/ -o /lustre/scratch110/sanger/sd11/epitax/battenberg/ --out_sample_subdir -t a -n b

### Generate dp input
> Note: this requires a sample and ismale list of just the tumours!

 python ~/repo/generateClonalityPipeline/generateDPInput.py -s /lustre/scratch110/sanger/sd11/epitax/samplelist.tumours.txt --ismale /lustre/scratch110/sanger/sd11/epitax/ismale.tumours.txt -r /lustre/scratch110/sanger/sd11/epitax/ -o /lustre/scratch110/sanger/sd11/epitax/dirichlet_input/
 
### Generate 1D Dirichlet Process pipeline
python ~/repo/generateClonalityPipeline/generateDirichlet1D.py -i /lustre/scratch110/sanger/sd11/epitax/samplelist.tumours.txt -b /lustre/scratch110/sanger/sd11/epitax/battenberg/ -r /lustre/scratch110/sanger/sd11/epitax/dirichlet_1d/ -d /lustre/scratch110/sanger/sd11/epitax/dirichlet_input/ --no_iters 1300 --no_iters_burn_in 300