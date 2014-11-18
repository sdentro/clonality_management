# generateClonalityPipeline
A script that will create an LSF based analysis pipeline

### Make project directory
mkdir /lustre/scratch110/sanger/sd11/dp_test

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


## Preprocessing
mkdir /lustre/scratch110/sanger/sd11/dp_test
python ~/repo/generateClonalityPipeline/setupProject.py -b /lustre/scratch110/sanger/sd11/dp_test

cd /lustre/scratch110/sanger/sd11/dp_test/bam
ln -s /nfs/cancer_trk0019/00000071/187863.bam PD7404a.bam
ln -s /nfs/cancer_trk0016/00000057/187864.bam.bai PD7404a.bam.bai
ln -s /nfs/cancer_trk0021/00000076/187866.bam PD7404b.bam
ln -s /nfs/cancer_trk0016/00000057/187867.bam.bai PD7404b.bam.bai
find $PWD/ | grep -v bai | grep "a\.bam" > ../samplesheet/input/bam_tumour.txt
find $PWD/ | grep -v bai | grep "b\.bam" > ../samplesheet/input/bam_tumour.txt
ls $PWD/ | grep -v bai | grep "a\.bam" | sed 's/.bam//' > ../samplesheet/input/id_tumour.txt
ls $PWD/ | grep -v bai | grep "b\.bam" | sed 's/.bam//' > ../samplesheet/input/id_normal.txt

cd ../battenberg
ln -s /lustre/scratch110/sanger/sd11/epitax/battenberg/PD7404a PD7404a
echo $PWD/PD7404a > ../samplesheet/input/bb_dirs.txt

cd ../variants
ln -s /lustre/scratch110/sanger/sd11/epitax/variants/filtered_vcf/PD7404a.filt.vcf.gz PD7404a.filt.vcf.gz

find $PWD | grep gz > ../samplesheet/input/variants.txt
echo "female" > ../samplesheet/input/gender.txt
echo "PD7404a" > ../samplesheet/input/samplelist.txt

## Create samplesheet
cd ../samplesheet
python ~/repo/generateClonalityPipeline/generateSamplesheet.py -s input/samplelist.txt --bt input/bam_tumour.txt --idt input/id_tumour.txt --bn input/bam_normal.txt --idn input/id_normal.txt -x input/gender.txt -v input/variants.txt -b input/bb_dirs.txt -o output/epitax_samplesheet.txt

## Setup BB pipeline (optional)
python ~/repo/generateClonalityPipeline/generateBattenbergPipeline.py --ss ../samplesheet/output/epitax_samplesheet.txt -r $PWD -t 5

TODO: work with output of this BB (perhaps generate a list of bb_dirs per sample, or second samplesheet as can already be read in by preprocessing pipeline creation script)


## Create preprocessing pipeline
cd ../
### With samplesheet - it's possible to overwrite the bb dirs column by supplying -b with a list of bb dirs!
python ~/repo/dirichlet_preprocessing/generate_inputfile.py --ss samplesheet/output/epitax_samplesheet.txt -o dirichlet_preprocessing/input/epitax.txt
### Without samplesheet
python ~/repo/dirichlet_preprocessing/generate_inputfile.py -s samplesheet/input/id_tumour.txt -b samplesheet/input/bb_dirs.txt -v samplesheet/input/variants.txt -x samplesheet/input/gender.txt --bam samplesheet/input/bam_tumour.txt -o dirichlet_preprocessing/input/epitax.txt

cd dirichlet_preprocessing
python ~/repo/dirichlet_preprocessing/generatePreprocessingPipeline.py -s input/epitax.txt -r output
./RunCommands.sh
python ~/repo/dirichlet_preprocessing/check_logs.py -s ../samplesheet/input/id_tumour.txt -r $PWD/output

cp dirichlet_preprocessing/output/*/*allDirichlet* dirichlet_input/

python ~/repo/generateClonalityPipeline/generateDPDataFile.py -p dp_test -i /lustre/scratch110/sanger/sd11/dp_test/samplesheet.txt -b /lustre/scratch110/sanger/sd11/dp_test/battenberg/ -d /lustre/scratch110/sanger/sd11/dp_test/dirichlet_input/ -r /lustre/scratch110/sanger/sd11/dp_test/dirichlet/

### Create the QC input script
python ~/repo/generateClonalityPipeline/generateQCrunScript.py -i /nfs/users/nfs_c/cgppipe/pancancer/workspace/sd11/variant_calling_pilot_50_broad/dirichlet/icgc_train1_broad.txt -d /nfs/users/nfs_c/cgppipe/pancancer/workspace/sd11/variant_calling_pilot_50_broad/dirichlet_input/ -q /nfs/users/nfs_c/cgppipe/pancancer/workspace/sd11/variant_calling_pilot_50_broad/qc/
### Create the wrappers
python ~/repo/generateClonalityPipeline/generateDPrunScript.py -i /lustre/scratch112/sanger/cgppipe/PanCancerDownloads/workspace/sd11/variant_calling_pilot_50_broad/dirichlet/icgc_train1_broad.txt -d /lustre/scratch112/sanger/cgppipe/PanCancerDownloads/workspace/sd11/variant_calling_pilot_50_broad/dirichlet_input/ -r /lustre/scratch112/sanger/cgppipe/PanCancerDownloads/workspace/sd11/variant_calling_pilot_50_broad/dirichlet -p icgc

Remove:
alleleFrequency scripts
generateDPInput

Add:
cgpBB
checks for required programs 