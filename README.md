<a href="https://www.cell.com/cell/fulltext/S0092-8674(21)00294-4">link</a>


# generateClonalityPipeline
A series of example scripts/commands to create the full project setup.

## Docker build
docker build -t cm:0.0.1 .

mkdir sandbox
cd sandbox
cp -R ../test/ .
docker run -it -v /home/centos/docker_containers/clonality_management/sandbox:/mnt/testing cm:0.0.1

## Dry run

cd /mnt/testing/
setupProject.py -b $PWD/
cp test/* samplesheet/input/

cd samplesheet/input
python ~/repo/generateClonalityPipeline/generateSamplesheet.py -s samplenames.lst  --bt tumourbam.lst --bn normalbam.lst --idt tumourid.lst --idn normalid.lst -v variants.lst -b bb_dirs.lst --vcf_indel indels.lst -x sex.lst -o ../output/test.txt

## Setup the directory
cd /lustre/scratch110/sanger/sd11
mkdir dp_test
cd dp_test

python ~/repo/generateClonalityPipeline/setupProject.py -b $PWD

## Setup the bam files
cp ../epitax/bam/*txt bam/
cp ../epitax/bam/*sh bam/
cd bam
./create_links.sh 
cd ../

## Setup the variant files
cp ../epitax/variants/filtered_vcf/*gz variants/

## Setup BB
rm battenberg
ln -s ../epitax/battenberg battenberg

## Create samplesheet - This step is different for each project
find $PWD/variants | grep -v tbi | grep -v sh | grep -v filtered | sort | grep gz > samplesheet/input/variants.txt
for item in `cat samplesheet/input/variants.txt`; do sample=`basename ${item} | sed 's/a.filt.vcf.gz//'g`; echo ${sample}; done > samplesheet/input/samplenames.txt

for item in `cat samplesheet/input/samplenames.txt`; do echo $PWD/bam/${item}"a.bam"; done > samplesheet/input/bam_tumour.txt
for item in `cat samplesheet/input/samplenames.txt`; do echo ${item}"a"; done > samplesheet/input/id_tumour.txt
for item in `cat samplesheet/input/samplenames.txt`; do echo $PWD/bam/${item}"b.bam"; done > samplesheet/input/bam_normal.txt
for item in `cat samplesheet/input/samplenames.txt`; do echo ${item}"b"; done > samplesheet/input/id_normal.txt

for item in `cat samplesheet/input/samplenames.txt`; do echo "female"; done > samplesheet/input/gender.txt

for item in `cat samplesheet/input/samplenames.txt`; do echo "${PWD}/battenberg/${item}a"; done > samplesheet/input/bb_dirs.txt

python ~/repo/generateClonalityPipeline/generateSamplesheet.py -s samplesheet/input/samplenames.txt --bt samplesheet/input/bam_tumour.txt --idt samplesheet/input/id_tumour.txt --bn samplesheet/input/bam_normal.txt --idn samplesheet/input/id_normal.txt -x samplesheet/input/gender.txt -v samplesheet/input/variants.txt -b samplesheet/input/bb_dirs.txt -o samplesheet/output/test_samplesheet.txt

## Generate the pre-processing pipeline
python ~/repo/dirichlet_preprocessing/generate_inputfile.py --ss samplesheet/output/test_samplesheet.txt -o dirichlet_preprocessing/input/test.txt
cd dirichlet_preprocessing
python ~/repo/generateClonalityPipeline/generatePreprocessingPipeline.py -s input/test.txt -r output -f ${LUSTRE}/Documents/GenomeFiles/refs_icgc_pancan/genome.fa.fai -i ${LUSTRE}/Documents/GenomeFiles/battenberg_ignore/ignore.txt --ip ${LUSTRE}/Documents/GenomeFiles/battenberg_ignore/ignore_mut_cn_phasing.txt

## Submit and check on the pre-processing pipeline
./RunCommands.sh
python ~/repo/dirichlet_preprocessing/check_logs.py -s input/test.txt -r $PWD/output

## After finish copy the completed output files
cd ../
find dirichlet_preprocessing/output/ | grep allDiri | xargs -i cp {} dirichlet_input

## Create the DP master file
python ~/repo/generateClonalityPipeline/generateDPDataFile.py -p test -i samplesheet/output/test_samplesheet.txt -d $PWD/dirichlet_input/ -b $PWD/battenberg/ -r $PWD/dirichlet/

## Create the QC run script
python ~/repo/generateClonalityPipeline/generateQCrunScript.py -i $PWD/dirichlet/test.txt -d $PWD/dirichlet_input/ -q $PWD/qc/

## Create the run scripts
python ~/repo/generateClonalityPipeline/generateDPrunScript.py -i $PWD/dirichlet/test.txt -d $PWD/dirichlet_input/ -r $PWD/dirichlet -p test
