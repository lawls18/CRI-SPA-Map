##############################################################################
################################# SETUP ######################################
##############################################################################

# user-input files/info
HOME=/projects/standard/albertf/lawle167/005_CSM_250828/ ## project directory
REF=/projects/standard/albertf/shared/genomes/sacCer_R64-5-1_240529/S288C_SGD_ref.fsa ## reference sequence
READDIR=/projects/standard/albertf/shared/MSI_data/2023-q2/230518_A00223_1060_AH7F7VDRX3/Albert_Project_027/ ## directory containing Illumina reads
READS='BYXW\|YFA\|YNL' ## unique identifier of samples
PAIR1='_R1' ## unique identifier of forward read file
PAIR2='_R2' ## unique identifier of reverse read file
READSUFF='_001.fastq.gz' ## ending of all read file names
CHROMS=(sacCer_chrI sacCer_chrII sacCer_chrIII sacCer_chrIV sacCer_chrV sacCer_chrVI sacCer_chrVII sacCer_chrVIII sacCer_chrIX sacCer_chrX sacCer_chrXI sacCer_chrXII sacCer_chrXIII sacCer_chrXIV sacCer_chrXV sacCer_chrXVI) ## names of nuclear chromosomes in reference fasta

# navigate to project directory
cd ${HOME}

# create output directories
ALIGNDIR=alignments/ ## name of directory for BAMS
COVDIR=coverage/ ## name of directory for read coverage files
VARDIR=variantCalls/ ## name of directory for VCFs
CNVDIR=CNVnator/ ## name of directory for CNVnator files
mkdir ${ALIGNDIR} ${COVDIR} ${VARDIR} ${CNVDIR}

# load module files
module load bwa/0.7.17
module load samtools/1.16.1
module load parallel/20210822
module load bcftools/1.16-gcc-8.2.0-5d4xg4y
module load cnvnator/0.4.1 

##############################################################################
############################# READ ALIGNMENT #################################
##############################################################################

# index reference sequence
bwa index -p ref ${REF}
samtools faidx ${REF}

# align reads to reference sequence, filter, collect coverage statistics, and index
SAMPLES=$(ls ${READDIR} | grep ${READS} | grep ${PAIR1} | sed "s/${PAIR1}.*//") ## extract sample names from read files
parallel "bwa mem -M -t 24 ref ${READDIR}{}${PAIR1}${READSUFF} ${READDIR}{}${PAIR2}${READSUFF} | 
samtools view -F 1284 -q 30 -d XS:0 -h - | 
samtools sort -n - | 
samtools fixmate -m - ${ALIGNDIR}{}_fixmate.bam 
samtools sort -@ 24 ${ALIGNDIR}{}_fixmate.bam | 
samtools markdup -r - ${ALIGNDIR}{}_rmdup.bam 
samtools index ${ALIGNDIR}{}_rmdup.bam 
samtools coverage -o ${COVDIR}{}_coverage.txt ${ALIGNDIR}{}_rmdup.bam
rm ${ALIGNDIR}{}_fixmate.bam" ::: ${SAMPLES}  ## remove intermediate file
rm ref* ## remove index files

###############################################################################
############################## VARIANT CALLING ################################
###############################################################################

# create file of file names for each VCF
find $PWD/${ALIGNDIR} | grep ".bam$" > ${VARDIR}FOFN.txt

# call biallelic variants, filter, normalize, and index
bcftools mpileup -a AD,DP -f ${REF} -Ou -b ${VARDIR}FOFN.txt | 
bcftools call -f GQ -mv -Ou --ploidy 1 | 
bcftools view -e 'QUAL<20' | 
bcftools norm -f ${REF} -m -both -Ob -o ${VARDIR}varCalls.bcf.gz 
bcftools index ${VARDIR}varCalls.bcf.gz 
bcftools view ${VARDIR}varCalls.bcf.gz > ${VARDIR}varCalls.txt ## convert VCF to .txt
rm ${VARDIR}FOFN.txt ## remove file of file names

###############################################################################
############################## READ COVERAGE ##################################
###############################################################################

# create individual chromosome fastas
mkdir ${CNVDIR}chr_fastas/
for chrom in ${CHROMS[@]}; do 
    samtools faidx -o ${CNVDIR}chr_fastas/${chrom}.fa ${REF} ${chrom}
done

# determine read coverage of 5000 bp regions using BAMs to identify any potential CNVs
parallel "cnvnator -root ${COVDIR}{}.root -chrom ${CHROMS[@]} -tree ${ALIGNDIR}{}_rmdup.bam
cnvnator -root ${CNVDIR}{}.root -his 5000 -d ${CNVDIR}chr_fastas
cnvnator -root ${CNVDIR}{}.root -stat 5000
cnvnator -root ${CNVDIR}{}.root -partition 5000
cnvnator -root ${CNVDIR}{}.root -call 5000 > ${CNVDIR}{}_cnv.txt" ::: ${SAMPLES} 
rm ${CNVDIR}*.root ## remove intermediate file
