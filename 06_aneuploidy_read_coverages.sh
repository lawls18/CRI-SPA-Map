##############################################################################
################################# SETUP ######################################
##############################################################################

# user-input files/info
HOME= ## PATH TO PROJECT DIRECTORY
ALIGNDIR=alignments/ ## name of directory for BAMS
CNVS='B3S1\|079C\|F7S2\|D3S3\|A3S1\|A3S4\|A3R2' ## unique identifiers of control and aneuploidy samples
BAMSUFF='_L001_rmdup.bam' ## ending of all read file names

# navigate to project directory
cd ${HOME}

# create output directories
READCOVDIR=aneuploid_coverage/ ## name of directory for read coverage files
mkdir ${READCOVDIR}

# load module files
module load parallel/20210822
module load bedtools/2.29.2

###############################################################################
############################## READ COVERAGE ##################################
###############################################################################

# extract read depth at all genome positions in a control and aneuploidy isolates
SAMPLES=$(ls ${ALIGNDIR} | grep ${CNVS} | grep "${BAMSUFF}$" | sed "s/${BAMSUFF}.*//") ## extract sample names from read files
parallel "bedtools genomecov -ibam ${ALIGNDIR}{}${BAMSUFF} -d > ${READCOVDIR}{}_covAll.txt" ::: ${SAMPLES}