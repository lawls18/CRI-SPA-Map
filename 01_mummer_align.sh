# user-input files
REF= ## PATH TO S288C REFERENCE SEQUENCE
BY= ## PATH TO BY4742 CHR XIV SEQUENCE
W303= ## PATH TO W303 CHR XIV SEQUENCE

# load module files
module load samtools/1.16.1
module load mummer/3.23

# create output directory
OUTDIR= ## DESIRED NAME FOR OUTPUT DIRECTORY
mkdir -p ${OUTDIR}

# extract chrXIV sequence from S288C and BY4742 files
samtools faidx -o ${OUTDIR}S288C_chrXIV.fasta ${REF} sacCer_chrXIV
samtools faidx -o ${OUTDIR}BY4742_chrXIV.fasta ${BY} chr14

# create Mummer delta files comparing BY4742 and W303 to S288C
nucmer --maxmatch -l 15 -g 10000 -b 10000 -p ${OUTDIR}S288C_BY_mummer ${OUTDIR}S288C_chrXIV.fasta ${OUTDIR}BY4742_chrXIV.fasta
nucmer --maxmatch -l 15 -g 10000 -b 10000 -p ${OUTDIR}S288C_W303_mummer ${OUTDIR}S288C_chrXIV.fasta ${W303}

# filter for non-overlapping alignments and set length requirement to 400
delta-filter -1 ${OUTDIR}S288C_BY_mummer.delta > ${OUTDIR}S288C_BY_mummer_filtered.delta
delta-filter -1 ${OUTDIR}S288C_W303_mummer.delta > ${OUTDIR}S288C_W303_mummer_filtered.delta

# output SNPs to tab-delimited tables
show-snps -ClrT ${OUTDIR}S288C_BY_mummer_filtered.delta > ${OUTDIR}S288C_BY_mummer.snps
show-snps -ClrT ${OUTDIR}S288C_W303_mummer_filtered.delta > ${OUTDIR}S288C_W303_mummer.snps