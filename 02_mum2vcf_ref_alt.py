# user-input file names and paths
refSeq_inFile = ## PATH TO S288C CHR XIV REFERENCE SEQUENCE
altSeq_inFile = ## PATH TO QUERY (BY4742 OR W303) CHR XIV SEQUENCE
var_inFile = ## PATH TO .SNPS FILE
vcf_outFile = ## PATH TO OUTPUT FILE

# initialize reference lists
refNames = []
refSeqs = []
cur_name = ''
cur_seq = ''

################## SEQUENCE PROCESSING ######################

# load in ref sequence
refSeq = open(refSeq_inFile, 'rt')

# separate ref chromosome names and sequences into parallel lists
for line in refSeq:
    line = line.strip()
    if line[0] == '>':
        if cur_name:
            refNames.append(cur_name)
            refSeqs.append(cur_seq)
        cur_name = line.split(' ')[0][1:]
        cur_seq = ''
    else:
        cur_seq += line

# append final ref chromosome name and sequence        
refNames.append(cur_name)
refSeqs.append(cur_seq)

# close ref sequence file
refSeq.close()

# initialize alt sequence lists
altNames = []
altSeqs = []
cur_name = ''
cur_seq = ''

# load in alt sequence
altSeq = open(altSeq_inFile, 'rt')

# separate alt chromosome names and sequences into parallel lists
for line in altSeq:
    line = line.strip()
    if line[0] == '>':
        if cur_name:
            altNames.append(cur_name)
            altSeqs.append(cur_seq)
        cur_name = line.split(' ')[0][1:]
        cur_seq = ''
    else:
        cur_seq += line

# append final alt chromosome name and sequence        
altNames.append(cur_name)
altSeqs.append(cur_seq)

# close alt sequence file
altSeq.close()
    
################## MUMMER TO VCF CONVERSION ###################

# load in mummer variant file
variants = open(var_inFile, 'rt')

# read in all variants and close file
var_lines = variants.readlines()
variants.close()
var_lines = var_lines[4:] # skips header lines

# initialize variables
vcf_vars = []
previous_cols = []

# for each variant
for line in var_lines:
   
    # separate information by columns
    columns = line.strip().split()
    
    # extract information desired for "VCF" file
    ref_pos = int(columns[0])
    ref_allele = columns[1]
    alt_allele = columns[2]
    alt_pos = int(columns[3])
    ref_chr = columns[10]
    alt_chr = columns[11]
    
    # if this is the first variant, store info
    if not previous_cols:
        previous_cols = [ref_pos, ref_allele, alt_allele, alt_pos, ref_chr, alt_chr]
    
    # if this is a ref deletion, concatenate alt alleles and track alt position
    elif ref_pos == previous_cols[0] and ref_allele == '.' and ref_chr == previous_cols[4]:
        previous_cols[2] += alt_allele
        previous_cols[3] = alt_pos
    
    # if this is a ref insertion, concatenate ref alleles and track ref position
    elif ref_pos == previous_cols[0] + 1 and alt_allele == '.' and ref_chr == previous_cols[4]:
        previous_cols[1] += ref_allele
        previous_cols[0] = ref_pos
    
    else:
        
        # if previous variant was a ref insertion
        if previous_cols[2] == '.':
            ref_chr_ind = refNames.index(previous_cols[4]) # determine which ref chromosome variant is on
            previous_cols[0] = previous_cols[0] - len(previous_cols[1]) # determine the ref position of the first inserted base 
            previous_cols[1] = refSeqs[ref_chr_ind][previous_cols[0] - 1] + previous_cols[1] # extract the base right before the insertion from the ref sequence and add it to the beginning of the ref allele
            alt_chr_ind = altNames.index(previous_cols[5]) # determine which alt chromosome variant is on
            previous_cols[2] = altSeqs[alt_chr_ind][previous_cols[3] - 1] # extract base before insertion from alt reference and list as alt allele
       
        # if previous variant was a ref deletion
        if previous_cols[1] == '.':
            alt_chr_ind = altNames.index(previous_cols[5]) # determine which alt chromosome variant is on
            previous_cols[3] = previous_cols[3] - len(previous_cols[2]) # determine the alt position of the first inserted base 
            previous_cols[2] = altSeqs[alt_chr_ind][previous_cols[3] - 1] + previous_cols[2] # extract the base right before the insertion from the alt reference and add it to the beginning of the alt allele
            ref_chr_ind = refNames.index(previous_cols[4]) # determine which ref chromosome variant is on
            previous_cols[1] = refSeqs[ref_chr_ind][previous_cols[0] - 1] # extract base before insertion from ref sequence and list as ref allele
       
        previous_cols = map(str, previous_cols) # convert all elements to strings for writing
        
        vcf_vars.append(previous_cols) # add previous variant to compiled "VCF" variant list
        
        previous_cols = [ref_pos, ref_allele, alt_allele, alt_pos, ref_chr, alt_chr] # reset tracker with current variant information

# add final variant information to compiled "VCF" variant list
previous_cols = map(str, previous_cols)
vcf_vars.append(previous_cols)

# open output file for writing
vcf = open(vcf_outFile, 'wt')

# write a header line for column names
vcf.write('ref_pos\tref_allele\talt_allele\talt_pos\tref_chr\talt_chr\n')

# write "VCF" variants to file
for variant in vcf_vars:
    vcf.writelines('\t'.join(variant) + '\n')

# close output file
vcf.close()