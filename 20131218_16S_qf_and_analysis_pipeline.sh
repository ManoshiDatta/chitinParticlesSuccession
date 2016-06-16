#!/bin/bash

# 16S demultiplexing/quality filtering/analysis pipeline
# Manoshi Sen Datta
# 13 November 2014

##################################################################################################################################################
## PATHS FOR PROGRAMS/DATA FOR RUNNING SCRIPTS ##

FILE_PATH="/Users/manoshi/Dropbox (MIT)/Pooled_particles_paper/Data/16S"

## PROGRAMS REQUIRED FOR ANALYSIS
# Path for Python scripts
PYTHONPATH="${FILE_PATH}"/scripts

# Path for Drive5 Python scripts
DRIVE5PATH=/Users/manoshi/Documents/drive5_py/

# Path for USEARCH
UPATH=/Users/manoshi/Documents/usearch

# Path for MOTHUR
MOTHURPATH=/Users/manoshi/Documents/Mothur.source

# Path for QIIME
QIIMEPATH=/Users/manoshi/Documents/MacQIIME_1.8.0-20140103_OS10.6/macqiime/configs/bash_profile.txt

# Taxonomic classification data path
TAXPATH=/Users/manoshi/Documents/Tax.classification

# RDP classifier path
RDPPATH=/Users/manoshi/Documents/rdp_classifier_2.3

## FILES FROM SEQUENCING RUN

# Path for forward reads (fastq file)
ILLFILE_F="${FILE_PATH}"/131218Alm_D13-6845_1_sequence.fastq

# Path for reverse reads (fastq file)
ILLFILE_R="${FILE_PATH}"/131218Alm_D13-6845_2_sequence.fastq

## SEQUENCING METADATA FILES

# Mapping file path
MAPPATH="${FILE_PATH}"/map_file_all.txt

# Barcode file path
BARCODEPATH="${FILE_PATH}"/barcodes_all.txt

OLIGOS_F="${FILE_PATH}"/oligos_F.txt
OLIGOS_R="${FILE_PATH}"/oligos_R.txt

FWDBARCODE=TGGTC

## OUTPUT INFORMATION

# SampleID
SAMPLEID=bulkPartSea

# OTUID
OTUID=97

#OTUPERCENT
OTUPERCENT=0.97

# Output path
OUTPUTPATH="${FILE_PATH}"/output
<<COMMENT
###################################################################################################################################################
# Demultiplex by forward barcode (separate my samples from everyone else's)
echo "Demultiplexing samples by forward barcode"

python "${PYTHONPATH}"/demultiplex_by_fwd_barcode_2.py -f "${ILLFILE_F}" -r "${ILLFILE_R}" -fo "${OUTPUTPATH}"/${SAMPLEID}_1.fastq -ro "${OUTPUTPATH}"/${SAMPLEID}_2.fastq -b ${FWDBARCODE}
COMMENT
###################################################################################################################################################
# Merge paired-end reads
echo "Merging paired-end reads"

${UPATH} -fastq_mergepairs "${OUTPUTPATH}"/${SAMPLEID}_1.fastq -reverse "${OUTPUTPATH}"/${SAMPLEID}_2.fastq -fastq_minmergelen 295 -fastq_maxmergelen 298 -fastq_maxdiffs 0 -fastq_qmin 2 -fastq_qmax 41 -fastq_ascii 64 -fastqout "${OUTPUTPATH}"/${SAMPLEID}.merged.fastq
###################################################################################################################################################
# Look at fastq statistics
echo "Look at fastq statistics"

${UPATH} -fastq_stats "${OUTPUTPATH}"/${SAMPLEID}.merged.fastq -log "${OUTPUTPATH}"/fastq.stats.log -fastq_qmin 2 -fastq_qmax 41 -fastq_ascii 64
###################################################################################################################################################
# Quality filter reads, convert to FASTA
echo "Quality filtering reads"

${UPATH} -fastq_filter "${OUTPUTPATH}"/${SAMPLEID}.merged.fastq -fastq_maxee 0.1 -fastaout "${OUTPUTPATH}"/${SAMPLEID}.merged.filtered.fasta
##################################################################################################################################################
# Trim primer sequences from merged sequences
echo "Trimming primer sequences"

cd ${MOTHURPATH}
# Trim primer sequences
./mothur "#trim.seqs(fasta=${OUTPUTPATH}/${SAMPLEID}.merged.filtered.fasta, oligos=${OLIGOS_F})"
# Summary of trimmed sequences
./mothur "#summary.seqs(fasta=${OUTPUTPATH}/${SAMPLEID}.merged.filtered.trim.fasta)"
##################################################################################################################################################
# Convert sequence headers to format compatible with USEARCH
echo "Converting sequence headers to USEARCH format"

python "${PYTHONPATH}"/changeHeader.py -b "${BARCODEPATH}" -i "${OUTPUTPATH}"/${SAMPLEID}.merged.filtered.trim.fasta -o "${OUTPUTPATH}"/${SAMPLEID}.merged.filtered.trim.header.fasta
##################################################################################################################################################
# Dereplicate reads
echo "Dereplicating reads"

${UPATH} -derep_fulllength "${OUTPUTPATH}"/${SAMPLEID}.merged.filtered.trim.header.fasta -output "${OUTPUTPATH}"/${SAMPLEID}.merged.filtered.trim.header.derep.fasta -sizeout
##################################################################################################################################################
# Abundance sort and discard singletons

echo "Abundance sorting and discarding singletons"

${UPATH} -sortbysize "${OUTPUTPATH}"/${SAMPLEID}.merged.filtered.trim.header.derep.fasta -output "${OUTPUTPATH}"/${SAMPLEID}.merged.filtered.trim.header.derep.sorted.fasta -minsize 1
##################################################################################################################################################
# Cluster OTUs

echo "Clustering OTUs"

${UPATH} -cluster_otus "${OUTPUTPATH}"/${SAMPLEID}.merged.filtered.trim.header.derep.sorted.fasta -otus "${OUTPUTPATH}"/${SAMPLEID}.otus.${OTUID}.fasta -otuid ${OTUPERCENT}
##################################################################################################################################################
# Filter out chimeras with reference database
echo "Filtering out chimeras"

${UPATH} -uchime_ref "${OUTPUTPATH}"/${SAMPLEID}.otus.${OTUID}.fasta -db ${TAXPATH}/gold.fa -strand plus -nonchimeras "${OUTPUTPATH}"/${SAMPLEID}.otus.${OTUID}.nonchimeras.fasta
##################################################################################################################################################
# Label OTU sequences (OTU_1, OTU_2, etc.)

python ${DRIVE5PATH}/fasta_number.py "${OUTPUTPATH}"/${SAMPLEID}.otus.${OTUID}.nonchimeras.fasta OTU_ > "${OUTPUTPATH}"/${SAMPLEID}.otus.${OTUID}.nonchimeras.labeled.fasta
##################################################################################################################################################
# Find taxonomic classification for OTUs

echo "Finding taxonomic classification for OTUs" >> ${OUTPUTINFOPATH}

java -jar ${RDPPATH}/rdp_classifier-2.3.jar -q "${OUTPUTPATH}"/${SAMPLEID}.otus.${OTUID}.nonchimeras.labeled.fasta -o "${OUTPUTPATH}"/${SAMPLEID}.otus.${OTUID}.nonchimeras.labeled.rdp.fasta
##################################################################################################################################################
# Clean taxonomic classification

python "${PYTHONPATH}"/clean_rdp.py -i "${OUTPUTPATH}"/${SAMPLEID}.otus.${OTUID}.nonchimeras.labeled.rdp.fasta -t 0.8 > "${OUTPUTPATH}"/${SAMPLEID}.otus.${OTUID}.nonchimeras.labeled.rdp.cleaned.fasta
##################################################################################################################################################
# Map reads (including singletons) back to OTUs

${UPATH} -usearch_global "${OUTPUTPATH}"/${SAMPLEID}.merged.filtered.trim.header.fasta -db "${OUTPUTPATH}"/${SAMPLEID}.otus.${OTUID}.nonchimeras.labeled.fasta -strand plus -id ${OTUPERCENT} -uc "${OUTPUTPATH}"/${SAMPLEID}.otus.${OTUID}.map.uc
##################################################################################################################################################
# Create OTU table

python ${DRIVE5PATH}/uc2otutab.py "${OUTPUTPATH}"/${SAMPLEID}.otus.${OTUID}.map.uc > "${OUTPUTPATH}"/${SAMPLEID}.otus.${OTUID}.otu_table.txt
##################################################################################################################################################
# Change headers on OTU table to reflect samples

python "${PYTHONPATH}"/change_OTU_table_headers.py -b "${BARCODEPATH}" -i "${OUTPUTPATH}"/${SAMPLEID}.otus.${OTUID}.otu_table.txt -o "${OUTPUTPATH}"/${SAMPLEID}.otus.${OTUID}.otu_table.headers.txt
##################################################################################################################################################