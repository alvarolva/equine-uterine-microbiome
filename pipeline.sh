#!/bin/bash

#############################
###1.IMPORT data in qimme2###
#############################
#1.1 Demultiplex
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest.csv \
  --output-path paired-end-demux.qza \
  --input-format PairedEndFastqManifestPhred33 ;

#1.2 Remove primers and adapters before runing
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences paired-end-demux.qza \
  --p-cores 7 \
  --p-front-f CCTACGGGNBGCASCAG \
  --p-front-r NGACTACNVGGGTATCTAATCC \
  --p-match-adapter-wildcards \
  --p-match-read-wildcards \
  --p-discard-untrimmed \
  --o-trimmed-sequences trimmed-demux.qza ;

#1.3. Sequence quality verification and number of reads per samples
qiime demux summarize \
  --i-data trimmed-demux.qza \
  --p-n 100000 \
  --o-visualization trimmed-demux.qzv ;

############################
###2. DENOISED WITH DADA2###
###########################
mkdir dada2_files
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs trimmed-demux.qza \
  --p-trunc-len-f 276 \
  --p-trunc-len-r 214 \
  --o-table dada2_files/table.qza \
  --o-representative-sequences dada2_files/rep-seqs.qza \
  --o-denoising-stats dada2_files/denoising-stats.qza \
  --p-n-threads 7 ;

# Denoising stats visualization
qiime metadata tabulate \
  --m-input-file dada2_files/denoising-stats.qza \
  --o-visualization dada2_files/denoising-stats.qzv ;

###################################
####3. Taxonomic classification####
###################################
#3.1 Summary of merged data
qiime feature-table tabulate-seqs \
  --i-data dada2_files/rep-seqs.qza \
  --o-visualization dada2_files/rep-seqs.qzv ;

#3.2. Non prokariotic sequence removal
qiime quality-control exclude-seqs \
  --i-query-sequences dada2_files/rep-seqs.qza \
  --i-reference-sequences silva_97_sequences.qza \
  --p-method vsearch \
  --p-perc-identity 0.8 \
  --p-perc-query-aligned 0.5 \
  --p-threads 7 \
  --o-sequence-hits dada2_files/filtered-seqs.qza \
  --o-sequence-misses dada2_files/misses.qza ;

qiime feature-table filter-features \
  --i-table dada2_files/table.qza \
  --m-metadata-file dada2_files/filtered-seqs.qza \
  --o-filtered-table dada2_files/filtered-table.qza ;

#3.3 ASVs TAXONOMY. Using ML classifier (sklearn).
qiime feature-classifier classify-sklearn \
  --i-classifier 2022.10.backbone.v4.nb.qza \
  --i-reads dada2_files/filtered-seqs.qza \
  --o-classification taxonomy.qza \
  --p-n-jobs 0 ;


qiime tools export \
  --input-path taxonomy.qza \
  --output-path exported_taxonomy ;

# visualize the taxonomy assignments
qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy_summary.qzv ;

#3.3. Chloroplast and mitochondrial 16S removal
qiime taxa filter-table \
  --i-table dada2_files/filtered-table.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude Chloroplast,mitochondria \
  --o-filtered-table final-table.qza ;

qiime tools export \
  --input-path final-table.qza \
  --output-path exported_final-table ;

biom convert \
  -i exported_final-table/feature-table.biom \
  -o exported_final-table/initial.feature-table.tsv \
  --to-tsv


####################################
####4.Filtering neglectable ASVs####
####################################
#4.1: Delete ASVs with min frequency < 5 in at least 5 samples
qiime feature-table filter-features \
  --i-table final-table.qza \
  --p-min-frequency 5 \
  --p-min-samples 5 \
  --o-filtered-table dada2_files/final-table-5fmin-5smin.qza ;

#4.2: Get table with all deleted ASVs
qiime tools export \
  --input-path dada2_files/final-table-5fmin-5smin.qza \
  --output-path dada2_files/exported_final-table-5fmin-5smin
biom convert \
  -i dada2_files/exported_final-table-5fmin-5smin/feature-table.biom \
  -o dada2_files/exported_final-table-5fmin-5smin/final-table-5fmin-5smin.tsv \
  --to-tsv ;

#4.3: Generate file of sequences to exclude
cut -f 1 dada2_files/exported_final-table-5fmin-5smin/final-table-5fmin-5smin.tsv | tail -n +2 > dada2_files/filtered_sequences_to_exclude.txt
qiime feature-table filter-features \
  --i-table final-table.qza \
  --m-metadata-file dada2_files/filtered_sequences_to_exclude.txt \
  --p-exclude-ids \
  --o-filtered-table dada2_files/deleted-table-5fmin-5smin.qza ;

#4.4: Select those previously excluded ASVs with min frequency >= 20 in at least 1 sample
qiime feature-table filter-features \
  --i-table dada2_files/deleted-table-5fmin-5smin.qza \
  --p-min-samples 1 \
  --p-min-frequency 20 \
  --o-filtered-table dada2_files/filtered-table-20fmin-1smin-from5-5.qza ;

#4.5: Transform the data and merge the files
qiime feature-table merge \
  --i-tables dada2_files/final-table-5fmin-5smin.qza dada2_files/filtered-table-20fmin-1smin-from5-5.qza \
  --p-overlap-method sum \
  --o-merged-table dada2_files/final-filtered-table-5-5-20-1.qza ;

#4.6: After all the filtering, file with all the sequences found
qiime feature-table filter-seqs \
  --i-data dada2_files/filtered-seqs.qza \
  --i-table dada2_files/final-filtered-table-5-5-20-1.qza \
  --p-no-exclude-ids \
  --o-filtered-data dada2_files/final-rep-seqs.qza ;

qiime tools export \
  --input-path dada2_files/final-filtered-table-5-5-20-1.qza \
  --output-path exported_final-table ;

biom convert \
  -i exported_final-table/feature-table.biom \
  -o exported_final-table/final-table.tsv --to-tsv


#Table visualization
qiime feature-table summarize \
  --i-table dada2_files/final-filtered-table-5-5-20-1.qza\
  --o-visualization final-table-summary.qzv ;
qiime feature-table tabulate-seqs \
  --i-data dada2_files/final-rep-seqs.qza  \
  --o-visualization final-table-counts.qzv ;

########################
####6. Phylogenetic ####
########################
mkdir phylogenytree
#6.1. Multiple sequence alignment with MAFFT
qiime alignment mafft \
  --i-sequences dada2_files/final-rep-seqs.qza  \
  --o-alignment phylogenytree/final-aligned-rep-seqs.qza \
  --p-n-threads 7 ;
#6.2. Masking the hypevariable positions in the alignement
qiime alignment mask \
  --i-alignment phylogenytree/final-aligned-rep-seqs.qza \
  --o-masked-alignment phylogenytree/final-masked-aligned-rep-seqs.qza ;
#6.3. Build phylogenetic tree with FastTree
qiime phylogeny fasttree \
  --i-alignment phylogenytree/final-masked-aligned-rep-seqs.qza \
  --o-tree phylogenytree/unrooted-tree.qza \
  --p-n-threads 7 ;
#6.4. Tree rooting
qiime phylogeny midpoint-root \
  --i-tree phylogenytree/unrooted-tree.qza \
  --o-rooted-tree phylogenytree/rooted-tree.qza ;
#6.5Tree exportation
qiime tools export \
  --input-path phylogenytree/unrooted-tree.qza \
  --output-path phylogenytree/exported-unrooted-tree/
qiime tools export \
  --input-path phylogenytree/rooted-tree.qza \
  --output-path phylogenytree/exported-rooted-tree/

######################################
####5. Alpha and Beta biodiversity####
######################################
mkdir diversity
depth=1118227

#5.1 Alpha biodiversitymetrics and rarefaction plot
qiime diversity alpha-rarefaction \
  --i-table dada2_files/final-filtered-table-5-5-20-1.qza \
  --i-phylogeny phylogenytree/rooted-tree.qza \
  --p-max-depth $depth \
  --m-metadata-file metadata.tsv \
  --o-visualization alpha-rarefaction.qzv
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny phylogenytree/rooted-tree.qza  \
  --i-table dada2_files/final-filtered-table-5-5-20-1.qza \
  --p-sampling-depth $depth \
  --m-metadata-file metadata.tsv \
  --output-dir diversity/core-metrics \
  --p-n-jobs-or-threads 4 ;

#Shanon
qiime diversity alpha \
  --i-table dada2_files/final-filtered-table-5-5-20-1.qza \
  --p-metric shannon \
  --o-alpha-diversity diversity/core-metrics/shannon_vector.qza
qiime tools export \
  --input-path diversity/core-metrics/shannon_vector.qza \
  --output-path diversity/core-metrics/alpha.shannon_vector_exported
qiime diversity alpha-group-significance \
  --i-alpha-diversity diversity/core-metrics/shannon_vector.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization diversity/core-metrics/shannon-group-significance.qzv

#Simpson
qiime diversity alpha \
  --i-table dada2_files/final-filtered-table-5-5-20-1.qza \
  --p-metric simpson \
  --o-alpha-diversity diversity/core-metrics/simpson_vector.qza
qiime tools export \
  --input-path diversity/core-metrics/simpson_vector.qza \
  --output-path diversity/core-metrics/alpha.simpson_vector_exported
qiime diversity alpha-group-significance \
  --i-alpha-diversity diversity/core-metrics/simpson_vector.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization diversity/core-metrics/simpson-group-significance.qzv

#Chao1
qiime diversity alpha \
  --i-table dada2_files/final-filtered-table-5-5-20-1.qza \
  --p-metric chao1 \
  --o-alpha-diversity diversity/core-metrics/chao1_vector.qza
qiime tools export \
  --input-path diversity/core-metrics/chao1_vector.qza \
  --output-path diversity/core-metrics/alpha.schao1_vector_exported
qiime diversity alpha-group-significance \
 --i-alpha-diversity diversity/core-metrics/chao1_vector.qza \
 --m-metadata-file metadata.tsv \
 --o-visualization diversity/core-metrics/chao1-group-significance.qzv

#Summary all the alpha ALPHA_DIVERSITYSTATS
ALPHA=$(ls diversity/core-metrics/alpha.*/*.tsv)
for file in $ALPHA; do
    sort -o "$file" -k1,1 "$file"
done
paste -d '\t' $ALPHA | cut -f1,2,4,6- | tail -n +2 | (echo -e "sample-id\tchao1\tshannon_entropy\tsimpson"; cat) > merged_alpha-diversity.tsv
head -n 1 metadata.tsv > temp_header.tsv
tail -n +2 metadata.tsv | grep -v "^#" | sort -o temp_sorted_content.tsv -k1,1
cat temp_header.tsv temp_sorted_content.tsv > metadata.sorted.tsv
rm temp_header.tsv temp_sorted_content.tsv
join metadata.sorted.tsv merged_alpha-diversity.tsv > ALPHA_DIVERSITYSTATS.txt


#Beta diversity (JACCARD)
qiime diversity beta-group-significance \
  --i-distance-matrix diversity/core-metrics/jaccard_distance_matrix.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column dna \
  --o-visualization diversity/core-metrics/jaccard_significance.qzv \
  --p-pairwise ;
qiime diversity adonis \
  --i-distance-matrix diversity/core-metrics/jaccard_distance_matrix.qza \
  --m-metadata-file metadata.tsv \
  --p-formula dna \
  --o-visualization diversity/core-metrics/jaccard_ADONIS.qzv ;
qiime tools export \
  --input-path diversity/core-metrics/jaccard_pcoa_results.qza \
  --output-path diversity/core-metrics/jaccard_pcoa_exported_pcoa

#Relative frequency table generation
mkdir diversity/rel-freq ;
levels="2 3 4 5 6 7"
for i in $levels; do
 qiime taxa collapse \
  --i-table dada2_files/final-filtered-table-5-5-20-1.qza \
  --i-taxonomy taxonomy.qza \
  --p-level ${i} \
  --o-collapsed-table diversity/rel-freq/table-l${i}.qza
  qiime feature-table relative-frequency \
    --i-table diversity/rel-freq/table-l${i}.qza \
    --o-relative-frequency-table diversity/rel-freq/table-rel-freq-l${i}.qza
  qiime tools export \
    --input-path diversity/rel-freq/table-rel-freq-l${i}.qza \
    --output-path diversity/rel-freq/l${i}-rel-freq-table
  biom convert \
  -i diversity/rel-freq/l${i}-rel-freq-table/feature-table.biom \
  -o diversity/rel-freq/l${i}-rel-freq-table.tsv --to-tsv
done ;

#Absolute frequency table generation
mkdir diversity/abs-freq
for i in $levels; do
  qiime composition add-pseudocount \
   --i-table diversity/rel-freq/table-l${i}.qza \
   --o-composition-table diversity/abs-freq/composition-l${i}.qza
  qiime tools export \
    --input-path diversity/abs-freq/composition-l${i}.qza \
    --output-path diversity/abs-freq/composition-l${i}.csv
  biom convert \
    -i diversity/abs-freq/composition-l${i}.csv/feature-table.biom \
    -o diversity/abs-freq/l${i}-abs-freq-table.tsv --to-tsv
done ;
