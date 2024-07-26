# Modified Tanya's code for build38

library(vroom)
library(dplyr)
library(tidyr)

# gencode annotation downloaded - Fri Jul 12 09:55:28 CDT 2024
# https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz
data = vroom("data/build38_genes/gencode.v46.annotation.gtf", col_names = F, comment = '#', quote = '')

# Column names acquired from: https://www.gencodegenes.org/pages/data_format.html
names(data) = c("chromosome_name", "annotation_source", "feature_type", "genomic_start", "genomic_end", "score", "genomic_strand", "genomic_phase", "tosplit")

# Convert final column into multiple columns (per annotation)
data$index = 1:nrow(data)

# Split list of annotations - new annotation per row (into label & value as
# a single column)
gene_ids = data %>% separate_rows(tosplit, sep = "; ")

# Split the row-separated tosplit column and remove the ";"
gene_ids = gene_ids %>% separate(tosplit, sep = ' ', into = c('label', 'value'))
gene_ids$value = gsub(';', '', gene_ids$value)

# Filter for relevant information first, then pivot the table to wide format:
gene_ids = gene_ids %>% filter(label %in% c("gene_id", "transcript_id", "protein_id", "gene_name", "gene_type", "transcript_type"))
tmp = gene_ids %>% pivot_wider(id_cols = index, names_from = label, values_from = value)

gene_ids = gene_ids %>% select(chromosome_name:genomic_phase, index) %>% distinct %>% left_join(., tmp, by = 'index')

# Create unique gene list for locuszooms etc
# Keep the relevant columns
gene_ids = gene_ids %>% select("chromosome_name", "genomic_start", "genomic_end", "genomic_strand", "gene_id", "transcript_id", "protein_id", "gene_name", "gene_type", "transcript_type", "feature_type")

# Keep lines for genes, transcripts, and CDS only
gene_ids = gene_ids[gene_ids$feature_type %in% c("gene", "transcript", "CDS"), ]

# Remove quotes from some of the string columns:
gene_ids$gene_id = gsub('"', '', gene_ids$gene_id)
gene_ids$transcript_id = gsub('"', '', gene_ids$transcript_id)
gene_ids$protein_id = gsub('"', '', gene_ids$protein_id)
gene_ids$gene_name = gsub('"', '', gene_ids$gene_name)
gene_ids$gene_type = gsub('"', '', gene_ids$gene_type)
gene_ids$transcript_type = gsub('"', '', gene_ids$transcript_type)

# Relabel gene_type categories (definitions from: https://www.gencodegenes.org/pages/biotypes.html)
# NOTE: "TEC" category is lefted as is - will be an option in the main LZ function
gene_type_proteincoding = c("protein_coding", "IG_C_gene", "IG_J_gene", "IG_D_gene", "IG_V_gene", "TR_C_gene", "TR_V_gene", "TR_J_gene", "TR_D_gene")
gene_type_psuedogene = c("transcribed_unprocessed_pseudogene", "unprocessed_pseudogene", "processed_pseudogene", "pseudogene", "transcribed_processed_pseudogene", "transcribed_unitary_pseudogene", "unitary_pseudogene", "translated_processed_pseudogene", "IG_C_pseudogene", "IG_pseudogene", "TR_J_pseudogene", "TR_V_pseudogene", "IG_V_pseudogene", "IG_J_pseudogene", "polymorphic_pseudogene", "rRNA_pseudogene")
gene_type_longnoncodingRNA = c("lncRNA",  "lincRNA", "sense_overlapping", "sense_intronic", "antisense", "processed_transcript")
gene_type_noncodingRNA = c("snRNA", "miRNA", "misc_RNA", "rRNA","scRNA", "Mt_rRNA", "Mt_tRNA", "snoRNA", "vault_RNA", "ribozyme", "scaRNA", "sRNA")

gene_ids$gene_type[gene_ids$gene_type %in% gene_type_proteincoding] = "proteincoding"
gene_ids$gene_type[gene_ids$gene_type %in% gene_type_psuedogene] = "psuedogene"
gene_ids$gene_type[gene_ids$gene_type %in% gene_type_longnoncodingRNA] = "lncRNA"
gene_ids$gene_type[gene_ids$gene_type %in% gene_type_noncodingRNA] = "ncRNA"

# Change the gene type of all the genes that have "ENSG" in the name;
gene_ids$gene_type[which(grepl('^ENSG', gene_ids$gene_name))] = "TEC"

# Remove those with no transcript_id:
gene_ids = gene_ids %>% filter(!is.na(transcript_id))

# Given a (character) vector, collapse the unique values into a single string
converge_str = function(vec, sep = '; ') {
	res = vec[!is.na(vec)]
	res = paste(unique(res), collapse = sep)
	return(res)
}

# Make a unique list of genes by summarising some of the genetic information
# from different trascripts
gene_ids = gene_ids %>% group_by(chromosome_name, gene_name, gene_type) %>% mutate(genomic_start = min(genomic_start),
												   genomic_end = max(genomic_end),
												   cds_start = ifelse('CDS' %in% feature_type, genomic_start, NA),
												   cds_end = ifelse('CDS' %in% feature_type, genomic_end, NA),
												   genomic_strand = converge_str(genomic_strand),
												   gene_id = converge_str(gene_id),
												   transcript_id = converge_str(transcript_id),
												   protein_id = converge_str(protein_id),
												   gene_name = converge_str(gene_name),
												   genomic_length = genomic_end - genomic_start,
												   cds_length = cds_end - cds_start,
												   coding = gene_type
												   ) %>% ungroup

# Select relevant columns and rename it to match the other data sets:
gene_ids = gene_ids %>% select(chromosome_name:genomic_end, cds_start:cds_end, genomic_strand, gene_id:protein_id, gene_name, genomic_length:coding)
colnames(gene_ids) = c("Chrom", "Start", "End", "cdsStart", "cdsEnd", "Strand", "ensemblGeneID", "ensemblTranscriptID", "ensemblProteinID", "Gene", "GeneLength", "cdsLength", "Coding")

gene_ids$Chrom = gsub('chr', '', gene_ids$Chrom)

gene_ids = gene_ids %>% distinct

# Sort data based on chr/pos
gene_ids = gene_ids %>% arrange(Chrom, Start)

# Deal with duplicated gene names
# List up duplicated genes
dup_genes = gene_ids %>% group_by(Gene) %>% summarize(count = n()) %>% arrange(desc(count)) %>% filter(count > 1) %>% ungroup

# Remove those that have 3 or more duplicates:
rm_list = dup_genes %>% filter(count > 2) %>% pull(Gene)
gene_ids = gene_ids %>% filter(!(Gene %in% rm_list))

# For those that are "true" duplicates, leave them in, since the "Coding" will
# be different between the duplicates and the LZ function can keep/remove
# different coding types

# Save
write.table(gene_ids, "data/build38_genes/Gencode_GRCh38_Genes_UniqueList2024.txt", quote = FALSE, sep = "\t", na = "", row.names = FALSE)

