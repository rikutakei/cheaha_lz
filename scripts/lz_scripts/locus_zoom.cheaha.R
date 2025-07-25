# Usage example:
#
# Rscript src/major_lz.R data_name EUR gene   ABCG2     [optional bp offset]
# Rscript src/major_lz.R data_name EUR region 4         89052323  89062323
# Rscript src/major_lz.R data_name EUR snp    rs2231142 [optional bp offset]

library(vroom)

whoami = Sys.getenv('USER')

# Load in LZ functios:
source(gsub('whoami', whoami, '/data/user/home/whoami/handy_scripts/locuszooms/scripts/lz_scripts/locus_zoom.function.R'))

# Start processing arguments:
args = commandArgs(trailingOnly = T)

# Load data;
data_name = args[1]
dat = vroom(data_name)

# Set ancestry:
ancestry = toupper(args[2])

# Change the base pair column name for Major TAMA GWAS from POS to BP
if (ancestry == 'TAMA') {
	colnames(dat)[3] = 'BP'
	p_type = 'BF'
} else {
	p_type = 'P'
}

# Set output directory:
out_dir = "/scratch/USER/lz_outputs/lz_plots/"
out_dir = gsub(pattern = 'USER', replacement = whoami, out_dir)
out_prefix = paste(out_dir, ancestry, sep = '')

# Deal with character chromosomes:
dat$CHR[dat$CHR == 'X'] = 23
dat$CHR[dat$CHR == 'Y'] = 24
dat$CHR[dat$CHR == 'MT'] = 25
dat$CHR = as.numeric(dat$CHR)

# Make option to load different build version:
gen_build = as.numeric(gsub('[a-zA-Z]', '', args[3]))
if (gen_build == 38) {
	gene = read.table(gsub('whoami', whoami, '/data/user/home/whoami/handy_scripts/locuszooms/data/build38_genes/Gencode_GRCh38_Genes_UniqueList2024.txt'), sep = '\t', header = T, stringsAsFactors = F)
} else {
	gene = read.table(gsub('whoami', whoami, '/data/user/home/whoami/handy_scripts/locuszooms/data/build37_genes/Gencode_GRCh37_Genes_UniqueList2021.txt'), sep = '\t', header = T, stringsAsFactors = F)
}

# Deal with SNP/gene/region argument better:
type = args[4]

input_snp = NA
input_gene = NA
input_region = NA
lead_ignore = F
offset = 0

if (type == 'gene') {
	input_gene = args[5]
	input_str = input_gene
	offset = as.numeric(ifelse(length(args) == 5, 100000, args[6]))
} else if (type == 'region') {
	chr = as.numeric(args[5])
	start = as.numeric(args[6])
	end = as.numeric(args[7])
	input_region = c(chr, start, end)
	input_str = paste(chr, ':', start, '-', end, sep = '')
} else if (type == 'snp'){
	input_snp = args[5]
	input_str = input_snp
	lead_ignore = T
	offset = as.numeric(ifelse(length(args) == 5, 500000, args[6]))
}

out_name = paste(c(out_prefix, input_str, gsub('/', '', format(Sys.time(), '%D')), 'lz.jpg'), collapse = '_')

locus.zoom(dat,
		   snp = input_snp,
		   region = input_region,
		   gene = input_gene,
		   offset_bp = offset,
		   genes.data = gene,
		   plot.title = paste(input_str, gsub('POP', ancestry, '(POP LD)')),
		   population = toupper(ancestry),
		   psuedogenes = F,
		   RNAs = F,
		   TEC = F,
		   sig.type = p_type,
		   rsid.check = F,
		   build = gen_build,
		   ignore.lead = lead_ignore,
		   file.name = out_name
		   )

