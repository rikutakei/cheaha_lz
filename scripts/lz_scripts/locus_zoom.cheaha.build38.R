# Usage example:
#
# Rscript src/major_lz.R data_name EUR gene   ABCG2     [optional bp offset]
# Rscript src/major_lz.R data_name EUR region 4         89052323  89062323
# Rscript src/major_lz.R data_name EUR snp    rs2231142 [optional bp offset]

library(vroom)

source('/data/user/home/rtakei/handy_scripts/locuszooms/scripts/lz_scripts/locus_zoom.function.R')

args = commandArgs(trailingOnly = T)

data_name = args[1]
ancestry = toupper(args[2])

dat = vroom(data_name)

# Change the base pair column name for Major TAMA GWAS from POS to BP
if (ancestry == 'TAMA') {
	colnames(dat)[3] = 'BP'
	p_type = 'BF'
} else {
	p_type = 'P'
}

# Load mapping data:
if (ancestry == 'UKB') {
	map_data = vroom('/data/user/home/rtakei/handy_scripts/locuszooms/data/build_mapping/all_ukbb_variants.txt')
} else {
	map_data = vroom('/data/user/home/rtakei/handy_scripts/locuszooms/data/build_mapping/all_1kgp_variants.txt')
}

gene = read.table('/data/user/home/rtakei/handy_scripts/locuszooms/data/build38_genes/Gencode_GRCh38_Genes_UniqueList2024.txt', sep = '\t', header = T, stringsAsFactors = F)

type = args[3]

out_dir = "/scratch/USER/lz_outputs/lz_plots/"
out_dir = gsub(pattern = 'USER', replacement = Sys.getenv('USER'), out_dir)
out_prefix = paste(out_dir, ancestry, sep = '')

if (type == 'gene') {
	out_name = paste(c(out_prefix, args[4], gsub('/', '', format(Sys.time(), '%D')), 'lz.jpg'), collapse = '_')
	locus.zoom(dat,
			   gene = args[4],
			   offset_bp = ifelse(length(args) == 4, 100000, args[5]),
			   genes.data = gene,
			   plot.title = paste(args[4], gsub('POP', ancestry, '(POP LD)')),
			   population = toupper(ancestry),
			   psuedogenes = F,
			   RNAs = F,
			   TEC = F,
			   sig.type = p_type,
			   rsid.check = F,
			   build38 = T,
			   map.data = map_data,
			   file.name = out_name)
} else if (type == 'region') {
	chr = as.numeric(args[4])
	start = as.numeric(args[5])
	end = as.numeric(args[6])
	region = c(chr, start, end)
	region_str = paste(chr, ':', start, '-', end, sep = '')
	out_name = paste(c(out_prefix, region, gsub('/', '', format(Sys.time(), '%D')), 'lz.jpg'), collapse = '_')
	locus.zoom(dat,
			   region = region,
			   offset_bp = 0,
			   genes.data = gene,
			   population = toupper(ancestry),
			   plot.title = paste(region_str, 'region', gsub('POP', ancestry, '(POP LD)')),
			   sig.type = p_type,
			   psuedogenes = F,
			   RNAs = F,
			   TEC = F,
			   rsid.check = F,
			   build38 = T,
			   map.data = map_data,
			   file.name = out_name)
} else if (type == 'snp'){
	out_name = paste(c(out_prefix, args[4], gsub('/', '', format(Sys.time(), '%D')), 'lz.jpg'), collapse = '_')
	locus.zoom(dat,
			   snp = args[4],
			   offset_bp = ifelse(length(args) == 4, 500000, args[5]),
			   genes.data = gene,
			   population = toupper(ancestry),
			   plot.title = paste(args[4], gsub('POP', ancestry, '(POP LD)')),
			   ignore.lead = T,
			   psuedogenes = F,
			   RNAs = F,
			   TEC = F,
			   sig.type = p_type,
			   rsid.check = F,
			   build38 = T,
			   map.data = map_data,
			   file.name = out_name)
}
