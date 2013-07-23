#! /usr/bin/Rscript --vanilla --default-packages=utils

INSTALL_DIR = Sys.getenv("MS_PIPELINE_PATH") 
source(paste(INSTALL_DIR,"src/Conversion/AnnotateWithUniprot_lib.R",sep=""))

suppressMessages(library(optparse))
suppressMessages(library(compiler))

option_list <- list(
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print extra output [default]"),
  make_option(c("-q", "--quietly"), action="store_false",
              dest="verbose", help="Print little output"),
  make_option(c("-d", "--data_file"),
              help="Data file with limma-like output (Protein, Peptide, LFC, p-value)"),
  make_option(c("-o", "--output_file"),
              help="Output file"),
  make_option(c("-u", "--uniprot_dir"),
              help="Organisms for which descriptions are available in the uniprot_dir"),
  make_option(c("-s", "--species"),
              help="Species for which descriptions are available in the uniprot_dir"),
  make_option(c("-k", "--key"),
              help="column containing the key in the data file")
)

parsedArgs = parse_args(OptionParser(option_list = option_list), args = commandArgs(trailingOnly=T))
annotate_with_uniprot <- cmpfun(annotate_with_uniprot)
data = read.delim(parsedArgs$data_file, stringsAsFactors=F)
annotate_with_uniprot(data, species=parsedArgs$species, key=parsedArgs$key, output_file=parsedArgs$output_file, uniprot_dir=parsedArgs$uniprot_dir)  
