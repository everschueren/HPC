#! /usr/bin/Rscript --vanilla --default-packages=utils

suppressMessages(library(sqldf))
suppressMessages(library(optparse))
suppressMessages(library(compiler))

## bio-conductor
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(GO.db))
suppressMessages(library(KEGG.db))
suppressMessages(library(PFAM.db))
suppressMessages(library(GSEABase))

## for hyperG tests
suppressMessages(library(annotate))
suppressMessages(library(genefilter))
suppressMessages(library(GOstats))

Enrichment.CUTOFF = 0.05
Enrichment.GO.EVIDENCE_TYPES = c('EXP','IDA','IPI','IMP','IGI','IEP','ISS','ISO','ISA','ISM','IGC','IBA','IBD','IKR','IRD','RCA','TAS','NAS','IC','ND','IEA')

##################
## annotations

# Experimental Evidence Codes
  # EXP: Inferred from Experiment
  # IDA: Inferred from Direct Assay
  # IPI: Inferred from Physical Interaction
  # IMP: Inferred from Mutant Phenotype
  # IGI: Inferred from Genetic Interaction
  # IEP: Inferred from Expression Pattern
# Computational Analysis Evidence Codes
  # ISS: Inferred from Sequence or Structural Similarity
  # ISO: Inferred from Sequence Orthology
  # ISA: Inferred from Sequence Alignment
  # ISM: Inferred from Sequence Model
  # IGC: Inferred from Genomic Context
  # IBA: Inferred from Biological aspect of Ancestor
  # IBD: Inferred from Biological aspect of Descendant
  # IKR: Inferred from Key Residues
  # IRD: Inferred from Rapid Divergence
  # RCA: inferred from Reviewed Computational Analysis
# Author Statement Evidence Codes
  # TAS: Traceable Author Statement
  # NAS: Non-traceable Author Statement
# Curator Statement Evidence Codes
  # IC: Inferred by Curator
  # ND: No biological Data available
# Automatically-assigned Evidence Codes
  # IEA: Inferred from Electronic Annotation
# Obsolete Evidence Codes
  # NR: Not Recorded

## BEGIN DEPRECATED? ###############################################################################################

Enrichment.firstGeneName = function(gene_name){
  unlist(str_split(gene_name, pattern=" "))[1]
}

Enrichment.uniprotKeys = function(){
  ## load uniprot set
  x = org.Hs.egUNIPROT
  mapped_genes = mappedkeys(x)
  unique(as.data.frame(x[mapped_genes])[,2])
}

Enrichment.getGODescription = function(df, id_idx=4){
  df = cbind(df, apply(df, 1, function(x) Term(x[id_idx])))
  colnames(df)[ncol(df)]='TERM'
  df
}

Enrichment.getPFAMDescription = function(df, id_idx=4){
  PF = PFAMDE
  mapped_keys = mappedkeys(PF)
  PF = as.data.frame(PF[mapped_keys])
  colname = colnames(df)[id_idx]
  res = sqldf(sprintf("select D.*, PF.de as 'TERM' from df D join PF PF on PF.ac = D.%s",colname))
  res
}

Enrichment.getKEGGDescription = function(df, id_idx=4){
  KEGG = KEGGPATHID2NAME
  mapped_keys = mappedkeys(KEGG)
  KEGG = as.data.frame(KEGG[mapped_keys])
  colname = colnames(df)[id_idx]
  res = sqldf(sprintf("select D.*, KG.path_name as 'TERM' from df D join KEGG KG on KG.path_id = D.%s",colname))
  res
}

Enrichment.uniprotToEntrez = function(ids){
  entrez_uniprot_map = org.Hs.egUNIPROT
  entrez_uniprot_map
}

## END DEPRECATED? ###############################################################################################

Enrichment.annotate = function(ids, db=org.Hs.eg.db,type="UNIPROT", cols=c("UNIPROT","GO","PFAM","PROSITE","PATH")){
  keys_ann = select(db, keys=ids, cols=cols, keytype=type)
  keys_ann
}


# Enrichment.getEntrezUniverse = function(mapping=org.Hs.egUNIPROT){
#   ## get the uniprot universe
#   Lkeys(mapping)
# }

Enrichment.getEntrezUniverse = function(){
  ## get the uniprot universe
  entrez_uniprot_map = org.Hs.egUNIPROT
  entrez_uniprot_map_keys = mappedkeys(entrez_uniprot_map)
  entrez_uniprot_keys = as.data.frame(entrez_uniprot_map[entrez_uniprot_map_keys])
  entrez_uniprot_keys_unique = unique(entrez_uniprot_keys$gene_id)  
  entrez_uniprot_keys_unique
}

Aux.listToDataFrame = function(l){
  res = c()
  for(t in names(l)){
    key = t
    val = l[[t]]
    for(v in val){
      res = rbind(res, c(t,v))
    }
  }
  as.data.frame(res) 
}

Aux.listToDataFrame = cmpfun(Aux.listToDataFrame)
  
Enrichment.getGOEnrichment = function(ids, db="org.Hs.eg.db", ontology=c("CC","BP","MF"), p_cutoff=Enrichment.CUTOFF){
  params = new("GOHyperGParams", 
               geneIds=ids,  
               universeGeneIds=Enrichment.getEntrezUniverse(),
               annotation=db,
               ontology=ontology, 
               pvalueCutoff=p_cutoff,
               conditional=F,
               testDirection="over")
  over_represented = hyperGTest(params)
  tmp = geneIdsByCategory(over_represented)
  tmp = Aux.listToDataFrame(tmp)
  colnames(tmp) = c("GO","ENTREZ")
  summ = summary(over_represented)
  goid = colnames(summ)[1]
  res = sqldf(sprintf("select S.*, T.ENTREZ from summ S left join tmp T on S.%s = T.GO", goid))
  res
}

Enrichment.getKEGGEnrichment = function(ids, db="org.Hs.eg.db", p_cutoff=Enrichment.CUTOFF){
  params = new("KEGGHyperGParams", 
               geneIds=ids,  
               universeGeneIds=Enrichment.getEntrezUniverse(),
               annotation=db,
               pvalueCutoff=p_cutoff,
               testDirection="over")
  over_represented = hyperGTest(params)
  tmp = geneIdsByCategory(over_represented)
  tmp = Aux.listToDataFrame(tmp)
  colnames(tmp) = c("KEGG","ENTREZ")
  summ = summary(over_represented)
  res = sqldf("select S.*, T.ENTREZ from summ S left join tmp T on S.KEGGID = T.KEGG")
  res
}

Enrichment.getPFAMEnrichment = function(ids, db="org.Hs.eg.db", p_cutoff=Enrichment.CUTOFF){
  params = new("PFAMHyperGParams", 
               geneIds=ids,  
               universeGeneIds=Enrichment.getEntrezUniverse(),
               annotation=db,
               pvalueCutoff=p_cutoff,
               testDirection="over")
  over_represented = hyperGTest(params)
  tmp = geneIdsByCategory(over_represented)
  tmp = Aux.listToDataFrame(tmp)
  colnames(tmp) = c("PFAM","ENTREZ")
  summ = summary(over_represented)
  
  PF = PFAMDE
  mapped_keys = mappedkeys(PF)
  PF = as.data.frame(PF[mapped_keys])
  res = sqldf("select S.PFAMID, S.Pvalue, S.OddsRatio, S.ExpCount, S.Count, S.Size, PF.de as 'Term', T.ENTREZ from summ S left join tmp T on PFAMID = T.PFAM left join PF PF on PF.ac = T.PFAM")
  res
}

Enrichment.group = function(enrichment_out, term_id_idx=2){
  out = enrichment_out
  for(e in names(enrichment_out)){
    enrichment = enrichment_out[[e]]
    enrichment_id_name = colnames(enrichment)[term_id_idx]
    enrichment_grouped = sqldf(sprintf("select setid, %s, Pvalue, OddsRatio, ExpCount, Count, Size, Term, group_concat(entrez_id) as 'entrez_id', group_concat(uniprot_ac) as 'uniprot_ac', group_concat(symbol) as 'symbol' from enrichment E group by %s order by Pvalue asc", enrichment_id_name, enrichment_id_name))
    out[[e]] = enrichment_grouped
  }
  out
}

Enrichment.main = function(data_file, output_dir, set_idx=1, prey_idx=2, grouped=T, enrichment_p_cutoff=enrichment_p_cutoff){
  
  df = read.delim(data_file, stringsAsFactors=F)

  GO_CC_groups = c()
  GO_BP_groups = c()
  GO_MF_groups = c()
  PFAM_groups = c()
  KEGG_groups = c()
  
  sets = unique(df[,set_idx])
  preys = unique(df[,prey_idx])
  
  for(i in 1:length(sets)){
    set = sets[i]
    print(set)
    
    ## GET ALL THE STUFF
    d = df[df[,set_idx]==set,]
    da = Enrichment.annotate(ids=unique(d[,prey_idx]), cols=c("UNIPROT","SYMBOL","ENTREZID"))
    da_entrez_ids = unique(da$ENTREZID)
    
    ###################
    ## ENRICHMENT WITH GO
    
    da_GO_CC_enriched = Enrichment.getGOEnrichment(da_entrez_ids, ontology="CC", p_cutoff=enrichment_p_cutoff)
    da_GO_CC_groups = sqldf(sprintf("select '%s' as 'setid', DAE.*, UNIPROT as 'uniprot_ac', ENTREZ as 'entrez_id', SYMBOL as 'symbol' from da_GO_CC_enriched DAE join da DA on (DA.ENTREZID = DAE.ENTREZ) order by DAE.Pvalue asc",set))
    GO_CC_groups = rbind(GO_CC_groups, da_GO_CC_groups)
    
    da_GO_MF_enriched = Enrichment.getGOEnrichment(da_entrez_ids, ontology="MF", p_cutoff=enrichment_p_cutoff)
    da_GO_MF_groups = sqldf(sprintf("select '%s' as 'setid', DAE.*, UNIPROT as 'uniprot_ac', ENTREZ as 'entrez_id', SYMBOL as 'symbol' from da_GO_MF_enriched DAE join da DA on (DA.ENTREZID = DAE.ENTREZ) order by DAE.Pvalue asc",set))
    GO_MF_groups = rbind(GO_MF_groups, da_GO_MF_groups)
    
    da_GO_BP_enriched = Enrichment.getGOEnrichment(da_entrez_ids, ontology="BP", p_cutoff=enrichment_p_cutoff)
    da_GO_BP_groups = sqldf(sprintf("select '%s' as 'setid', DAE.*, UNIPROT as 'uniprot_ac', ENTREZ as 'entrez_id', SYMBOL as 'symbol' from da_GO_BP_enriched DAE join da DA on (DA.ENTREZID = DAE.ENTREZ) order by DAE.Pvalue asc",set))
    GO_BP_groups = rbind(GO_BP_groups, da_GO_BP_groups)
    
    #####################
    ## GET THE PFAM STUFF
    
    da_PFAM_enriched = Enrichment.getPFAMEnrichment(da_entrez_ids, p_cutoff=enrichment_p_cutoff)
    da_PFAM_groups = sqldf(sprintf("select '%s' as 'setid', DAE.*, UNIPROT as 'uniprot_ac', ENTREZ as 'entrez_id', SYMBOL as 'symbol' from da_PFAM_enriched DAE join da DA on (DA.ENTREZID = DAE.ENTREZ) order by DAE.Pvalue asc",set))
    PFAM_groups = rbind(PFAM_groups, da_PFAM_groups)
    
    #####################
    ## GET THE KEGG STUFF
    da_KEGG_enriched = Enrichment.getKEGGEnrichment(da_entrez_ids, p_cutoff=enrichment_p_cutoff)
    da_KEGG_groups = sqldf(sprintf("select '%s' as 'setid', DAE.*, UNIPROT as 'uniprot_ac', ENTREZ as 'entrez_id', SYMBOL as 'symbol' from da_KEGG_enriched DAE join da DA on (DA.ENTREZID = DAE.ENTREZ) order by DAE.Pvalue asc",set))
    KEGG_groups = rbind(KEGG_groups, da_KEGG_groups)
  }
  
  res = list(GO_CC_groups=GO_CC_groups, GO_BP_groups=GO_BP_groups, GO_MF_groups=GO_MF_groups, PFAM_groups=PFAM_groups, KEGG_groups=KEGG_groups)
  
  if(grouped){
    res = Enrichment.group(res)
    suffix = "_grouped"
  }else{
    suffix = ""
  }
  
  for(l in names(res)){
    e = res[[l]] 
    write.table(e, file=sprintf("%s/%s%s.txt",output_dir,l,suffix), eol="\n", sep="\t", quote=F, row.names=F, col.names=T)
  }
}

##############
## MAIN ACTION

option_list <- list(
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print extra output [default]"),
  make_option(c("-q", "--quietly"), action="store_false",
              dest="verbose", help="Print little output"),
  make_option(c("-g", "--grouped"), default=TRUE,
              help="Flag to format output grouped per annotation term"),
  make_option(c("-d", "--data_file"),
              help="data file containing values"),
  make_option(c("-o", "--output_dir"),
              help="output file for enrichment"),
  make_option(c("-s", "--set_idx"), default=1,
              help="column with set descriptions"),
  make_option(c("-p", "--prey_idx"), default=2,
              help="column with genes/proteins as uniprot ID"),
  make_option(c("-e", "--enrichment_p_cutoff"), default=Enrichment.CUTOFF,
              help="p-value cutoff for enrichment")
)

parsedArgs = parse_args(OptionParser(option_list = option_list), args = commandArgs(trailingOnly=T))
Enrichment.main(data_file=parsedArgs$data_file, output_dir=parsedArgs$output_dir, set_idx=parsedArgs$set_idx, prey_idx=parsedArgs$prey_idx, grouped=parsedArgs$grouped, enrichment_p_cutoff=parsedArgs$enrichment_p_cutoff)  
