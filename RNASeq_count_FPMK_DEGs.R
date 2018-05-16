setwd('/Users/gloria/Documents/Project/pan-cancer-2018/analysis/mrna/Diff_Expression_Genes')
Cancer = 'PRAD'
#Inputs
Counts = read.csv(paste("/Volumes/Data/TCGA/geneexpression/counts/",Cancer,"_counts_tp_nt_paird.csv.data.txt",sep=''),sep=',',header = TRUE)
Counts_data = data.frame(Counts[,2:dim(Counts)[2]])
rownames(Counts_data) = Counts[,1]

FPKM = read.csv(paste("/Volumes/Data/TCGA/geneexpression/FPKM/",Cancer,"_FPKM_tp_nt_paird.csv.data.txt",sep=''),sep=',',header = TRUE)
FPKM_data = data.frame(FPKM[,2:dim(FPKM)[2]])
rownames(FPKM_data) = FPKM[,1]

#Sample Classification
samples = read.table(paste("/Volumes/Data/TCGA/geneexpression/FPKM/",Cancer,"_FPKM_tp_nt_paird.csv.sample.txt",sep=''),sep = '\t')
rownames(samples) = samples[,1]
tumor_sample = as.vector(samples[which(samples[,2]=='Tumor'),1])
nontumor_sample = as.vector(samples[which(samples[,2]=='Nontumor'),1])
tumor_group = gsub('-','.',tumor_sample)
nontumor_group = gsub('-','.',nontumor_sample)

#IDconversion
#source("https://bioconductor.org/biocLite.R")
#biocLite("biomaRt")
#library(biomaRt)
#ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
ids= substr(rownames(FPKM_data), 1, 15)
fids <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol','gene_biotype'), filters = 'ensembl_gene_id', values = ids, mart = ensembl)
protein_coding_genes = fids[which(fids[,'gene_biotype'] == 'protein_coding'),'ensembl_gene_id']
FPKM_protein_coding = FPKM_data[which(substr(rownames(FPKM_data),1,15) %in% protein_coding_genes),]
symbols = c(unlist(fids[,'hgnc_symbol']))
names(symbols) = c(unlist(fids[,'ensembl_gene_id']))

#Filtering
matrix_FPKM = FPKM_protein_coding
  H=0
  L=0
  sample_num = dim(matrix_FPKM)[2]
  detect=c()
  trans_expr=c()
  count=0
  data= matrix_FPKM
  namelist=rownames(matrix_FPKM)
  for(i in seq(1,dim(matrix_FPKM)[1])){
    count=0
    for (j in seq(1,dim(matrix_FPKM)[2]))
    {
      if (data[i,j] > 0){count = count + 1}
    }
    if (count >= 0.25* sample_num & count >= 6)
    {	
    H=H+1
    trans_expr[H]=namelist[i]
    }
    if (count > 0) 
    {
      L = L + 1
      detect[L]=namelist[i]
    }
  }
  length(detect)
  length(trans_expr)
  FPKM_expr_order= matrix_FPKM[trans_expr,]
  Count_expr_order = Counts_data[trans_expr,]


#DEGs
#library(edgeR)
  x= Count_expr_order
  rownames(x)= rownames(Count_expr_order)
  group <-factor(samples[,2])
  y <- DGEList(counts=x, group=group)
  design <- model.matrix(~0+group, data=y$samples)
  my.contrasts <- makeContrasts(
    TvsN = groupTumor - groupNontumor,
    levels=design)
  
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design)
  fit <- glmFit(y, design)
  lrt_TvsN <- glmLRT(fit, contrast=my.contrasts[,"TvsN"])
  result=(topTags(lrt_TvsN,p.value=0.01,n=10000))
  dim(result)
  t=unlist((result[,'logFC']))
  up_TvsN=result[which(as.numeric(t[1:(length(t)-3)]) > 1),]
  down_TvsN=result[which(as.numeric(t[1:(length(t)-3)]) < -1),]
  gene_names_up = as.vector(symbols[substr(rownames(up_TvsN),1,15)])
  gene_names_down = as.vector(symbols[substr(rownames(down_TvsN),1,15)])
  gene_names_deg = c(gene_names_up,gene_names_down)
  write.csv(gene_names_up,file=paste(Cancer,"up_genes_deg.csv",sep = ''))
  write.csv(gene_names_down,file = paste(Cancer,"down_genes_deg.csv",sep=''))
  write.csv(gene_names_deg,file = paste(Cancer,"deg_genes_deg.csv",sep=''))
  write.csv(up_TvsN,file=paste(Cancer,"RNASeq_TvsN_up.csv",sep=''))
  write.csv(down_TvsN,file=paste(Cancer,"RNASeq_TvsN_down.csv",sep=''))

