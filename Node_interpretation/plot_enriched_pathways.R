GO_file = './Node_interpretation/pseudomonas_GO_terms.txt'
KEGG_file = './Node_interpretation/pseudomonas_KEGG_terms.txt'
GO = read.table(GO_file,sep='\t',header=F,row.names=1,stringsAsFactors =F)
KEGG = read.table(KEGG_file,sep='\t',header=F,row.names=1,stringsAsFactors =F)
net_size = 50
gene_num = 5549

enrich_table = matrix(,nrow=0,ncol=5)
plot_folder = './Node_interpretation/GO_KEGG_plots/'
dir.create(plot_folder)
for (i in 1:net_size){
  HW_file = paste('./Train_test_DAs/high_weight_gene/2std/Node',i,'.txt',sep='')
  HW_genes = read.table(HW_file, header=T,sep='\t',stringsAsFactors =F)
  #GO enrichement analysis
  or_list = c()
  for (term in 1:nrow(GO)){
    term_genes = unlist(strsplit(GO[term,2],';'))
    high_in = length(intersect(HW_genes$gene, term_genes))
    low_in = length(term_genes) - high_in
    high_out = length(HW_genes$gene) - high_in
    low_out = gene_num - high_in - high_out - low_in
    odds_ratio = round((1.0*high_in / low_in) / (1.0*high_out / low_out),3)
    or_list = c(or_list, odds_ratio)
  }
  or_list_name = setNames(or_list, rownames(GO))
  GO_or_rank = sort(or_list_name,decreasing=T)
  #Plot top 10 enriched GO terms for each node
  plot_file = paste(plot_folder,'GO_Node',i,'.pdf',sep='')
  pdf(plot_file, height = 5, width = 15)
  par(mar = c(4,40,2,2.1))
  max_xlim = ifelse(max(GO_or_rank[10:1])>100,500,100)
  bp = barplot(GO_or_rank[10:1],cex.names=0.8, las=2, horiz=T,main = paste('Node',i,sep=''),xlab='Odds Ratio',xlim=c(1,max_xlim),log = "x")
  text(0,bp,round(GO_or_rank[10:1],2),cex=1,pos=4,size=3)
  dev.off()
  
  #KEGG enrichement analysis
  or_list = c()
  for (term in 1:nrow(KEGG)){
    term_genes = unlist(strsplit(KEGG[term,2],';'))
    high_in = length(intersect(HW_genes$gene, term_genes))
    low_in = length(term_genes) - high_in
    high_out = length(HW_genes$gene) - high_in
    low_out = gene_num - high_in - high_out - low_in
    odds_ratio = round((1.0*high_in / low_in) / (1.0*high_out / low_out),3)
    or_list = c(or_list, odds_ratio)
  }
  or_list_name = setNames(or_list, rownames(KEGG))
  KEGG_or_rank = sort(or_list_name,decreasing=T)
  #Plot top 10 enriched KEGG terms for each node
  plot_file = paste(plot_folder,'KEGG_Node',i,'.pdf',sep='')
  pdf(plot_file, height = 5, width = 15)
  par(mar = c(4,40,2,2.1))
  max_xlim = ifelse(max(KEGG_or_rank[10:1])>100,500,100)
  bp = barplot(KEGG_or_rank[10:1],cex.names=0.8, las=2, horiz=T,main = paste('Node',i,sep=''),xlab='Odds Ratio',xlim=c(1,max_xlim),log = "x")
  text(0,bp,round(KEGG_or_rank[10:1],2),cex=1,pos=4,size=3)
  dev.off()
  
  #combine GO and KEGG results into one table
  combined_table=(cbind(rep(paste('Node',i),10),names(GO_or_rank)[1:10],GO_or_rank[1:10],names(KEGG_or_rank)[1:10],KEGG_or_rank[1:10],row.names=NULL))
  rownames(combined_table) = NULL
  enrich_table = rbind(enrich_table,combined_table)
}

colnames(enrich_table) = c('Node','GO terms','Odds Ratio','KEGG terms','Odds Ratio')
write.table(enrich_table,'./Node_interpretation/GO_KEGG_enrichment.txt',quote=F,sep='\t',col.names=T,row.names=F)
