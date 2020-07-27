library(GNET2)

library(clusterProfiler)
library(STRINGdb)
library(minet)
library(GENIE3)
library(PRROC)


# library(igraph)
# library(xgboost)
# library(matrixStats)
# library(reshape2)
# library(dplyr)

r_seed = 121

exp_data <- read.csv('data/GSE79676_rpkm_table_ercc.csv',as.is = T,row.names = 1)
tf_list <- read.delim('data/Transcription_Factor_Trisomy_1.txt',
                      as.is = T,header = T,sep = '\t')
tf_list1 <- tf_list$gene_id




exp_data1 <- log2(exp_data[apply(exp_data, 1, function(x)sum(x==0))<9,]+1)
tf_list2 <- intersect(rownames(exp_data1),tf_list1)



regulators_list <- tf_list2[1:10]
string_db <- STRINGdb$new( version="11", species=3702, score_threshold=400, input_directory="")

for (i in regulators_list) {
  example1_mapped <- string_db$map( data.frame('gene'=regulators_list),'gene',removeUnmappedRows = T)
  
}


xxx <- string_db$get_neighbors( unique(example1_mapped$STRING_id) )
all_genes_network_string <- unique(c(xxx,example1_mapped$STRING_id))
yyy <- string_db$get_interactions( all_genes_network_string )
yyy <- yyy[yyy$from %in% example1_mapped$STRING_id,]


el_true <- cbind.data.frame('from'=gsub('.+\\.(.+)\\..+','\\1',yyy$from),
                            'to'=gsub('.+\\.(.+)\\..+','\\1',yyy$to),'score'=yyy$combined_score)

el_true <- unique(el_true)
el_true <- el_true[el_true$from %in% regulators_list,]
el_true <- el_true[el_true$to %in% rownames(exp_data1),]



true_network_regulators <- intersect(tf_list2,el_true$from)
all_genes_network <- intersect(unique(c(el_true$from,el_true$to)),rownames(exp_data1))


true_network_matrix <- matrix(0,nrow = length(true_network_regulators),ncol = length(all_genes_network))
rownames(true_network_matrix) <- true_network_regulators
colnames(true_network_matrix) <- all_genes_network

for (i in 1:nrow(el_true)) {
  true_network_matrix[el_true$from[i],el_true$to[i]] <- 1
}
exp_data2 <- exp_data[all_genes_network,]
print(dim(true_network_matrix))
print(sum(true_network_matrix))


#Run benchmark


set.seed(r_seed)
kmeans_start <- Sys.time()
gnet_results_kmeans <- gnet(input = exp_data2,reg_names = true_network_regulators,init_method = 'kmeans',
                            init_group_num = 10)
kmeans_end <- Sys.time()
print(kmeans_end-kmeans_start)

gbdt_start <- Sys.time()
gnet_results_gbdt <- gnet(input = exp_data2,reg_names = true_network_regulators,init_method = 'boosting' ,
                          init_group_num = 10)
gbdt_end <- Sys.time()
print(gbdt_end-gbdt_start)


result_gnet_boosting2 <- extract_edges(gnet_results_gbdt)[rownames(true_network_matrix),colnames(true_network_matrix)]
result_gnet_kmeans2 <- extract_edges(gnet_results_kmeans)[rownames(true_network_matrix),colnames(true_network_matrix)]


mi <- build.mim(t(exp_data2))

result_genie <- GENIE3(as.matrix(exp_data2))
result_genie[is.na(result_genie)] <- 0
result_aracne <- aracne(mi, eps = 0.4)
result_minet <- mrnetb(mi)
result_clr <- clr(mi)


network_genie <- as.matrix(result_genie)[rownames(true_network_matrix),colnames(true_network_matrix)]
network_aracne <- result_aracne[rownames(true_network_matrix),colnames(true_network_matrix)]
network_minet <- result_minet[rownames(true_network_matrix),colnames(true_network_matrix)]
network_clr <- result_clr[rownames(true_network_matrix),colnames(true_network_matrix)]



roc_aracne <- roc.curve(network_aracne[true_network_matrix==1],network_aracne[true_network_matrix==0],curve=TRUE)
roc_minet <- roc.curve(network_minet[true_network_matrix==1],network_minet[true_network_matrix==0],curve=TRUE)
roc_clr <- roc.curve(network_clr[true_network_matrix==1],network_clr[true_network_matrix==0],curve=TRUE)
roc_genie <- roc.curve(network_genie[true_network_matrix==1],network_genie[true_network_matrix==0],curve=TRUE)
roc_gnet_kmeans <- roc.curve(result_gnet_kmeans2[true_network_matrix==1],result_gnet_kmeans2[true_network_matrix==0],curve=TRUE)
roc_gnet_boosting <- roc.curve(result_gnet_boosting2[true_network_matrix==1],result_gnet_boosting2[true_network_matrix==0],curve=TRUE)


o <- c(roc_aracne$auc,roc_minet$auc,roc_clr$auc,roc_genie$auc,roc_gnet_kmeans$auc,roc_gnet_boosting$auc)
names(o) <- c('ARACNe','MRNET','CLR','GENIE3','GNET2-Kmeans','GNET2-boosting')
print(o)


plot(roc_aracne$curve[,1],roc_aracne$curve[,2], type="n",xlim=c(0,1),ylim=c(0,1),xlab='FPR',ylab = 'TPR')
lines(roc_aracne$curve[,1],roc_aracne$curve[,2],col=2)
lines(roc_minet$curve[,1],roc_minet$curve[,2],col=3)
lines(roc_clr$curve[,1],roc_clr$curve[,2],col=4)
lines(roc_genie$curve[,1],roc_genie$curve[,2],col=5)
lines(roc_gnet_kmeans$curve[,1],roc_gnet_kmeans$curve[,2],col=6)
lines(roc_gnet_boosting$curve[,1],roc_gnet_boosting$curve[,2],col=7)
lines(c(0,1),c(0,1),col=1,lty = 2)

method_names = c('ARACNe','MRNET','CLR','GENIE3','GNET2-Kmeans','GNET2-boosting')
legend('bottomright',method_names,col=2:7,lty = 1,bty = "n")
