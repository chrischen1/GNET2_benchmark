library(GNET2)

library(minet)
library(GENIE3)
library(PRROC)
library(reshape2)
library(igraph)

# library(xgboost)
# library(matrixStats)
# library(dplyr)
# library(reshape2)

exp_type = 'Simulated E. coli dataset'
data_path = 'data/results_ecoli_seed2020/'
network_file = list.files(data_path,pattern = '.+.sif',full.names = T)
regulator_file = list.files(data_path,pattern = '.+Add_external.txt',full.names = T)
expression_file = list.files(data_path,pattern = '.+Add_unnormalized_dataset.txt',full.names = T)

r_seed = 121

network_data <- read.delim(network_file,  sep=' ',header = F,as.is = T)
el_ecoli <- cbind(network_data$V1,network_data$V3)
el_ecoli <- el_ecoli[grep('^bgr_.+',network_data$V1,invert = T),]
g_ecoli <- igraph::simplify(graph_from_edgelist(el_ecoli,directed = F))
adj_ecoli <- as.matrix(as_adjacency_matrix(g_ecoli))
diag(adj_ecoli) <- 0


ecoli_tf_list <- read.csv(regulator_file,skip=1,as.is = T,header = F)$V1
exp_data <- read.delim(expression_file, as.is = T,sep = '\t',row.names = 1)
exp_data <- exp_data[apply(exp_data, 1, sd) != 0,]

gene_names <- intersect(rownames(exp_data),rownames(adj_ecoli))

exprMatr <- exp_data[gene_names,]
true_network <- adj_ecoli[gene_names,gene_names]
tf_list <- intersect(ecoli_tf_list,gene_names)
true_network_reg <- true_network[tf_list,]

exprMatr <- log2(exprMatr+1)


set.seed(r_seed)
init_method = 'boosting'
gnet_results_boosting <- gnet(exprMatr,tf_list,init_method = init_method)

init_method = 'kmeans'
gnet_results_kmeans <- gnet(exprMatr,tf_list,init_method = init_method)


result_gnet_boosting2 <- extract_edges(gnet_results_boosting)[rownames(true_network_reg),colnames(true_network_reg)]
result_gnet_kmeans2 <- extract_edges(gnet_results_kmeans)[rownames(true_network_reg),colnames(true_network_reg)]


mi <- build.mim(t(exprMatr))

result_genie <- GENIE3(as.matrix(exprMatr))
result_aracne <- aracne(mi, eps = 0.1)
result_minet <- minet(mi,method = 'mrnetb')
result_clr <- clr(mi)


result_aracne_reg <- result_aracne[tf_list,]
result_minet_reg <- result_minet[tf_list,]
result_clr_reg <- result_clr[tf_list,]
result_genie_reg <- result_genie[tf_list,colnames(true_network)]


roc_gnet_kmeans <- roc.curve(result_gnet_kmeans2[true_network_reg==1],result_gnet_kmeans2[true_network_reg==0],curve=TRUE)
roc_gnet_boosting <- roc.curve(result_gnet_boosting2[true_network_reg==1],result_gnet_boosting2[true_network_reg==0],curve=TRUE)

roc_aracne <- roc.curve(result_aracne_reg[true_network_reg==1],result_aracne_reg[true_network_reg==0],curve=TRUE)
roc_minet <- roc.curve(result_minet_reg[true_network_reg==1],result_minet_reg[true_network_reg==0],curve=TRUE)
roc_clr <- roc.curve(result_clr_reg[true_network_reg==1],result_clr_reg[true_network_reg==0],curve=TRUE)
roc_genie <- roc.curve(result_genie_reg[true_network_reg==1],result_genie_reg[true_network_reg==0],curve=TRUE)

o <- c(roc_aracne$auc,roc_minet$auc,roc_clr$auc,roc_genie$auc,roc_gnet_kmeans$auc,roc_gnet_boosting$auc)
names(o) <- c('ARACNe','MRNET','CLR','GENIE3','GNET2-Kmeans','GNET2-boosting')
print(o)

plot(roc_aracne$curve[,1],roc_aracne$curve[,2], type="n",xlim=c(0,1),ylim=c(0,1),xlab='FPR',ylab = 'TPR',main=exp_type)
lines(roc_aracne$curve[,1],roc_aracne$curve[,2],col=2)
lines(roc_minet$curve[,1],roc_minet$curve[,2],col=3)
lines(roc_clr$curve[,1],roc_clr$curve[,2],col=4)
lines(roc_genie$curve[,1],roc_genie$curve[,2],col=5)
lines(roc_gnet_kmeans$curve[,1],roc_gnet_kmeans$curve[,2],col=6)
lines(roc_gnet_boosting$curve[,1],roc_gnet_boosting$curve[,2],col=7)
lines(c(0,1),c(0,1),col=1,lty = 2)

method_names = c('ARACNe','MRNET','CLR','GENIE3','GNET2-Kmeans','GNET2-boosting')
legend('bottomright',method_names,col=2:7,lty = 1,bty = "n")
