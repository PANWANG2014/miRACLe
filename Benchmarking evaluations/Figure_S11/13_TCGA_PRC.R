###################################################################################################
## This Rscript generates benchmarking evaluation of TCGA - precision-recall curve

###################################################################################################
## Loading required packages
library(Roleswitch)
library(glmnet)
library(ROCR)

###################################################################################################
######################################## FUNCTIONS ################################################
###################################################################################################

valid_trans = function(vset, trans_ID){
  vsett = vset[which(as.numeric(vset[, 1]) %in% as.numeric(trans_ID[, 2])),]
  vsett[, 1] = trans_ID[match(as.numeric(vsett[, 1]), as.numeric(trans_ID[, 2])), 1]
  return(vsett)
}

matrix_transfer = function(input_matrix){
  mirna_name = unique(input_matrix[, 2])
  mrna_name = unique(input_matrix[, 1])
  mod_matrix = matrix(data = 0, ncol = length(mirna_name), nrow = length(mrna_name))
  pos1 = match(input_matrix[, 2], mirna_name)
  pos2 = match(input_matrix[, 1], mrna_name)
  for(i in 1:nrow(input_matrix)){
    mod_matrix[pos2[i], pos1[i]] = mod_matrix[pos2[i], pos1[i]] + abs(as.numeric(input_matrix[i, 3]))
  }
  return(list(mirna_name, mrna_name, mod_matrix))
}

name_func = function(name_base, mirna_base, mrna_base){
  x = match(name_base[, "miRNA"], mirna_base[1, ])
  y = match(name_base[, "mRNA"], mrna_base[1, ])
  return(list(x,y))
}

mirna_matrix = function(name_base, mirna_base, mirna_name, mirna){  
  mirna_name1 = rep(0,nrow(mirna_base))
  for(i in 1:nrow(mirna_base)){
    mirna_name1[i] = unlist(strsplit(mirna_base[i,1],"\\|"))[1]
  }
  mirna_use = mirna_base[mirna_name1%in%mirna, mirna_name]
  mirna_name1 = mirna_name1[mirna_name1%in%mirna]
  return(list(mirna_use, mirna_name1))
}

mrna_matrix = function(name_base, mrna_base, mrna_name, mrna){
  mrna_exist = rep(0,nrow(mrna_base))
  mrna_name1 = rep(0, nrow(mrna_base))
  for(i in 1 : nrow(mrna_base)){
    mrna_name1[i] = unlist(strsplit(mrna_base[i, 1], "\\|"))[1]
    if(mrna_name1[i] %in% mrna == 1){
      mrna_exist[i] = i
    }
  }
  mrna_use = mrna_base[mrna_exist, mrna_name]
  mrna_name1 = mrna_name1[mrna_exist]
  return(list(mrna_use, mrna_name1))
}

mirna_mrna_data = function(mirna_use, mrna_use, mirna_name, mrna_name, cutoff){
  mirna_use[is.na(mirna_use)] = 0
  mrna_use[is.na(mrna_use)] = 0
  mirna_use = matrix(as.numeric(mirna_use), nrow = nrow(mirna_use))
  mrna_use = matrix(as.numeric(mrna_use), nrow = nrow(mrna_use))
  mirna_use[mirna_use < 0] = 0
  mrna_use[mrna_use < 0] = 0
  mirna_sgn = seq(1, nrow(mirna_use), 1)
  mrna_sgn = seq(1, nrow(mrna_use), 1)
  for (i in 1:nrow(mirna_use)){
    if(sum(mirna_use[i,] == 0) >= ncol(mirna_use) * cutoff){
      mirna_sgn[i] = 0
    }
  }
  for (i in 1:nrow(mrna_use)){
    if(sum(mrna_use[i,] == 0) >= ncol(mrna_use) * cutoff){
      mrna_sgn[i] = 0
    }
  }  
  mirna_use1 = mirna_use[mirna_sgn,]
  mrna_use1 = mrna_use[mrna_sgn,]
  mirna_name1 = mirna_name[mirna_sgn]
  mrna_name1 = mrna_name[mrna_sgn]
  rownames(mirna_use1) = mirna_name1
  rownames(mrna_use1) = mrna_name1
  return(list(mirna_use1, mrna_use1, mirna_name1, mrna_name1))
}

mirna_mrna_loc = function(mirna_name, mrna_name, mirna, mrna, CWCS){
  mirna_loc = match(mirna_name, mirna)
  mrna_loc = match(mrna_name, mrna)
  CWCS = t(CWCS)
  CWCS_use = CWCS[mirna_loc, mrna_loc]
  rownames(CWCS_use) = mirna_name
  colnames(CWCS_use) = mrna_name
  return(list(mirna_name, mrna_name, CWCS_use))
}

sum_miracle = function(mirna, mrna, CWCS){
  res = rep(0, ncol(mirna))
  for(i in 1:ncol(mirna)){
    mirna_use1 = as.numeric(mirna[, i])
    mrna_use1 = as.numeric(mrna[, i])
    temp = mirna_use1 %*% t(mrna_use1)
    temp2 = CWCS*temp
    res[i] = sum(temp2)
  }
  return(res)
}

miracle_pop = function(mirna_use1, mrna_use1, mirna_name, mrna_name, CWCS, sumup){
  mirna_use = matrix(as.numeric(mirna_use1), nrow = nrow(mirna_use1))
  mrna_use = matrix(as.numeric(mrna_use1), nrow = nrow(mrna_use1))
  pro_matrix_temp = matrix(data = 0, ncol = nrow(mrna_use), nrow = nrow(mirna_use), byrow = TRUE)
  pro_matrix_agg = matrix(data = 0, ncol = nrow(mrna_use), nrow = nrow(mirna_use), byrow = TRUE)
  for(m in 1:ncol(mirna_use)){
    pro_matrix_temp = (mirna_use[, m] %*% t(mrna_use[, m])) * CWCS / sumup[m]
    pro_matrix_agg = pro_matrix_agg + pro_matrix_temp
  }
  return(list(mirna_name, mrna_name, pro_matrix_agg))
}

promise_pop = function(mirna_use1, mrna_use1, mirna_name, mrna_name, CWCS){
  mirna_use = matrix(as.numeric(mirna_use1), nrow = nrow(mirna_use1))
  mrna_use = matrix(as.numeric(mrna_use1), nrow = nrow(mrna_use1))
  pro_matrix_mrna = matrix(data = 0, ncol = nrow(mrna_use), nrow = nrow(mirna_use), byrow = TRUE)
  CWCS_new = t(CWCS)
  for(m in 1:ncol(mirna_use)){
    x = matrix(mrna_use[, m])
    rownames(x) = rownames(CWCS_new)
    z = matrix(mirna_use[, m])
    rownames(z) = colnames(CWCS_new)
    rs = roleswitch(x, z, CWCS_new)
    pro_matrix_temp1 = t(rs$p.x)
    pro_matrix_mrna = pro_matrix_mrna + pro_matrix_temp1
  }
  return(list(mirna_name, mrna_name, pro_matrix_mrna))
}

pearson_cal = function(mirna_use1, mrna_use1, mirna_name, mrna_name, score_matrix){
  mirna_use = matrix(as.numeric(mirna_use1),nrow = nrow(mirna_use1))
  mrna_use = matrix(as.numeric(mrna_use1),nrow = nrow(mrna_use1))
  corMat = cor(t(mirna_use), t(mrna_use))
  corMat = corMat * score_matrix
  return(list(mirna_name, mrna_name, corMat))
}

LASSO_cal = function(mirna_use1, mrna_use1, mirna_name, mrna_name, score_matrix){
  mirna_use = matrix(as.numeric(mirna_use1), nrow = nrow(mirna_use1))
  mrna_use = matrix(as.numeric(mrna_use1), nrow = nrow(mrna_use1))
  mirna_use = t(scale(t(mirna_use)))
  mrna_use = t(scale(t(mrna_use)))
  
  res = matrix(data = 0, ncol = nrow(mrna_use), nrow = nrow(mirna_use), byrow = TRUE)
  for(i in 1:nrow(mrna_use)) {
    aGene = glmnet(t(mirna_use), t(mrna_use)[,i], alpha = 1)$beta
    aGene = rowMeans(as.matrix(aGene))    
    res[,i] = aGene 
  }
  res = res * score_matrix
  return(list(mirna_name, mrna_name, res))
}

ELASTIC_cal = function(mirna_use1, mrna_use1, mirna_name, mrna_name, score_matrix){ 
  mirna_use = matrix(as.numeric(mirna_use1), nrow = nrow(mirna_use1))
  mrna_use = matrix(as.numeric(mrna_use1), nrow = nrow(mrna_use1))
  mirna_use = t(scale(t(mirna_use)))
  mrna_use = t(scale(t(mrna_use)))
  
  res = matrix(data = 0, ncol = nrow(mrna_use), nrow = nrow(mirna_use), byrow = TRUE)
  for(i in 1:nrow(mrna_use)) {
    aGene = glmnet(t(mirna_use), t(mrna_use)[,i], alpha = 0.5)$beta
    aGene = rowMeans(as.matrix(aGene))     
    res[,i] = aGene 
  }
  res = res * score_matrix

  return(list(mirna_name, mrna_name, res))
}

Zscore_cal = function(mirna_use1, mrna_use1, mirna_name, mrna_name, score_matrix){
  mirna_use = matrix(as.numeric(mirna_use1),nrow=nrow(mirna_use1))
  mrna_use = matrix(as.numeric(mrna_use1),nrow=nrow(mrna_use1))
  mirna_use = t(scale(t(mirna_use)))
  mrna_use = t(scale(t(mrna_use)))
  
  res = matrix(data=0,ncol=nrow(mrna_use),nrow=nrow(mirna_use),byrow=TRUE)
  for (i in 1:nrow(mirna_use)){
    for (j in 1:nrow(mrna_use)){
      indexminA = which(mirna_use[i, ] == min(mirna_use[i, ]))
      xBminA = mrna_use[j, indexminA]
      res[i,j] = abs(median(xBminA))
    }
  }
  res = res * score_matrix
  return(list(mirna_name, mrna_name, res))
}

result_summary = function(name_cancer, mirna_cancer, mrna_cancer, score_matrix, method){
  score_list = matrix_transfer(score_matrix)
  x = name_func(name_cancer, mirna_cancer, mrna_cancer)
  z1 = mirna_matrix(name_cancer, mirna_cancer, x[[1]], score_list[[1]])
  z2 = mrna_matrix(name_cancer, mrna_cancer, x[[2]], score_list[[2]])
  z3 = mirna_mrna_data(z1[[1]], z2[[1]], z1[[2]], z2[[2]], 1)
  z7 = mirna_mrna_loc(z3[[3]], z3[[4]], score_list[[1]], score_list[[2]], score_list[[3]])
  if(method == "miracle"){
    z8 = sum_miracle(z3[[1]], z3[[2]], z7[[3]])
    z9 = miracle_pop(z3[[1]], z3[[2]], z3[[3]], z3[[4]], z7[[3]], z8) #result of miracle
    return(z9)
  }
  else if(method == "promise"){
    z10 = promise_pop(z3[[1]], z3[[2]], z3[[3]], z3[[4]], z7[[3]]) #result of promise
    return(z10)
  }
  else if(method == "target"){
    z4 = mirna_mrna_data(z1[[1]], z2[[1]], z1[[2]], z2[[2]], 0.8)
    z6 = mirna_mrna_loc(z4[[3]], z4[[4]], score_list[[1]], score_list[[2]], score_list[[3]])
    return(z6)
  }
  else if(method == "microT"){
    return(z7)
  }
  else if(method == "pearson"){
    z11 = pearson_cal(z3[[1]], z3[[2]], z3[[3]], z3[[4]], z7[[3]])
    return(z11)
  }
  else if(method == "lasso"){
    z12 = LASSO_cal(z3[[1]], z3[[2]], z3[[3]], z3[[4]], z7[[3]]) #result of LASSO
    return(z12)
  }
  else if(method == "elastic"){
    z13 = ELASTIC_cal(z3[[1]], z3[[2]], z3[[3]], z3[[4]], z7[[3]]) #result of elastic net
    return(z13)
  }
  else if(method == "Zscore"){
    z14 = Zscore_cal(z3[[1]], z3[[2]], z3[[3]], z3[[4]], z7[[3]]) #result of Zscore score
    return(z14)
  }
}

vali_pos = function(matrix, pos_list){
  miRNA_list = matrix[[1]]
  mRNA_list = matrix[[2]]
  score_matrix = matrix[[3]]
  temp1_pos = match(pos_list[, 1], mRNA_list)
  temp2_pos = match(pos_list[, 2], miRNA_list)
  temp3_pos = intersect(which(is.na(temp1_pos) == FALSE), which(is.na(temp2_pos) == FALSE))
  pos_list_mod = pos_list[temp3_pos, ]
  temp4_pos = match(pos_list_mod[, 1], mRNA_list)
  temp5_pos = match(pos_list_mod[, 2], miRNA_list)
  positive_position = (temp4_pos - 1) * length(miRNA_list) + temp5_pos
  
  pos = which(score_matrix > 0)
  poss = positive_position[which(positive_position %in% pos)]
  labels = rep(0, length(pos))
  labels[which(pos %in% poss)] = 1
  predictions = score_matrix[pos]
  simple = list(predictions, labels)
  names(simple) = c("predictions", "labels")
  df = data.frame(simple)
  pred = prediction(df$predictions, df$labels)
  return(pred)
}

###################################################################################################
####################################### DATA INPUTS ###############################################
###################################################################################################

symbol_to_ID = as.matrix(read.table("Symbol to ID.txt", head = TRUE, sep = "\t"))

CWCS_all = as.matrix(read.table("TargetScan7_CWCS.txt", head=TRUE, sep = "\t"))
CWCS_conserved = as.matrix(read.table("TargetScan7_CWCS_cons.txt", head=TRUE, sep = "\t"))
qMRE_conserved = as.matrix(read.table("TargetScan7_qMRE_cons.txt", head = TRUE, sep = "\t"))
microT = as.matrix(read.table("DIANA_microT_CDS.txt", head = TRUE, sep = "\t"))
combine = as.matrix(read.table("Combine_MMIs.txt", head = TRUE, sep = "\t"))

Vset1 = as.matrix(read.table("Vset_all.txt", head = TRUE, sep = "\t"))

name_cancer = as.matrix(read.table("TCGA_BRCA_sampleMatch.txt", head = TRUE, sep = "\t"))
mirna_cancer = as.matrix(read.table("BRCA_miRNA_expression.txt", head = FALSE, sep = "\t"))
mrna_cancer = as.matrix(read.table("BRCA_mRNA_expression.txt", head = FALSE, sep = "\t"))

###################################################################################################
######################################### MAIN code ###############################################
###################################################################################################

V1_trans = valid_trans(Vset1, symbol_to_ID)

miracle_all = result_summary(name_cancer, mirna_cancer, mrna_cancer, CWCS_all, "miracle")
miracle_cons = result_summary(name_cancer, mirna_cancer, mrna_cancer, CWCS_conserved, "miracle")
promise_cons = result_summary(name_cancer, mirna_cancer, mrna_cancer, qMRE_conserved, "promise")
targetscan = result_summary(name_cancer, mirna_cancer, mrna_cancer, CWCS_conserved, "target")
microT = result_summary(name_cancer, mirna_cancer, mrna_cancer, microT, "microT")
pearson = result_summary(name_cancer, mirna_cancer, mrna_cancer, combine, "pearson")
lasso = result_summary(name_cancer, mirna_cancer, mrna_cancer, combine, "lasso")
elastic = result_summary(name_cancer, mirna_cancer, mrna_cancer, combine, "elastic")
Zscore = result_summary(name_cancer, mirna_cancer, mrna_cancer, combine, "Zscore")

miracle_all_vali = vali_pos(miracle_all, V1_trans)
miracle_cons_vali = vali_pos(miracle_cons, V1_trans)
promise_cons_vali = vali_pos(promise_cons, V1_trans)
targetscan_vali = vali_pos(targetscan, V1_trans)
microT_vali = vali_pos(microT, V1_trans)
pearson_vali = vali_pos(pearson, V1_trans)
lasso_vali = vali_pos(lasso, V1_trans)
elastic_vali = vali_pos(elastic, V1_trans)
Zscore_vali = vali_pos(Zscore, V1_trans)

perf_miracle_all = performance(miracle_all_vali, "prec", "rec")
perf_miracle_cons = performance(miracle_cons_vali, "prec", "rec")
perf_promise_cons = performance(promise_cons_vali, "prec", "rec")
perf_targetscan = performance(targetscan_vali, "prec", "rec")
perf_microT = performance(microT_vali, "prec", "rec")
perf_pearson = performance(pearson_vali, "prec", "rec")
perf_lasso = performance(lasso_vali, "prec", "rec")
perf_elastic = performance(elastic_vali, "prec", "rec")
perf_Zscore = performance(Zscore_vali, "prec", "rec")

pdf(file = "PR_HNSC_Vall.pdf", height = 5, width = 5)
plot(perf_miracle_cons, col = rgb(r=228, g=26, b=28, maxColorValue = 255), cex.lab = 1, ylim = c(0, 1),
     lwd = 1, main = "Precision-Recall plot for HNSC in Vall", cex.main = 1, cex.axis = 1)
plot(perf_miracle_all, add = TRUE, col = rgb(r=55, g=126, b=184, maxColorValue = 255), lwd = 1)
plot(perf_promise_cons, add = TRUE, col = rgb(r=0, g=68, b=27, maxColorValue = 255), lwd = 1)
plot(perf_targetscan, add = TRUE, col = rgb(r=152, g=78, b=163, maxColorValue = 255), lwd = 1)
plot(perf_microT, add = TRUE, col = rgb(r=53, g=151, b=143, maxColorValue = 255), lwd = 1)
plot(perf_pearson, add = TRUE, col = rgb(r=255, g=127, b=0, maxColorValue = 255), lwd = 1)
plot(perf_lasso, add = TRUE, col = rgb(r=153, g=153, b=153, maxColorValue = 255), lwd = 1)
plot(perf_elastic, add = TRUE, col = rgb(r=166, g=86, b=40, maxColorValue = 255), lwd = 1)
plot(perf_Zscore, add = TRUE, col = rgb(r=247, g=129, b=191, maxColorValue = 255), lwd = 1)
dev.off()

