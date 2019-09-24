###################################################################################################
## This Rscript generates benchmarking evaluation of TCGA enrichment

###################################################################################################
## Loading required packages
library(Roleswitch)
library(glmnet)

###################################################################################################
######################################## FUNCTIONS ################################################
###################################################################################################

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
  mrna_name_all = matrix(data = 0, ncol = 2, nrow = nrow(mrna_base), byrow = TRUE)
  mrna_name1 = rep(0, nrow(mrna_base))
  for(i in 1 : nrow(mrna_base)){
    mrna_name1[i] = unlist(strsplit(mrna_base[i, 1], "\\|"))[1]
    mrna_name_all[i, 1] = unlist(strsplit(mrna_base[i, 1], "\\|"))[1]
    mrna_name_all[i, 2] = unlist(strsplit(mrna_base[i, 1], "\\|"))[2]
    if(mrna_name1[i] %in% mrna == 1){
      mrna_exist[i] = i
    }
  }
  mrna_use = mrna_base[mrna_exist, mrna_name]
  mrna_name1 = mrna_name1[mrna_exist]
  mrna_name_sp = mrna_name_all[mrna_exist, ]
  mrna_fullname = mrna_base[mrna_exist, 1]
  return(list(mrna_use, mrna_name1, mrna_name_sp, mrna_fullname))
}

mirna_mrna_data = function(mirna_use, mrna_use, mirna_name, mrna_name, mrna_name_sp, mrna_fullname, cutoff){
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
  mrna_name_sp1 = mrna_name_sp[mrna_sgn,]
  mrna_fullname1 = mrna_fullname[mrna_sgn]
  rownames(mirna_use1) = mirna_name1
  rownames(mrna_use1) = mrna_fullname1
  return(list(mirna_use1, mrna_use1, mirna_name1, mrna_name1, mrna_name_sp1, mrna_fullname1))
}

MMI_location = function(Vset, mirna_name, mrna_name){
  Vset_loc = c(0, nrow(Vset))
  temp1 = as.numeric(mrna_name[, 2])
  temp2 = match(Vset[, 1], temp1)
  temp3 = match(Vset[, 2], mirna_name)
  for(i in 1 : nrow(Vset)){
    if(is.na(temp2[i]) == FALSE && is.na(temp3[i]) == FALSE){
      Vset_loc[i] = length(mirna_name) * (temp2[i] - 1) + temp3[i]
    }
  }  
  return(Vset_loc)
}

mirna_mrna_loc = function(mirna_name, mrna_name, mirna, mrna, CWCS, mrna_fullname){
  mirna_loc = match(mirna_name, mirna)
  mrna_loc = match(mrna_name, mrna)
  CWCS = t(CWCS)
  CWCS_use = CWCS[mirna_loc, mrna_loc]
  rownames(CWCS_use) = mirna_name
  colnames(CWCS_use) = mrna_fullname
  return(list(mirna_loc, mrna_loc, CWCS_use))
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

miracle_enrich = function(mirna_use1, mrna_use1, mirna_name, mrna_name, CWCS, sumup, thres_num){
  mirna_use = matrix(as.numeric(mirna_use1), nrow = nrow(mirna_use1))
  mrna_use = matrix(as.numeric(mrna_use1), nrow = nrow(mrna_use1))
  pro_matrix_agg = matrix(data = 0, ncol = nrow(mrna_use), nrow = nrow(mirna_use), byrow = TRUE)
  for(m in 1:ncol(mirna_use)){
    pro_matrix_temp = (mirna_use[, m] %*% t(mrna_use[, m])) * CWCS / sumup[m]
    pro_matrix_agg = pro_matrix_agg + pro_matrix_temp
  }
  #thres_num = floor(length(which(pro_matrix_agg != 0)) * alpha)
  order_agg = order(pro_matrix_agg, decreasing = TRUE)[1 : thres_num]
  value_agg = sort(pro_matrix_agg, decreasing = TRUE)[1 : thres_num]
  mrna_list = ceiling(order_agg/nrow(CWCS))
  mirna_list = order_agg - (mrna_list - 1) * nrow(CWCS)
  mirna_unique = unique(mirna_list)
  mrna_unique = unique(mrna_list)
  mirna_total_sum = rep(0, length(mirna_unique))
  mrna_total_sum = rep(0, length(mrna_unique))
  for(k in 1:length(mirna_unique)){
    mirna_total_sum[k] = sum(value_agg[which(mirna_list == mirna_unique[k])])
  }
  for(k in 1:length(mrna_unique)){
    mrna_total_sum[k] = sum(value_agg[which(mrna_list == mrna_unique[k])])
  }
  mirna_ranksum = sort(mirna_total_sum, decreasing = TRUE) 
  mirna_rank = mirna_unique[order(mirna_total_sum, decreasing = TRUE)]
  mrna_ranksum = sort(mrna_total_sum, decreasing = TRUE) 
  mrna_rank = mrna_unique[order(mrna_total_sum, decreasing = TRUE)]
  mirna_ranklist = mirna_name[mirna_rank] 
  mrna_ranklist = mrna_name[mrna_rank] 
  return(list(mirna_ranklist, mrna_ranklist, mirna_ranksum, mrna_ranksum))
}

promise_enrich = function(mirna_use1, mrna_use1, mirna_name, mrna_name, qMRE, thres_num){
  mirna_use = matrix(as.numeric(mirna_use1), nrow = nrow(mirna_use1))
  mrna_use = matrix(as.numeric(mrna_use1), nrow = nrow(mrna_use1))
  pro_matrix_mrna = matrix(data = 0, ncol = nrow(mrna_use), nrow = nrow(mirna_use), byrow = TRUE)
  qMRE_new = t(qMRE)
  for(m in 1:ncol(mirna_use)){
    x = matrix(mrna_use[, m])
    rownames(x) = rownames(qMRE_new)
    z = matrix(mirna_use[, m])
    rownames(z) = colnames(qMRE_new)
    rs = roleswitch(x, z, qMRE_new)
    pro_matrix_temp = t(rs$p.x)
    pro_matrix_mrna = pro_matrix_mrna + pro_matrix_temp
  }
  #thres_num = floor(length(which(pro_matrix_mrna != 0)) * alpha)
  order_mrna = order(pro_matrix_mrna, decreasing = TRUE)[1 : thres_num]
  value_mrna = sort(pro_matrix_mrna, decreasing = TRUE)[1 : thres_num]
  mrna_list = ceiling(order_mrna/nrow(qMRE))
  mirna_list = order_mrna - (mrna_list - 1) * nrow(qMRE)
  mirna_unique = unique(mirna_list)
  mrna_unique = unique(mrna_list)
  mirna_total_sum = rep(0, length(mirna_unique))
  mrna_total_sum = rep(0, length(mrna_unique))
  for(k in 1:length(mirna_unique)){
    mirna_total_sum[k] = sum(value_mrna[which(mirna_list == mirna_unique[k])])
  }
  for(k in 1:length(mrna_unique)){
    mrna_total_sum[k] = sum(value_mrna[which(mrna_list == mrna_unique[k])])
  }
  mirna_ranksum = sort(mirna_total_sum, decreasing = TRUE) 
  mirna_rank = mirna_unique[order(mirna_total_sum, decreasing = TRUE)]
  mrna_ranksum = sort(mrna_total_sum, decreasing = TRUE) 
  mrna_rank = mrna_unique[order(mrna_total_sum, decreasing = TRUE)]
  mirna_ranklist = mirna_name[mirna_rank] 
  mrna_ranklist = mrna_name[mrna_rank] 
  return(list(mirna_ranklist, mrna_ranklist, mirna_ranksum, mrna_ranksum))
}

targetscan_enrich = function(CWCS, mirna_name, mrna_name, thres_num){
  #thres_num = floor(length(which(CWCS != 0)) * alpha)
  order_tar = order(CWCS, decreasing = TRUE)[1 : thres_num]
  value_tar = sort(CWCS, decreasing = TRUE)[1 : thres_num]
  mrna_list = ceiling(order_tar/nrow(CWCS))
  mirna_list = order_tar - (mrna_list - 1) * nrow(CWCS)
  mirna_unique = unique(mirna_list)
  mrna_unique = unique(mrna_list)
  mirna_total_sum = rep(0, length(mirna_unique))
  mrna_total_sum = rep(0, length(mrna_unique))
  for(k in 1:length(mirna_unique)){
    mirna_total_sum[k] = sum(value_tar[which(mirna_list == mirna_unique[k])])
  }
  for(k in 1:length(mrna_unique)){
    mrna_total_sum[k] = sum(value_tar[which(mrna_list == mrna_unique[k])])
  }
  mirna_ranksum = sort(mirna_total_sum, decreasing = TRUE) 
  mirna_rank = mirna_unique[order(mirna_total_sum, decreasing = TRUE)]
  mrna_ranksum = sort(mrna_total_sum, decreasing = TRUE) 
  mrna_rank = mrna_unique[order(mrna_total_sum, decreasing = TRUE)]
  mirna_ranklist = mirna_name[mirna_rank] 
  mrna_ranklist = mrna_name[mrna_rank] 
  return(list(mirna_ranklist, mrna_ranklist, mirna_ranksum, mrna_ranksum)) 
}

pearson_enrich = function(mirna_use1, mrna_use1, mirna_name, mrna_name, score_matrix, thres_num){
  mirna_use = matrix(as.numeric(mirna_use1),nrow = nrow(mirna_use1))
  mrna_use = matrix(as.numeric(mrna_use1),nrow = nrow(mrna_use1))
  corMat = cor(t(mirna_use), t(mrna_use))
  corMat = corMat * score_matrix
  #thres_num = floor(length(which(corMat != 0))* alpha)
  order_agg = order(corMat, decreasing = TRUE)[1 : thres_num]
  value_agg = sort(corMat, decreasing = TRUE)[1 : thres_num]
  mrna_list = ceiling(order_agg/nrow(score_matrix))
  mirna_list = order_agg - (mrna_list - 1) * nrow(score_matrix)
  mirna_unique = unique(mirna_list)
  mrna_unique = unique(mrna_list)
  mirna_total_sum = rep(0, length(mirna_unique))
  mrna_total_sum = rep(0, length(mrna_unique))
  for(k in 1:length(mirna_unique)){
    mirna_total_sum[k] = sum(value_agg[which(mirna_list == mirna_unique[k])])
  }
  for(k in 1:length(mrna_unique)){
    mrna_total_sum[k] = sum(value_agg[which(mrna_list == mrna_unique[k])])
  }
  mirna_ranksum = sort(mirna_total_sum, decreasing = TRUE) 
  mirna_rank = mirna_unique[order(mirna_total_sum, decreasing = TRUE)]
  mrna_ranksum = sort(mrna_total_sum, decreasing = TRUE) 
  mrna_rank = mrna_unique[order(mrna_total_sum, decreasing = TRUE)]
  mirna_ranklist = mirna_name[mirna_rank] 
  mrna_ranklist = mrna_name[mrna_rank] 
  return(list(mirna_ranklist, mrna_ranklist, mirna_ranksum, mrna_ranksum))
}

LASSO_enrich = function(mirna_use1, mrna_use1, mirna_name, mrna_name, score_matrix, thres_num){
  mirna_use = matrix(as.numeric(mirna_use1), nrow = nrow(mirna_use1))
  mrna_use = matrix(as.numeric(mrna_use1), nrow = nrow(mrna_use1))
  mirna_use = t(scale(t(mirna_use)))
  mrna_use = t(scale(t(mrna_use)))
  res = matrix(data = 0, ncol = nrow(mrna_use), nrow = nrow(mirna_use), byrow = TRUE)
  for(i in 1:nrow(mrna_use)) {
    aGene = glmnet(t(mirna_use), t(mrna_use)[,i], alpha = 1, nlambda = 500)$beta
    aGene = rowMeans(as.matrix(aGene))    
    res[,i] = aGene 
  }
  res = res * score_matrix
  #thres_num = floor(length(which(res != 0))* alpha)
  order_agg = order(res, decreasing = TRUE)[1 : thres_num]
  value_agg = sort(res, decreasing = TRUE)[1 : thres_num]
  mrna_list = ceiling(order_agg/nrow(score_matrix))
  mirna_list = order_agg - (mrna_list - 1) * nrow(score_matrix)
  mirna_unique = unique(mirna_list)
  mrna_unique = unique(mrna_list)
  mirna_total_sum = rep(0, length(mirna_unique))
  mrna_total_sum = rep(0, length(mrna_unique))
  for(k in 1:length(mirna_unique)){
    mirna_total_sum[k] = sum(value_agg[which(mirna_list == mirna_unique[k])])
  }
  for(k in 1:length(mrna_unique)){
    mrna_total_sum[k] = sum(value_agg[which(mrna_list == mrna_unique[k])])
  }
  mirna_ranksum = sort(mirna_total_sum, decreasing = TRUE) 
  mirna_rank = mirna_unique[order(mirna_total_sum, decreasing = TRUE)]
  mrna_ranksum = sort(mrna_total_sum, decreasing = TRUE) 
  mrna_rank = mrna_unique[order(mrna_total_sum, decreasing = TRUE)]
  mirna_ranklist = mirna_name[mirna_rank] 
  mrna_ranklist = mrna_name[mrna_rank] 
  return(list(mirna_ranklist, mrna_ranklist, mirna_ranksum, mrna_ranksum))
}

elastic_enrich = function(mirna_use1, mrna_use1, mirna_name, mrna_name, score_matrix, thres_num){
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
  #thres_num = floor(length(which(res != 0))* alpha)
  order_agg = order(res, decreasing = TRUE)[1 : thres_num]
  value_agg = sort(res, decreasing = TRUE)[1 : thres_num]
  mrna_list = ceiling(order_agg/nrow(score_matrix))
  mirna_list = order_agg - (mrna_list - 1) * nrow(score_matrix)
  mirna_unique = unique(mirna_list)
  mrna_unique = unique(mrna_list)
  mirna_total_sum = rep(0, length(mirna_unique))
  mrna_total_sum = rep(0, length(mrna_unique))
  for(k in 1:length(mirna_unique)){
    mirna_total_sum[k] = sum(value_agg[which(mirna_list == mirna_unique[k])])
  }
  for(k in 1:length(mrna_unique)){
    mrna_total_sum[k] = sum(value_agg[which(mrna_list == mrna_unique[k])])
  }
  mirna_ranksum = sort(mirna_total_sum, decreasing = TRUE) 
  mirna_rank = mirna_unique[order(mirna_total_sum, decreasing = TRUE)]
  mrna_ranksum = sort(mrna_total_sum, decreasing = TRUE) 
  mrna_rank = mrna_unique[order(mrna_total_sum, decreasing = TRUE)]
  mirna_ranklist = mirna_name[mirna_rank] 
  mrna_ranklist = mrna_name[mrna_rank] 
  return(list(mirna_ranklist, mrna_ranklist, mirna_ranksum, mrna_ranksum))
}

Zscore_enrich = function(mirna_use1, mrna_use1, mirna_name, mrna_name, score_matrix, thres_num){
  mirna_use = matrix(as.numeric(mirna_use1),nrow=nrow(mirna_use1))
  mrna_use = matrix(as.numeric(mrna_use1),nrow=nrow(mrna_use1))
  mirna_use = t(scale(t(mirna_use)))
  mrna_use = t(scale(t(mrna_use)))
  res = matrix(data=0,ncol=nrow(mrna_use),nrow=nrow(mirna_use),byrow=TRUE)
  #zAB=(xBminA-meanB)/sdB
  for(i in 1:nrow(mirna_use)){
    for(j in 1:nrow(mrna_use)){
      indexminA = which(mirna_use[i, ] == min(mirna_use[i, ]))
      xBminA = mrna_use[j, indexminA]
      res[i,j] = abs(median(xBminA))
    }
  }
  res = res * score_matrix
  #thres_num = floor(length(which(res != 0))* alpha)
  order_agg = order(res, decreasing = TRUE)[1 : thres_num]
  value_agg = sort(res, decreasing = TRUE)[1 : thres_num]
  mrna_list = ceiling(order_agg/nrow(score_matrix))
  mirna_list = order_agg - (mrna_list - 1) * nrow(score_matrix)
  mirna_unique = unique(mirna_list)
  mrna_unique = unique(mrna_list)
  mirna_total_sum = rep(0, length(mirna_unique))
  mrna_total_sum = rep(0, length(mrna_unique))
  for(k in 1:length(mirna_unique)){
    mirna_total_sum[k] = sum(value_agg[which(mirna_list == mirna_unique[k])])
  }
  for(k in 1:length(mrna_unique)){
    mrna_total_sum[k] = sum(value_agg[which(mrna_list == mrna_unique[k])])
  }
  mirna_ranksum = sort(mirna_total_sum, decreasing = TRUE) 
  mirna_rank = mirna_unique[order(mirna_total_sum, decreasing = TRUE)]
  mrna_ranksum = sort(mrna_total_sum, decreasing = TRUE) 
  mrna_rank = mrna_unique[order(mrna_total_sum, decreasing = TRUE)]
  mirna_ranklist = mirna_name[mirna_rank] 
  mrna_ranklist = mrna_name[mrna_rank] 
  return(list(mirna_ranklist, mrna_ranklist, mirna_ranksum, mrna_ranksum))
}

fc_summary = function(mrna_list, mrna_name, mrna_valid){
  mrna_num = sum(mrna_valid[, 1] %in% mrna_name)
  result_fc = matrix(data = 0, nrow = 100, ncol = 3, byrow = TRUE)
  colnames(result_fc) = c("mrna_valid1_number", "mrna_valid_fc", "mrna_valid_pvalue")
  
  a = matrix(data = 0, ncol = 2, nrow = 2, byrow = TRUE)
  for(i in 1 : 100){
    result_fc[i, 1] = sum(mrna_list[1:i] %in% mrna_valid[, 1])
    result_fc[i, 2] = (sum(mrna_list[1:i] %in% mrna_valid[, 1])/i)/(mrna_num/length(mrna_name))
    a[2, 1] = result_fc[i, 1]
    a[2, 2] = i - a[2, 1]
    a[1, 1] = mrna_num - a[2, 1]
    a[1, 2] = length(mrna_name) - i - a[1, 1]
    result_fc[i, 3] = -log10(fisher.test(a, alternative = "less")$p.value)
  }
  return(result_fc)
}

result_summary = function(name_cancer, mirna_cancer, mrna_cancer, score_matrix, alpha, method){
  score_list = matrix_transfer(score_matrix)
  x = name_func(name_cancer, mirna_cancer, mrna_cancer)
  z1 = mirna_matrix(name_cancer, mirna_cancer, x[[1]], score_list[[1]])
  z2 = mrna_matrix(name_cancer, mrna_cancer, x[[2]], score_list[[2]])
  if(method == "miracle"){
    z3 = mirna_mrna_data(z1[[1]], z2[[1]], z1[[2]], z2[[2]], z2[[3]], z2[[4]], 1)
    #z5 = MMI_location(Vset, z3[[3]], z3[[5]]) #this is used to locate the validation set
    z7 = mirna_mrna_loc(z3[[3]], z3[[4]], score_list[[1]], score_list[[2]], score_list[[3]], z3[[6]])
    z8 = sum_miracle(z3[[1]], z3[[2]], z7[[3]]) 
    z9 = miracle_enrich(z3[[1]], z3[[2]], z3[[3]], z3[[4]], z7[[3]], z8, alpha) #result of miracle
    return(list(z9, z3))
  }
  else if(method == "promise"){
    z3 = mirna_mrna_data(z1[[1]], z2[[1]], z1[[2]], z2[[2]], z2[[3]], z2[[4]], 1)
    #z5 = MMI_location(Vset, z3[[3]], z3[[5]]) #this is used to locate the validation set
    z6 = mirna_mrna_loc(z3[[3]], z3[[4]], score_list[[1]], score_list[[2]], score_list[[3]], z3[[6]])
    z10 = promise_enrich(z3[[1]], z3[[2]], z3[[3]], z3[[4]], z6[[3]], alpha) #result of promise all
    return(list(z10, z3))
  }
  else if(method == "sequence"){
    z4 = mirna_mrna_data(z1[[1]], z2[[1]], z1[[2]], z2[[2]], z2[[3]], z2[[4]], 0.8)
    #z13 = MMI_location(Vset, z4[[3]], z4[[5]]) #this is used to locate the validation set
    z7_tar = mirna_mrna_loc(z4[[3]], z4[[4]], score_list[[1]], score_list[[2]], score_list[[3]], z4[[6]])
    z12 = targetscan_enrich(z7_tar[[3]], z4[[3]], z4[[4]], alpha) #get the rankers of targetscan_qMRE
    return(list(z12, z4))
  }
  else if(method == "pearson"){
    z3 = mirna_mrna_data(z1[[1]], z2[[1]], z1[[2]], z2[[2]], z2[[3]], z2[[4]], 1)
    #z13 = MMI_location(Vset, z4[[3]], z4[[5]]) #this is used to locate the validation set
    z7 = mirna_mrna_loc(z3[[3]], z3[[4]], score_list[[1]], score_list[[2]], score_list[[3]], z3[[6]])
    z12 = pearson_enrich(z3[[1]], z3[[2]], z3[[3]], z3[[4]], z7[[3]], alpha) #get the rankers of pearson
    return(list(z12, z3))
  }
  else if(method == "lasso"){
    z3 = mirna_mrna_data(z1[[1]], z2[[1]], z1[[2]], z2[[2]], z2[[3]], z2[[4]], 1)
    #z13 = MMI_location(Vset, z4[[3]], z4[[5]]) #this is used to locate the validation set
    z7 = mirna_mrna_loc(z3[[3]], z3[[4]], score_list[[1]], score_list[[2]], score_list[[3]], z3[[6]])
    z12 = LASSO_enrich(z3[[1]], z3[[2]], z3[[3]], z3[[4]], z7[[3]], alpha) #get the rankers of lasso
    return(list(z12, z3))
  }
  else if(method == "elastic"){
    z3 = mirna_mrna_data(z1[[1]], z2[[1]], z1[[2]], z2[[2]], z2[[3]], z2[[4]], 1)
    #z13 = MMI_location(Vset, z4[[3]], z4[[5]]) #this is used to locate the validation set
    z7 = mirna_mrna_loc(z3[[3]], z3[[4]], score_list[[1]], score_list[[2]], score_list[[3]], z3[[6]])
    z12 = elastic_enrich(z3[[1]], z3[[2]], z3[[3]], z3[[4]], z7[[3]], alpha) #get the rankers of elastic net
    return(list(z12, z3))
  }
  else if(method == "Zscore"){
    z3 = mirna_mrna_data(z1[[1]], z2[[1]], z1[[2]], z2[[2]], z2[[3]], z2[[4]], 1)
    #z13 = MMI_location(Vset, z4[[3]], z4[[5]]) #this is used to locate the validation set
    z7 = mirna_mrna_loc(z3[[3]], z3[[4]], score_list[[1]], score_list[[2]], score_list[[3]], z3[[6]])
    z12 = Zscore_enrich(z3[[1]], z3[[2]], z3[[3]], z3[[4]], z7[[3]], alpha) #get the rankers of Zscore
    return(list(z12, z3))
  }
}

###################################################################################################
####################################### DATA INPUTS ###############################################
###################################################################################################

CWCS_all = as.matrix(read.table("TargetScan7_CWCS.txt", head=TRUE, sep = "\t"))
CWCS_conserved = as.matrix(read.table("TargetScan7_CWCS_cons.txt", head=TRUE, sep = "\t"))
qMRE_conserved = as.matrix(read.table("TargetScan7_qMRE_cons.txt", head = TRUE, sep = "\t"))
microT = as.matrix(read.table("DIANA_microT_CDS.txt", head = TRUE, sep = "\t"))
combine = as.matrix(read.table("Combine_MMIs.txt", head = TRUE, sep = "\t"))

mrna_valid = as.matrix(read.table("Cancer_gene_set.txt", head = TRUE, sep = "\t"))

name_cancer = as.matrix(read.table("TCGA_HNSC_sampleMatch.txt", head = TRUE, sep = "\t"))
mirna_cancer = as.matrix(read.table("HNSC_miRNA_expression.txt", head = FALSE, sep = "\t"))
mrna_cancer = as.matrix(read.table("HNSC_mRNA_expression.txt", head = FALSE, sep = "\t"))

###################################################################################################
######################################### MAIN code ###############################################
###################################################################################################

#the following code is used for preparation
thres_num = 10000

z9 = result_summary(name_cancer, mirna_cancer, mrna_cancer, CWCS_all, thres_num, "miracle")
z9_con = result_summary(name_cancer, mirna_cancer, mrna_cancer, CWCS_conserved, thres_num, "miracle")
z10 = result_summary(name_cancer, mirna_cancer, mrna_cancer, qMRE_conserved, thres_num, "promise")
z11 = result_summary(name_cancer, mirna_cancer, mrna_cancer, CWCS_conserved, thres_num, "sequence")
z12 = result_summary(name_cancer, mirna_cancer, mrna_cancer, microT, thres_num, "sequence")
z13 = result_summary(name_cancer, mirna_cancer, mrna_cancer, combine, thres_num, "pearson")
z14 = result_summary(name_cancer, mirna_cancer, mrna_cancer, combine, thres_num, "lasso")
z15 = result_summary(name_cancer, mirna_cancer, mrna_cancer, combine, thres_num, "elastic")
z16 = result_summary(name_cancer, mirna_cancer, mrna_cancer, combine, thres_num, "Zscore")

z20 = fc_summary(z9[[1]][[2]], z9[[2]][[4]], mrna_valid)
z21 = fc_summary(z9_con[[1]][[2]], z9_con[[2]][[4]], mrna_valid)
z22 = fc_summary(z10[[1]][[2]], z10[[2]][[4]], mrna_valid)
z23 = fc_summary(z11[[1]][[2]], z11[[2]][[4]], mrna_valid)
z24 = fc_summary(z12[[1]][[2]], z12[[2]][[4]], mrna_valid)
z25 = fc_summary(z13[[1]][[2]], z13[[2]][[4]], mrna_valid)
z26 = fc_summary(z14[[1]][[2]], z14[[2]][[4]], mrna_valid)
z27 = fc_summary(z15[[1]][[2]], z15[[2]][[4]], mrna_valid)
z28 = fc_summary(z16[[1]][[2]], z16[[2]][[4]], mrna_valid)

Cancer_gene_pvalue = cbind(z20[, 3], z21[, 3], z22[, 3], z23[, 3], z24[, 3], z25[, 3], z26[, 3], z27[, 3], z28[, 3])
Cancer_gene_fc = cbind(z20[, 2], z21[, 2], z22[, 2], z23[, 2], z24[, 2], z25[, 2], z26[, 2], z27[, 2], z28[, 2])

colnames(Cancer_gene_pvalue) = c('miRACLe', 'miRACLe_cons', 'ProMISe_mrna_cons', 'TargetScan7_CWCS_cons', 'DIANA-microT-CDS', 'PearsonCorrelation', 
                             'LASSO', 'Elastic-net', 'Zscore')
colnames(Cancer_gene_fc) = c('miRACLe', 'miRACLe_cons', 'ProMISe_mrna_cons', 'TargetScan7_CWCS_cons', 'DIANA-microT-CDS', 'PearsonCorrelation', 
                             'LASSO', 'Elastic-net', 'Zscore')

write.table(z20, file = "Enrichment for miRACLe.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(z21, file = "Enrichment for miRACLe_cons.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(z22, file = "Enrichment for ProMISe_mrna_cons.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(z23, file = "Enrichment for TargetScan_CWCS_cons.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(z24, file = "Enrichment for DIANA_microT_CDS.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(z25, file = "Enrichment for PearsonCorrelation.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(z26, file = "Enrichment for LASSO.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(z27, file = "Enrichment for Elastic_net.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(z28, file = "Enrichment for Zscore.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(Cancer_gene_pvalue, file = "Cancer_gene_pvalue.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(Cancer_gene_fc, file = "Cancer_gene_fc.txt", quote = FALSE, sep = "\t", row.names = FALSE)



