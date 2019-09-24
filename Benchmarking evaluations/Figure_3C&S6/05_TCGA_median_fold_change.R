###################################################################################################
## This Rscript generates benchmarking evaluation of TCGA - median fold change

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

miracle_medfc = function(mirna_use1, mrna_use1, mrna_name, CWCS, sumup, mirna_name, vset, num){
  mirna_use = matrix(as.numeric(mirna_use1), nrow = nrow(mirna_use1))
  mrna_use = matrix(as.numeric(mrna_use1), nrow = nrow(mrna_use1))
  pro_matrix_agg = matrix(data = 0, ncol = nrow(mrna_use), nrow = nrow(mirna_use))
  for(m in 1:ncol(mirna_use)){
    pro_matrix_temp = (mirna_use[, m] %*% t(mrna_use[, m])) * CWCS / sumup[m]
    pro_matrix_agg = pro_matrix_agg + pro_matrix_temp
  }
  #thres_num = ceiling(length(which(pro_matrix_agg != 0))/100)
  vset_name = as.vector(vset[1,-1])
  temp = match(vset_name, mirna_name)
  mirna_pos = temp[which(is.na(temp) != TRUE)]
  pro_matrix_sel = pro_matrix_agg[mirna_pos, ]
  sel_mrna_name = matrix(data = 0, ncol = length(mirna_pos), nrow = num, byrow = TRUE)
  sel_mrna_score = matrix(data = 0, ncol = length(mirna_pos), nrow = num, byrow = TRUE)
  tar_num = rep(0, length(mirna_pos))
  for(i in 1: length(mirna_pos)){
    sel_mrna_score[,i] = sort(pro_matrix_sel[i,], decreasing = TRUE)[1:num]
    temp = order(pro_matrix_sel[i,], decreasing = TRUE)[1:num]
    sel_mrna_name[,i] = as.numeric(mrna_name[temp,2])
    tar_num[i] = length(which(pro_matrix_sel[i,] != 0))
  }
  colnames(sel_mrna_name) = mirna_name[mirna_pos]
  colnames(sel_mrna_score) = mirna_name[mirna_pos]
  return(list(sel_mrna_name, sel_mrna_score, tar_num, pro_matrix_agg))
}

promise_medfc = function(mirna_use1, mrna_use1, mrna_name, CWCS, mirna_name, vset, num){
  mirna_use = matrix(as.numeric(mirna_use1), nrow = nrow(mirna_use1))
  mrna_use = matrix(as.numeric(mrna_use1), nrow = nrow(mrna_use1))
  pro_matrix_mrna = matrix(data = 0, ncol = nrow(mrna_use), nrow = nrow(mirna_use), byrow = TRUE)
  pro_matrix_full = matrix(data = 0, ncol = nrow(mrna_use), nrow = nrow(mirna_use), byrow = TRUE)
  CWCS_new = t(CWCS)
  for(m in 1:ncol(mirna_use)){
    x = matrix(mrna_use[, m])
    rownames(x) = rownames(CWCS_new)
    z = matrix(mirna_use[, m])
    rownames(z) = colnames(CWCS_new)
    rs = roleswitch(x, z, CWCS_new)
    pro_matrix_temp1 = t(rs$p.x)
    pro_matrix_temp2 = t(rs$p.xz)
    pro_matrix_mrna = pro_matrix_mrna + pro_matrix_temp1
    pro_matrix_full = pro_matrix_full + pro_matrix_temp2
  }
  #thres_num1 = ceiling(length(which(pro_matrix_mrna != 0))/100)
  #thres_num2 = ceiling(length(which(pro_matrix_full != 0))/100)
  vset_name = as.vector(vset[1,-1])
  temp = match(vset_name, mirna_name)
  mirna_pos = temp[which(is.na(temp) != TRUE)]
  sel_matrix_mrna = pro_matrix_mrna[mirna_pos, ]
  sel_matrix_full = pro_matrix_full[mirna_pos, ]
  sel_name_mrna = matrix(data = 0, ncol = length(mirna_pos), nrow = num, byrow = TRUE)
  sel_score_mrna = matrix(data = 0, ncol = length(mirna_pos), nrow = num, byrow = TRUE)
  sel_name_full = matrix(data = 0, ncol = length(mirna_pos), nrow = num, byrow = TRUE)
  sel_score_full = matrix(data = 0, ncol = length(mirna_pos), nrow = num, byrow = TRUE)
  tar_num_mrna = rep(0, length(mirna_pos))
  tar_num_full = rep(0, length(mirna_pos))
  for(i in 1: length(mirna_pos)){
    sel_score_mrna[,i] = sort(sel_matrix_mrna[i,], decreasing = TRUE)[1:num]
    temp1 = order(sel_matrix_mrna[i,], decreasing = TRUE)[1:num]
    sel_name_mrna[,i] = as.numeric(mrna_name[temp1,2])
    tar_num_mrna[i] = length(which(sel_matrix_mrna[i,] != 0))
    sel_score_full[,i] = sort(sel_matrix_full[i,], decreasing = TRUE)[1:num]
    temp2 = order(sel_matrix_full[i,], decreasing = TRUE)[1:num]
    sel_name_full[,i] = as.numeric(mrna_name[temp2,2])
    tar_num_full[i] = length(which(sel_matrix_full[i,] != 0))
  }
  return(list(sel_name_mrna, sel_score_mrna, sel_name_full, sel_score_full, tar_num_mrna, tar_num_full))
}

pearson_medfc = function(mirna_use1, mrna_use1, mrna_name, mirna_name, score_matrix, vset, num){
  mirna_use = matrix(as.numeric(mirna_use1),nrow = nrow(mirna_use1))
  mrna_use = matrix(as.numeric(mrna_use1),nrow = nrow(mrna_use1))
  corMat = cor(t(mirna_use), t(mrna_use))
  corMat = corMat * score_matrix
  vset_name = as.vector(vset[1,-1])
  temp = match(vset_name, mirna_name)
  mirna_pos = temp[which(is.na(temp) != TRUE)]
  pro_matrix_sel = corMat[mirna_pos, ]
  sel_mrna_name = matrix(data = 0, ncol = length(mirna_pos), nrow = num, byrow = TRUE)
  sel_mrna_score = matrix(data = 0, ncol = length(mirna_pos), nrow = num, byrow = TRUE)
  tar_num = rep(0, length(mirna_pos))
  for(i in 1: length(mirna_pos)){
    sel_mrna_score[,i] = sort(pro_matrix_sel[i,], decreasing = TRUE)[1:num]
    temp = order(pro_matrix_sel[i,], decreasing = TRUE)[1:num]
    sel_mrna_name[,i] = as.numeric(mrna_name[temp,2])
    tar_num[i] = length(which(pro_matrix_sel[i,] != 0))
  }
  return(list(sel_mrna_name, sel_mrna_score, tar_num))
}

LASSO_medfc = function(mirna_use1, mrna_use1, mrna_name, mirna_name, score_matrix, vset, num){
  mirna_use = matrix(as.numeric(mirna_use1), nrow = nrow(mirna_use1))
  mrna_use = matrix(as.numeric(mrna_use1), nrow = nrow(mrna_use1))
  res = matrix(data = 0, ncol = nrow(mrna_use), nrow = nrow(mirna_use), byrow = TRUE)
  for(i in 1:nrow(mrna_use)) {
    aGene = glmnet(t(mirna_use), t(mrna_use)[,i], alpha = 1, nlambda = 500)$beta
    aGene = rowMeans(as.matrix(aGene))    
    res[,i] = aGene 
  }
  res = res * score_matrix
  vset_name = as.vector(vset[1,-1])
  temp = match(vset_name, mirna_name)
  mirna_pos = temp[which(is.na(temp) != TRUE)]
  pro_matrix_sel = res[mirna_pos, ]
  sel_mrna_name = matrix(data = 0, ncol = length(mirna_pos), nrow = num, byrow = TRUE)
  sel_mrna_score = matrix(data = 0, ncol = length(mirna_pos), nrow = num, byrow = TRUE)
  tar_num = rep(0, length(mirna_pos))
  for(i in 1: length(mirna_pos)){
    sel_mrna_score[,i] = sort(pro_matrix_sel[i,], decreasing = TRUE)[1:num]
    temp = order(pro_matrix_sel[i,], decreasing = TRUE)[1:num]
    sel_mrna_name[,i] = as.numeric(mrna_name[temp,2])
    tar_num[i] = length(which(pro_matrix_sel[i,] != 0))
  }
  
  return(list(sel_mrna_name, sel_mrna_score, tar_num))
}

ELASTIC_medfc = function(mirna_use1, mrna_use1, mrna_name, mirna_name, score_matrix, vset, num){ 
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
  vset_name = as.vector(vset[1,-1])
  temp = match(vset_name, mirna_name)
  mirna_pos = temp[which(is.na(temp) != TRUE)]
  pro_matrix_sel = res[mirna_pos, ]
  sel_mrna_name = matrix(data = 0, ncol = length(mirna_pos), nrow = num, byrow = TRUE)
  sel_mrna_score = matrix(data = 0, ncol = length(mirna_pos), nrow = num, byrow = TRUE)
  tar_num = rep(0, length(mirna_pos))
  for(i in 1: length(mirna_pos)){
    sel_mrna_score[,i] = sort(pro_matrix_sel[i,], decreasing = TRUE)[1:num]
    temp = order(pro_matrix_sel[i,], decreasing = TRUE)[1:num]
    sel_mrna_name[,i] = as.numeric(mrna_name[temp,2])
    tar_num[i] = length(which(pro_matrix_sel[i,] != 0))
  }
  
  return(list(sel_mrna_name, sel_mrna_score, tar_num))
}

Zscore_medfc = function(mirna_use1, mrna_use1, mrna_name, mirna_name, score_matrix, vset, num){
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
  vset_name = as.vector(vset[1,-1])
  temp = match(vset_name, mirna_name)
  mirna_pos = temp[which(is.na(temp) != TRUE)]
  pro_matrix_sel = res[mirna_pos, ]
  sel_mrna_name = matrix(data = 0, ncol = length(mirna_pos), nrow = num, byrow = TRUE)
  sel_mrna_score = matrix(data = 0, ncol = length(mirna_pos), nrow = num, byrow = TRUE)
  tar_num = rep(0, length(mirna_pos))
  for(i in 1: length(mirna_pos)){
    sel_mrna_score[,i] = sort(pro_matrix_sel[i,], decreasing = TRUE)[1:num]
    temp = order(pro_matrix_sel[i,], decreasing = TRUE)[1:num]
    sel_mrna_name[,i] = as.numeric(mrna_name[temp,2])
    tar_num[i] = length(which(pro_matrix_sel[i,] != 0))
  }
  
  return(list(sel_mrna_name, sel_mrna_score, tar_num))
}

targetscan_medfc = function(CWCS, mrna_name, mirna_name, vset, num){
  vset_name = as.vector(vset[1,-1])
  temp = match(vset_name, mirna_name)
  mirna_pos = temp[which(is.na(temp) != TRUE)]
  pro_matrix_sel = CWCS[mirna_pos, ]
  sel_mrna_name = matrix(data = 0, ncol = length(mirna_pos), nrow = num, byrow = TRUE)
  sel_mrna_score = matrix(data = 0, ncol = length(mirna_pos), nrow = num, byrow = TRUE)
  tar_num = rep(0, length(mirna_pos))
  for(i in 1: length(mirna_pos)){
    sel_mrna_score[,i] = sort(pro_matrix_sel[i,], decreasing = TRUE)[1:num]
    temp = order(pro_matrix_sel[i,], decreasing = TRUE)[1:num]
    sel_mrna_name[,i] = as.numeric(mrna_name[temp,2])
    tar_num[i] = length(which(pro_matrix_sel[i,] != 0))
  }
  return(list(sel_mrna_name, sel_mrna_score, tar_num))
}

medfc_test = function(vset, mirna_name, sel_mrna_name1, tar_num, num){
  tar_num1 = tar_num[which(tar_num != 0)]
  thres = pmin(tar_num1, num)
  mirna_vali = as.vector(vset[1,-1])
  vset_name = as.vector(vset[-1,1])
  vset_score = vset[-1,-1]
  mirna_pos = match(mirna_vali, mirna_name)
  mirna_pos_use = mirna_pos[which(is.na(mirna_pos) != TRUE)]
  vset_score_use2 = vset_score[,which(is.na(mirna_pos) != TRUE)]
  vset_score_use = vset_score_use2[,which(tar_num != 0)]
  sel_mrna_name = sel_mrna_name1[,which(tar_num != 0)]
  cc = rep(0, ncol(vset_score_use))
  for(i in 1:ncol(vset_score_use)){
    temp = match(as.numeric(sel_mrna_name[1:thres[i],i]), as.numeric(vset_name))
    cc[i] = length(which(is.na(temp) != TRUE))
  }
  pos = which(cc != 0)
  vset_score_use = vset_score_use[, pos]
  sel_mrna_name = sel_mrna_name[, pos]
  thres = thres[pos]
  mrna_count = rep(0, ncol(vset_score_use))
  mrnafc_matrix = matrix(data = 0, ncol = ncol(vset_score_use), nrow = nrow(sel_mrna_name), byrow = TRUE)
  for(i in 1:ncol(vset_score_use)){
    temp = match(as.numeric(sel_mrna_name[1:thres[i],i]), as.numeric(vset_name))
    #mrna_count[i] = length(which(is.na(temp) != TRUE))
    temp2 = vset_score_use[temp[which(is.na(temp) != TRUE)], i]
    mrna_count[i] = length(which(is.na(temp2) != TRUE))
    mrnafc_matrix[1:mrna_count[i], i] = temp2[which(is.na(temp2) != TRUE)]
  }
  return(list(mrna_count, mrnafc_matrix))
}

medfc_vector = function(mrna_count, mrnafc_matrix){
  median_mat = matrix(data = 0, ncol = ncol(mrnafc_matrix), nrow = nrow(mrnafc_matrix), byrow = TRUE)
  for(i in 1:ncol(mrnafc_matrix)){
    for(j in 1:nrow(mrnafc_matrix)){
      median_mat[j,i] = median(as.numeric(mrnafc_matrix[1:j,i]))
    }
  }
  fc_vector = rep(0, max(mrna_count))
  for(i in 1:length(fc_vector)){
    temp1 = which(mrna_count >= i)
    fc_vector[i] = mean(as.numeric(median_mat[i,temp1]))
  }
  return(fc_vector)
}

result_summary = function(name_cancer, mirna_cancer, mrna_cancer, score_matrix, thres_num, Vset, method){
  score_list = matrix_transfer(score_matrix)
  x = name_func(name_cancer, mirna_cancer, mrna_cancer)
  z1 = mirna_matrix(name_cancer, mirna_cancer, x[[1]], score_list[[1]])
  z2 = mrna_matrix(name_cancer, mrna_cancer, x[[2]], score_list[[2]])
  if(method == "miracle"){
    z3 = mirna_mrna_data(z1[[1]], z2[[1]], z1[[2]], z2[[2]], z2[[3]], z2[[4]], 1)
    z7 = mirna_mrna_loc(z3[[3]], z3[[4]], score_list[[1]], score_list[[2]], score_list[[3]], z3[[6]])
    z8 = sum_miracle(z3[[1]], z3[[2]], z7[[3]])
    z9 = miracle_medfc(z3[[1]], z3[[2]], z3[[5]], z7[[3]], z8, z3[[3]], Vset, thres_num) #get the result of miracle
    return(list(z9, z3))
  }
  else if(method == "promise"){
    z3 = mirna_mrna_data(z1[[1]], z2[[1]], z1[[2]], z2[[2]], z2[[3]], z2[[4]], 1)
    z6 = mirna_mrna_loc(z3[[3]], z3[[4]], score_list[[1]], score_list[[2]], score_list[[3]], z3[[6]])
    z10 = promise_medfc(z3[[1]], z3[[2]], z3[[5]], z6[[3]], z3[[3]], Vset, thres_num) #result of promise all
    return(list(z10, z3))
  }
  else if(method == "sequence"){
    z4 = mirna_mrna_data(z1[[1]], z2[[1]], z1[[2]], z2[[2]], z2[[3]], z2[[4]], 0.8)
    z14 = mirna_mrna_loc(z4[[3]], z4[[4]], score_list[[1]], score_list[[2]], score_list[[3]], z4[[6]])
    z12 = targetscan_medfc(z14[[3]], z4[[5]], z4[[3]], Vset, thres_num) #get the rankers of targetscan_qMRE
    return(list(z12, z4))
  }
  else if(method == "pearson"){
    z3 = mirna_mrna_data(z1[[1]], z2[[1]], z1[[2]], z2[[2]], z2[[3]], z2[[4]], 1)
    z7 = mirna_mrna_loc(z3[[3]], z3[[4]], score_list[[1]], score_list[[2]], score_list[[3]], z3[[6]])
    z11 = pearson_medfc(z3[[1]], z3[[2]], z3[[5]], z3[[3]], z7[[3]], Vset, thres_num) #get the result of Pearson correlation
    return(list(z11, z3))
  }
  else if(method == "lasso"){
    z3 = mirna_mrna_data(z1[[1]], z2[[1]], z1[[2]], z2[[2]], z2[[3]], z2[[4]], 1)
    z7 = mirna_mrna_loc(z3[[3]], z3[[4]], score_list[[1]], score_list[[2]], score_list[[3]], z3[[6]])
    z12 = LASSO_medfc(z3[[1]], z3[[2]], z3[[5]], z3[[3]], z7[[3]], Vset, thres_num) #get the result of LASSO
    return(list(z12, z3))
  }
  else if(method == "elastic"){
    z3 = mirna_mrna_data(z1[[1]], z2[[1]], z1[[2]], z2[[2]], z2[[3]], z2[[4]], 1)
    z7 = mirna_mrna_loc(z3[[3]], z3[[4]], score_list[[1]], score_list[[2]], score_list[[3]], z3[[6]])
    z13 = ELASTIC_medfc(z3[[1]], z3[[2]], z3[[5]], z3[[3]], z7[[3]], Vset, thres_num) #get the result of elastic
    return(list(z13, z3))
  }
  else if(method == "Zscore"){
    z3 = mirna_mrna_data(z1[[1]], z2[[1]], z1[[2]], z2[[2]], z2[[3]], z2[[4]], 1)
    z7 = mirna_mrna_loc(z3[[3]], z3[[4]], score_list[[1]], score_list[[2]], score_list[[3]], z3[[6]])
    z13 = Zscore_medfc(z3[[1]], z3[[2]], z3[[5]], z3[[3]], z7[[3]], Vset, thres_num) #get the result of Zscore
    return(list(z13, z3))
  }
}

medfc_summary = function(vset){
  num = 800
  z9 = result_summary(name_cancer, mirna_cancer, mrna_cancer, CWCS_all, num, vset, "miracle")
  z9_con = result_summary(name_cancer, mirna_cancer, mrna_cancer, CWCS_conserved, num, vset, "miracle")
  z10 = result_summary(name_cancer, mirna_cancer, mrna_cancer, qMRE_conserved, num, vset, "promise") #promise conserved
  z11 = result_summary(name_cancer, mirna_cancer, mrna_cancer, combine, num, vset, "pearson") #Pearson
  z12 = result_summary(name_cancer, mirna_cancer, mrna_cancer, combine, num, vset, "lasso") #LASSO
  z13 = result_summary(name_cancer, mirna_cancer, mrna_cancer, combine, num, vset, "elastic") #Elastic net
  z14 = result_summary(name_cancer, mirna_cancer, mrna_cancer, combine, num, vset, "Zscore") #Zscore
  z15 = result_summary(name_cancer, mirna_cancer, mrna_cancer, CWCS_conserved, num, vset, "sequence") #CWCS conserved
  z16 = result_summary(name_cancer, mirna_cancer, mrna_cancer, microT, num, vset, "sequence") #microT
  
  z21 = medfc_test(vset, z9[[2]][[3]], z9[[1]][[1]], z9[[1]][[3]], num)
  z22 = medfc_vector(z21[[1]], z21[[2]]) #for miracle
  z21_1 = medfc_test(vset, z9_con[[2]][[3]], z9_con[[1]][[1]], z9_con[[1]][[3]], num)
  z22_1 = medfc_vector(z21_1[[1]], z21_1[[2]]) #for miracle conserved
  z23 = medfc_test(vset, z10[[2]][[3]], z10[[1]][[1]], z10[[1]][[5]], num)
  z24 = medfc_vector(z23[[1]], z23[[2]]) #for promise_mrna conserved
  z25 = medfc_test(vset, z11[[2]][[3]], z11[[1]][[1]], z11[[1]][[3]], num)
  z26 = medfc_vector(z25[[1]], z25[[2]]) #for pearson
  z27 = medfc_test(vset, z12[[2]][[3]], z12[[1]][[1]], z12[[1]][[3]], num)
  z28 = medfc_vector(z27[[1]], z27[[2]]) #for LASSO
  z29 = medfc_test(vset, z13[[2]][[3]], z13[[1]][[1]], z13[[1]][[3]], num)
  z30 = medfc_vector(z29[[1]], z29[[2]]) #for elastic net
  z31 = medfc_test(vset, z14[[2]][[3]], z14[[1]][[1]], z14[[1]][[3]], num)
  z32 = medfc_vector(z31[[1]], z31[[2]]) #for zscore
  z33 = medfc_test(vset, z15[[2]][[3]], z15[[1]][[1]], z15[[1]][[3]], num)
  z34 = medfc_vector(z33[[1]], z33[[2]]) #for CWCS conserved
  z35 = medfc_test(vset, z16[[2]][[3]], z16[[1]][[1]], z16[[1]][[3]], num)
  z36 = medfc_vector(z35[[1]], z35[[2]]) #for microT
  
  final_median = rbind(z22[1:500], z22_1[1:500], z24[1:500], z26[1:500], z28[1:500], z30[1:500], z32[1:500],
                       z34[1:500], z36[1:500])
  colnames(final_median) = seq(1, 500, 1)
  rownames(final_median) = c('miRACLe', 'miRACLe_cons', 'ProMISe_mrna_cons', 'PearsonCorrelation', 'LASSO', 
                             'Elastic-net', 'Zscore', 'TargetScan7_CWCS_cons', 'DIANA-microT-CDS')
  return(final_median)
}

wil_test = function(valid){
  result = matrix(data = 0, ncol = 2, nrow = 7)
  for(i in 3:9){
    result[(i - 2), 1] = wilcox.test(as.vector(valid[1, ]), as.vector(valid[i, ]), alternative = "less", paired = TRUE, correct = FALSE)$p.value
    result[(i - 2), 2] = wilcox.test(as.vector(valid[2, ]), as.vector(valid[i, ]), alternative = "less", paired = TRUE, correct = FALSE)$p.value
  }
  
  colnames(result) = c('miRACLe', 'miRACLe_cons')
  rownames(result) = c('ProMISe_mrna_cons', 'TargetScan7_CWCS_cons', 'PearsonCorrelation', 
                       'LASSO', 'Elastic-net', 'Zscore', 'DIANA-microT-CDS')
  return(result)
}

###################################################################################################
####################################### DATA INPUTS ###############################################
###################################################################################################

CWCS_all = as.matrix(read.table("TargetScan7_CWCS.txt", head=TRUE, sep = "\t"))
CWCS_conserved = as.matrix(read.table("TargetScan7_CWCS_cons.txt", head=TRUE, sep = "\t"))
qMRE_conserved = as.matrix(read.table("TargetScan7_qMRE_cons.txt", head = TRUE, sep = "\t"))
microT = as.matrix(read.table("DIANA_microT_CDS.txt", head = TRUE, sep = "\t"))
combine = as.matrix(read.table("Combine_MMIs.txt", head = TRUE, sep = "\t"))

vset = as.matrix(read.table("Transet_multi.txt", head = FALSE, sep = "\t"))

name_cancer = as.matrix(read.table("TCGA_HNSC_sampleMatch.txt", head = TRUE, sep = "\t"))
mirna_cancer = as.matrix(read.table("HNSC_miRNA_expression.txt", head = FALSE, sep = "\t"))
mrna_cancer = as.matrix(read.table("HNSC_mRNA_expression.txt", head = FALSE, sep = "\t"))

###################################################################################################
######################################### MAIN code ###############################################
###################################################################################################

T1_result = medfc_summary(vset)
wil_T1 = wil_test(T1_result)
write.table(T1_result, file = "median fold change.txt", quote = FALSE, sep = "\t")
write.table(wil_T1, file = "p-value of Wilcoxon test.txt", quote = FALSE, sep = "\t")

