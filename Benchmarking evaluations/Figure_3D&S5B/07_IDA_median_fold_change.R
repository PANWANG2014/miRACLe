###################################################################################################
## This Rscript generates benchmarking evaluation of MCC/NCI60 - median fold change

###################################################################################################
## Loading required packages
library(pcalg)

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

RNA_name_sub = function(mirna_name){
  mirna = c()
  for(i in 1:length(mirna_name)){
    temp = unlist(strsplit(mirna_name[i],"\\|"))[1]
    mirna = c(mirna, temp)
  }
  return(mirna)
}

miracle_medfc = function(mirna_use1, mrna_use1, mrna_name, mirna_name, CWCS, vset, num){
  mirna_use = matrix(as.numeric(mirna_use1), nrow = nrow(mirna_use1))
  mrna_use = matrix(as.numeric(mrna_use1), nrow = nrow(mrna_use1))
  pro_matrix_agg = matrix(data = 0, ncol = nrow(mrna_use), nrow = nrow(mirna_use))
  for(m in 1:ncol(mirna_use)){
    pro_matrix_temp1 = (mirna_use[, m] %*% t(mrna_use[, m])) * CWCS
    sumup = sum(pro_matrix_temp1)
    pro_matrix_temp = pro_matrix_temp1 / sumup
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

miracleDE_medfc = function(mirna_use1, mrna_use1, mrna_name, mirna_name, mrna_fullname, CWCS, diff_mirna, diff_mrna, vset, num){
  mirna_use = matrix(as.numeric(mirna_use1), nrow = nrow(mirna_use1))
  mrna_use = matrix(as.numeric(mrna_use1), nrow = nrow(mrna_use1))
  diff_mirna_sim = RNA_name_sub(diff_mirna)
  pos1 = match(diff_mirna_sim, mirna_name)
  pos2 = match(diff_mrna, mrna_fullname)
  
  pos3 = pos1[which(is.na(pos1) != TRUE)]
  pos4 = pos2[which(is.na(pos2) != TRUE)]
  diff_mirna_use = diff_mirna_sim[which(is.na(pos1) != TRUE)]
  diff_mrna_use = diff_mrna[which(is.na(pos2) != TRUE)]
  
  mirna_use_new = mirna_use[pos3, ]
  mrna_use_new = mrna_use[pos4, ]
  mrna_name_new = mrna_name[pos4, ]
  mirna_name_new = mirna_name[pos3]
  mrna_fullname_new = mrna_fullname[pos4]
  CWCS_use = CWCS[pos3, pos4]
  
  pro_matrix_agg = matrix(data = 0, ncol = nrow(mrna_use_new), nrow = nrow(mirna_use_new))
  for(m in 1:ncol(mirna_use_new)){
    pro_matrix_temp1 = (mirna_use_new[, m] %*% t(mrna_use_new[, m])) * CWCS_use
    sumup = sum(pro_matrix_temp1)
    pro_matrix_temp = pro_matrix_temp1 / sumup
    pro_matrix_agg = pro_matrix_agg + pro_matrix_temp
  }
  #thres_num = ceiling(length(which(pro_matrix_agg != 0))/100)
  vset_name = as.vector(vset[1,-1])
  temp = match(vset_name, mirna_name_new)
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
  mirna_name_col = mirna_name_new[which(is.na(match(mirna_name_new, vset_name)) != TRUE)]
  colnames(sel_mrna_name) = mirna_name_col
  colnames(sel_mrna_score) = mirna_name_col
  return(list(sel_mrna_name, sel_mrna_score, tar_num, pro_matrix_agg))
}

cIDA_medfc = function(mirna_use1, mrna_use1, mrna_name, mirna_name, mrna_fullname, diff_mirna, diff_mrna, alpha, vset, num){
  mirna_use = matrix(as.numeric(mirna_use1), nrow = nrow(mirna_use1))
  mrna_use = matrix(as.numeric(mrna_use1), nrow = nrow(mrna_use1))
  diff_mirna_sim = RNA_name_sub(diff_mirna)
  pos1 = match(diff_mirna_sim, mirna_name)
  pos2 = match(diff_mrna, mrna_fullname)
  
  pos3 = pos1[which(is.na(pos1) != TRUE)]
  pos4 = pos2[which(is.na(pos2) != TRUE)]
  diff_mirna_use = diff_mirna_sim[which(is.na(pos1) != TRUE)]
  diff_mrna_use = diff_mrna[which(is.na(pos2) != TRUE)]
  
  mirna_use_new = mirna_use[pos3, ]
  mrna_use_new = mrna_use[pos4, ]
  mrna_name_new = mrna_name[pos4, ]
  mirna_name_new = mirna_name[pos3]
  mrna_fullname_new = mrna_fullname[pos4]
  
  data = cbind(t(mirna_use_new), t(mrna_use_new))
  cause = 1:nrow(mirna_use_new)
  effect = (1 + nrow(mirna_use_new)): (nrow(mrna_use_new) + nrow(mirna_use_new))
  multiset = character(0)
  result = matrix(nrow = nrow(mrna_use_new), ncol = nrow(mirna_use_new))
  suffStat = list(C = cor(data), n = nrow(data))
  pcFit = pc(suffStat, indepTest = gaussCItest, p = ncol(data), alpha = alpha)
  for(l in cause){
    caef = idaFast(l, effect, cov(data), pcFit@graph)
    caef1 = matrix(nrow = length(effect), ncol = 1)
    for (k in 1:length(effect)){
      caefabs = abs(caef)
      index = which(caefabs == min(caefabs[k,]), arr.ind = TRUE)
      pos = index[1, 2]
      caef1[k, ] = caef[k, pos]
    }
    result[, l] = caef1
    print(l)
  }
  result = t(result)
  
  vset_name = as.vector(vset[1,-1])
  temp = match(vset_name, mirna_name_new)
  mirna_pos = temp[which(is.na(temp) != TRUE)]
  pro_matrix_sel = result[mirna_pos, ]
  sel_mrna_name = matrix(data = 0, ncol = length(mirna_pos), nrow = num, byrow = TRUE)
  sel_mrna_score = matrix(data = 0, ncol = length(mirna_pos), nrow = num, byrow = TRUE)
  tar_num = rep(0, length(mirna_pos))
  for(i in 1: length(mirna_pos)){
    sel_mrna_score[, i] = sort(pro_matrix_sel[i, ], decreasing = TRUE)[1:num]
    temp = order(pro_matrix_sel[i, ], decreasing = TRUE)[1:num]
    sel_mrna_name[, i] = as.numeric(mrna_name[temp, 2])
    tar_num[i] = length(which(pro_matrix_sel[i,] != 0))
  }
  mirna_name_col = mirna_name_new[which(is.na(match(mirna_name_new, vset_name)) != TRUE)]
  colnames(sel_mrna_name) = mirna_name_col
  colnames(sel_mrna_score) = mirna_name_col
  return(list(sel_mrna_name, sel_mrna_score, tar_num, result))
}

mIDA_medfc = function(mirna_use1, mrna_use1, mrna_name, mirna_name, mrna_fullname, score_matrix, diff_mirna, diff_mrna, alpha, vset, num){
  mirna_use = matrix(as.numeric(mirna_use1), nrow = nrow(mirna_use1))
  mrna_use = matrix(as.numeric(mrna_use1), nrow = nrow(mrna_use1))
  diff_mirna_sim = RNA_name_sub(diff_mirna)
  pos1 = match(diff_mirna_sim, mirna_name)
  pos2 = match(diff_mrna, mrna_fullname)
  
  pos3 = pos1[which(is.na(pos1) != TRUE)]
  pos4 = pos2[which(is.na(pos2) != TRUE)]
  diff_mirna_use = diff_mirna_sim[which(is.na(pos1) != TRUE)]
  diff_mrna_use = diff_mrna[which(is.na(pos2) != TRUE)]
  
  mirna_use_new = mirna_use[pos3, ]
  mrna_use_new = mrna_use[pos4, ]
  mrna_name_new = mrna_name[pos4, ]
  mirna_name_new = mirna_name[pos3]
  mrna_fullname_new = mrna_fullname[pos4]
  score_matrix_use = score_matrix[pos3, pos4]
  
  data = cbind(t(mirna_use_new), t(mrna_use_new))
  cause = 1:nrow(mirna_use_new)
  effect = (1 + nrow(mirna_use_new)): (nrow(mrna_use_new) + nrow(mirna_use_new))
  multiset = character(0)
  result = matrix(nrow = nrow(mrna_use_new), ncol = nrow(mirna_use_new))
  suffStat = list(C = cor(data), n = nrow(data))
  pcFit = pc(suffStat, indepTest = gaussCItest, p = ncol(data), alpha = alpha)
  for(l in cause){
    caef = idaFast(l, effect, cov(data), pcFit@graph)
    caef1 = matrix(nrow = length(effect), ncol = 1)
    for (k in 1:length(effect)){
      caefabs = abs(caef)
      index = which(caefabs == min(caefabs[k,]), arr.ind = TRUE)
      pos = index[1, 2]
      caef1[k, ] = caef[k, pos]
    }
    result[, l] = caef1
    print(l)
  }
  result = t(result)
  res = result * score_matrix_use
  
  vset_name = as.vector(vset[1,-1])
  temp = match(vset_name, mirna_name_new)
  mirna_pos = temp[which(is.na(temp) != TRUE)]
  pro_matrix_sel = res[mirna_pos, ]
  sel_mrna_name = matrix(data = 0, ncol = length(mirna_pos), nrow = num, byrow = TRUE)
  sel_mrna_score = matrix(data = 0, ncol = length(mirna_pos), nrow = num, byrow = TRUE)
  tar_num = rep(0, length(mirna_pos))
  for(i in 1: length(mirna_pos)){
    sel_mrna_score[, i] = sort(pro_matrix_sel[i, ], decreasing = TRUE)[1:num]
    temp = order(pro_matrix_sel[i, ], decreasing = TRUE)[1:num]
    sel_mrna_name[, i] = as.numeric(mrna_name[temp, 2])
    tar_num[i] = length(which(pro_matrix_sel[i,] != 0))
  }
  mirna_name_col = mirna_name_new[which(is.na(match(mirna_name_new, vset_name)) != TRUE)]
  colnames(sel_mrna_name) = mirna_name_col
  colnames(sel_mrna_score) = mirna_name_col
  return(list(sel_mrna_name, sel_mrna_score, tar_num, res))
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

result_summary = function(name_cancer, mirna_cancer, mrna_cancer, matrix, diff_mirna, diff_mrna, thres_num, alpha, Vset, method){
  score_list = matrix_transfer(matrix)
  x = name_func(name_cancer, mirna_cancer, mrna_cancer)
  z1 = mirna_matrix(name_cancer, mirna_cancer, x[[1]], score_list[[1]])
  z2 = mrna_matrix(name_cancer, mrna_cancer, x[[2]], score_list[[2]])
  if(method == "miracle"){
    z3 = mirna_mrna_data(z1[[1]], z2[[1]], z1[[2]], z2[[2]], z2[[3]], z2[[4]], 1)
    z7 = mirna_mrna_loc(z3[[3]], z3[[4]], score_list[[1]], score_list[[2]], score_list[[3]], z3[[6]])
    z9 = miracle_medfc(z3[[1]], z3[[2]], z3[[5]], z3[[3]], z7[[3]], Vset, thres_num) #miracle
    return(list(z9, z3))
  }
  if(method == "miracleDE"){
    z3 = mirna_mrna_data(z1[[1]], z2[[1]], z1[[2]], z2[[2]], z2[[3]], z2[[4]], 1)
    z7 = mirna_mrna_loc(z3[[3]], z3[[4]], score_list[[1]], score_list[[2]], score_list[[3]], z3[[6]])
    z9 = miracleDE_medfc(z3[[1]], z3[[2]], z3[[5]], z3[[3]], z3[[6]], z7[[3]], diff_mirna, diff_mrna, Vset, thres_num) #miracleDE
    return(list(z9, z3))
  }
  else if(method == "cIDA"){
    z3 = mirna_mrna_data(z1[[1]], z2[[1]], z1[[2]], z2[[2]], z2[[3]], z2[[4]], 0.8)
    z7 = mirna_mrna_loc(z3[[3]], z3[[4]], score_list[[1]], score_list[[2]], score_list[[3]], z3[[6]])
    z12 = cIDA_medfc(z3[[1]], z3[[2]], z3[[5]], z3[[3]], z3[[6]], diff_mirna, diff_mrna, alpha, Vset, thres_num) #cIDA
    return(list(z12, z3))
  }
  else if(method == "mIDA"){
    z3 = mirna_mrna_data(z1[[1]], z2[[1]], z1[[2]], z2[[2]], z2[[3]], z2[[4]], 0.8)
    z7 = mirna_mrna_loc(z3[[3]], z3[[4]], score_list[[1]], score_list[[2]], score_list[[3]], z3[[6]])
    z12 = mIDA_medfc(z3[[1]], z3[[2]], z3[[5]], z3[[3]], z3[[6]], z7[[3]], diff_mirna, diff_mrna, alpha, Vset, thres_num) #mIDA
    return(list(z12, z3))
  }
}

medfc_summary = function(alpha, vset){
  num = 800
  z10 = result_summary(name_cancer, mirna_cancer, mrna_cancer, CWCS_all, diff_mirna, diff_mrna, num, alpha, vset, "miracle")
  z11 = result_summary(name_cancer, mirna_cancer, mrna_cancer, CWCS_conserved, diff_mirna, diff_mrna, num, alpha, vset, "miracle")
  z12 = result_summary(name_cancer, mirna_cancer, mrna_cancer, CWCS_all, diff_mirna, diff_mrna, num, alpha, vset, "miracleDE")
  z13 = result_summary(name_cancer, mirna_cancer, mrna_cancer, CWCS_conserved, diff_mirna, diff_mrna, num, alpha, vset, "miracleDE")
  z14 = result_summary(name_cancer, mirna_cancer, mrna_cancer, combine, diff_mirna, diff_mrna, num, alpha, vset, "cIDA")
  z15 = result_summary(name_cancer, mirna_cancer, mrna_cancer, combine, diff_mirna, diff_mrna, num, alpha, vset, "mIDA")
  
  z21 = medfc_test(vset, z10[[2]][[3]], z10[[1]][[1]], z10[[1]][[3]], num)
  z22 = medfc_vector(z21[[1]], z21[[2]]) #for miracle
  z23 = medfc_test(vset, z11[[2]][[3]], z11[[1]][[1]], z11[[1]][[3]], num)
  z24 = medfc_vector(z23[[1]], z23[[2]]) #for miracle conserved
  z25 = medfc_test(vset, z12[[2]][[3]], z12[[1]][[1]], z12[[1]][[3]], num)
  z26 = medfc_vector(z25[[1]], z25[[2]]) #for miracleDE
  z27 = medfc_test(vset, z13[[2]][[3]], z13[[1]][[1]], z13[[1]][[3]], num)
  z28 = medfc_vector(z27[[1]], z27[[2]]) #for miracleDE conserved
  z29 = medfc_test(vset, z14[[2]][[3]], z14[[1]][[1]], z14[[1]][[3]], num)
  z30 = medfc_vector(z29[[1]], z29[[2]]) #for cIDA
  z31 = medfc_test(vset, z15[[2]][[3]], z15[[1]][[1]], z15[[1]][[3]], num)
  z32 = medfc_vector(z31[[1]], z31[[2]]) #for mIDA
  
  final_median = rbind(z22[1:500], z24[1:500], z26[1:500], z28[1:500], z30[1:500], z32[1:500])
  colnames(final_median) = seq(1, 500, 1)
  rownames(final_median) = c('miracle_all', 'miracle_conserved', 'miracleDE_all', 'miracleDE_conserved', 
                             'cIDA', 'mIDA')
  return(final_median)
}

wil_test = function(valid){
  num = length(which(is.na(valid[4, ]) != TRUE))
  a1 = wilcox.test(as.vector(valid[1, ]), as.vector(valid[3, ]), alternative = "less", paired = TRUE, correct = FALSE)$p.value
  a3 = wilcox.test(as.vector(valid[3, ]), as.vector(valid[5, ]), alternative = "less", paired = TRUE, correct = FALSE)$p.value
  a4 = wilcox.test(as.vector(valid[3, ]), as.vector(valid[6, ]), alternative = "less", paired = TRUE, correct = FALSE)$p.value
  a2 = wilcox.test(as.vector(valid[2, 1:num]), as.vector(valid[4, 1:num]), alternative = "less", paired = TRUE, correct = FALSE)$p.value
  a5 = wilcox.test(as.vector(valid[4, 1:num]), as.vector(valid[5, 1:num]), alternative = "less", paired = TRUE, correct = FALSE)$p.value
  a6 = wilcox.test(as.vector(valid[4, 1:num]), as.vector(valid[6, 1:num]), alternative = "less", paired = TRUE, correct = FALSE)$p.value
  
  a = c(a1, a2, a3, a4, a5, a6)
  result = matrix(a, ncol = 6)
  colnames(result) = c('miRACLe vs miRACLe_DE', 'miRACLe_cons vs miRACLe_cons_DE', 'miRACLe_DE vs cIDA', 'miRACLe_DE vs mIDA',
                       'miRACLe_cons_DE vs cIDA', 'miRACLe_cons_DE vs mIDA')
  return(result)
}

###################################################################################################
####################################### DATA INPUTS ###############################################
###################################################################################################

CWCS_all = as.matrix(read.table("TargetScan7_CWCS.txt", head = TRUE, sep = "\t"))
CWCS_conserved = as.matrix(read.table("TargetScan7_CWCS_cons.txt", head = TRUE, sep = "\t"))
combine = as.matrix(read.table("Combine_MMIs.txt", head = TRUE, sep = "\t"))
symbol_to_ID = as.matrix(read.table("Symbol_to_ID.txt", head = TRUE, sep = "\t"))

name_cancer = as.matrix(read.table("MCC_tumor_sampleMatch.txt", head = TRUE, sep = "\t"))
mirna_cancer = as.matrix(read.table("MCC_miRNA_expression.txt", head = FALSE, sep = "\t"))
mrna_cancer = as.matrix(read.table("MCC_mRNA_expression.txt", head = FALSE, sep = "\t"))
diff_mirna = as.matrix(read.table("MCC_diffmiRNAlist.txt", head = FALSE, sep = "\t"))
diff_mrna = as.matrix(read.table("MCC_diffmRNAlist.txt", head = FALSE, sep = "\t"))

vset = as.matrix(read.table("Transet_multi.txt", head = FALSE, sep = "\t"))

# name_cancer = as.matrix(read.table("NCI60_mesenchymal_sampleMatch.txt", head = TRUE, sep = "\t"))
# mirna_cancer = as.matrix(read.table("NCI60_miRNA_expression.txt", head = FALSE, sep = "\t"))
# mrna_cancer = as.matrix(read.table("NCI60_mRNA_expression.txt", head = FALSE, sep = "\t"))
# diff_mirna = as.matrix(read.table("NCI60_diffmiRNAlist.txt", head = FALSE, sep = "\t"))
# diff_mrna = as.matrix(read.table("NCI60_diffmRNAlist.txt", head = FALSE, sep = "\t"))

###################################################################################################
######################################### MAIN code ###############################################
###################################################################################################

alpha = 0.01

T1_result = medfc_summary(alpha, vset)
wil_T1 = wil_test(T1_result)
write.table(T1_result, file = "median fold change.txt", quote = FALSE, sep = "\t")
write.table(wil_T1, file = "p-value of Wilcoxon test.txt", quote = FALSE, sep = "\t", row.names = FALSE)

