###################################################################################################
## This Rscript generates benchmarking evaluation of TCGA - precision analysis

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

miracle_pop = function(mirna_use1, mrna_use1, CWCS, sumup, thres_num){
  mirna_use = matrix(as.numeric(mirna_use1), nrow = nrow(mirna_use1))
  mrna_use = matrix(as.numeric(mrna_use1), nrow = nrow(mrna_use1))
  pro_matrix_temp = matrix(data = 0, ncol = nrow(mrna_use), nrow = nrow(mirna_use), byrow = TRUE)
  pro_matrix_agg = matrix(data = 0, ncol = nrow(mrna_use), nrow = nrow(mirna_use), byrow = TRUE)
  for(m in 1:ncol(mirna_use)){
    pro_matrix_temp = (mirna_use[, m] %*% t(mrna_use[, m])) * CWCS / sumup[m]
    pro_matrix_agg = pro_matrix_agg + pro_matrix_temp
  }
  order_agg = order(pro_matrix_agg, decreasing = TRUE)[1 : thres_num]
  return(order_agg)
}

promise_pop = function(mirna_use1, mrna_use1, CWCS, thres_num){
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
  order_mrna = order(pro_matrix_mrna, decreasing = TRUE)[1 : thres_num]
  order_full = order(pro_matrix_full, decreasing = TRUE)[1 : thres_num]
  return(list(order_mrna, order_full))
}

pearson_cal = function(mirna_use1, mrna_use1, thres_num, score_matrix){
  mirna_use = matrix(as.numeric(mirna_use1),nrow = nrow(mirna_use1))
  mrna_use = matrix(as.numeric(mrna_use1),nrow = nrow(mrna_use1))
  corMat = cor(t(mirna_use), t(mrna_use))
  corMat = corMat * score_matrix
  sort_one_rank = order(corMat, decreasing=FALSE)[1:thres_num]
  sort_one_prob = sort(corMat, decreasing=FALSE)[1:thres_num]
  return(list(sort_one_rank, sort_one_prob, corMat))
}

LASSO_cal = function(mirna_use1, mrna_use1, thres_num, score_matrix){
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
  sort_one_rank = order(res, decreasing = FALSE)[1 : thres_num]
  sort_one_prob = sort(res, decreasing = FALSE)[1 : thres_num]
  return(list(sort_one_rank, sort_one_prob, res))
}

ELASTIC_cal = function(mirna_use1, mrna_use1, thres_num, score_matrix){ 
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
  sort_one_rank = order(res, decreasing = FALSE)[1 : thres_num]
  sort_one_prob = sort(res, decreasing = FALSE)[1 : thres_num]
  return(list(sort_one_rank, sort_one_prob, res))
}

Zscore_cal = function(mirna_use1, mrna_use1, thres_num, score_matrix){
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
  sort_one_rank = order(res, decreasing = TRUE)[1 : thres_num]
  sort_one_prob = sort(res, decreasing = TRUE)[1 : thres_num]
  return(list(sort_one_rank, sort_one_prob, res))
}

targetscan_rank = function(CWCS, thres_num){
  target_rank = order(CWCS, decreasing = TRUE)[1 : thres_num]
  target_prob = sort(CWCS, decreasing = TRUE)[1 : thres_num]
  return(list(target_rank, target_prob))
}

compare_result = function(rank_matrix, Vset, num){
  Vset = as.numeric(Vset)
  rank = as.numeric(rank_matrix)[1:num]
  rank_result = rank %in% Vset
  num1 = num/100
  rank_str = rep(0, num1)
  rank_str_percent = rep(0, num1)
  for(i in 1:length(rank_str)){
    temp1 = 100 * i
    rank_str[i] = sum(rank_result[1:temp1])
    rank_str_percent[i] = rank_str[i]/temp1
  }
  return(list(rank_str, rank_str_percent))
}

result_summary = function(name_cancer, mirna_cancer, mrna_cancer, score_matrix, thres_num, Vset, method){
  score_list = matrix_transfer(score_matrix)
  x = name_func(name_cancer, mirna_cancer, mrna_cancer)
  z1 = mirna_matrix(name_cancer, mirna_cancer, x[[1]], score_list[[1]])
  z2 = mrna_matrix(name_cancer, mrna_cancer, x[[2]], score_list[[2]])
  if(method == "miracle"){
    z3 = mirna_mrna_data(z1[[1]], z2[[1]], z1[[2]], z2[[2]], z2[[3]], z2[[4]], 1)
    z5 = MMI_location(Vset, z3[[3]], z3[[5]])
    z7 = mirna_mrna_loc(z3[[3]], z3[[4]], score_list[[1]], score_list[[2]], score_list[[3]], z3[[6]])
    z8 = sum_miracle(z3[[1]], z3[[2]], z7[[3]])
    z9 = miracle_pop(z3[[1]], z3[[2]], z7[[3]], z8,  thres_num) #result of miracle
    return(list(z9, z5))
  }
  else if(method == "promise"){
    z3 = mirna_mrna_data(z1[[1]], z2[[1]], z1[[2]], z2[[2]], z2[[3]], z2[[4]], 1)
    z5 = MMI_location(Vset, z3[[3]], z3[[5]]) 
    z6 = mirna_mrna_loc(z3[[3]], z3[[4]], score_list[[1]], score_list[[2]], score_list[[3]], z3[[6]])
    z10 = promise_pop(z3[[1]], z3[[2]], z6[[3]], thres_num) #result of promise
    return(list(z10, z5))
  }
  else if(method == "sequence"){
    z4 = mirna_mrna_data(z1[[1]], z2[[1]], z1[[2]], z2[[2]], z2[[3]], z2[[4]], 0.8)
    z5_tar = MMI_location(Vset, z4[[3]], z4[[5]]) #this is used to locate the validation set
    z14 = mirna_mrna_loc(z4[[3]], z4[[4]], score_list[[1]], score_list[[2]], score_list[[3]], z4[[6]])
    z12 = targetscan_rank(z14[[3]], thres_num) #result of sequence method
    return(list(z12, z5_tar))
  }
  else if(method == "pearson"){
    z3 = mirna_mrna_data(z1[[1]], z2[[1]], z1[[2]], z2[[2]], z2[[3]], z2[[4]], 1)
    z5 = MMI_location(Vset, z3[[3]], z3[[5]]) 
    z7 = mirna_mrna_loc(z3[[3]], z3[[4]], score_list[[1]], score_list[[2]], score_list[[3]], z3[[6]])
    z11 = pearson_cal(z3[[1]], z3[[2]], thres_num, z7[[3]]) #result of Pearson correlation
    return(list(z11, z5))
  }
  else if(method == "lasso"){
    z3 = mirna_mrna_data(z1[[1]], z2[[1]], z1[[2]], z2[[2]], z2[[3]], z2[[4]], 1)
    z5 = MMI_location(Vset, z3[[3]], z3[[5]])
    z7 = mirna_mrna_loc(z3[[3]], z3[[4]], score_list[[1]], score_list[[2]], score_list[[3]], z3[[6]])
    z12 = LASSO_cal(z3[[1]], z3[[2]], thres_num, z7[[3]]) #result of LASSO
    return(list(z12, z5))
  }
  else if(method == "elastic"){
    z3 = mirna_mrna_data(z1[[1]], z2[[1]], z1[[2]], z2[[2]], z2[[3]], z2[[4]], 1)
    z5 = MMI_location(Vset, z3[[3]], z3[[5]])
    z7 = mirna_mrna_loc(z3[[3]], z3[[4]], score_list[[1]], score_list[[2]], score_list[[3]], z3[[6]])
    z13 = ELASTIC_cal(z3[[1]], z3[[2]], thres_num, z7[[3]]) #result of elastic net
    return(list(z13, z5))
  }
  else if(method == "Zscore"){
    z3 = mirna_mrna_data(z1[[1]], z2[[1]], z1[[2]], z2[[2]], z2[[3]], z2[[4]], 1)
    z5 = MMI_location(Vset, z3[[3]], z3[[5]])
    z7 = mirna_mrna_loc(z3[[3]], z3[[4]], score_list[[1]], score_list[[2]], score_list[[3]], z3[[6]])
    z13 = Zscore_cal(z3[[1]], z3[[2]], thres_num, z7[[3]]) #result of Zscore score
    return(list(z13, z5))
  }
}

validation_summary = function(thres_num, Vset){
  miracle_all = result_summary(name_cancer, mirna_cancer, mrna_cancer, CWCS_all, thres_num, Vset, "miracle")
  miracle_conserved = result_summary(name_cancer, mirna_cancer, mrna_cancer, CWCS_conserved, thres_num, Vset, "miracle")
  promise_conserved = result_summary(name_cancer, mirna_cancer, mrna_cancer, qMRE_conserved, thres_num, Vset, "promise")
  targetscan_CWCS_conserved = result_summary(name_cancer, mirna_cancer, mrna_cancer, CWCS_conserved, thres_num, Vset, "sequence")
  microT_all = result_summary(name_cancer, mirna_cancer, mrna_cancer, microT, thres_num, Vset, "sequence")

  pearson = result_summary(name_cancer, mirna_cancer, mrna_cancer, combine, thres_num, Vset, "pearson")
  lasso = result_summary(name_cancer, mirna_cancer, mrna_cancer, combine, thres_num, Vset, "lasso")
  elastic = result_summary(name_cancer, mirna_cancer, mrna_cancer, combine, thres_num, Vset, "elastic")
  Zscore = result_summary(name_cancer, mirna_cancer, mrna_cancer, combine, thres_num, Vset, "Zscore")
  
  z20_1 = compare_result(miracle_all[[1]], miracle_all[[2]], thres_num)
  z20_1_con = compare_result(miracle_conserved[[1]], miracle_conserved[[2]], thres_num)
  z20_4 = compare_result(promise_conserved[[1]][[1]], promise_conserved[[2]], thres_num)
  z20_6 = compare_result(targetscan_CWCS_conserved[[1]][[1]], targetscan_CWCS_conserved[[2]], thres_num)
  z20_8 = compare_result(pearson[[1]][[1]], pearson[[2]], thres_num)
  z20_9 = compare_result(lasso[[1]][[1]], lasso[[2]], thres_num)
  z20_10 = compare_result(elastic[[1]][[1]], elastic[[2]], thres_num)
  z20_11 = compare_result(Zscore[[1]][[1]], Zscore[[2]], thres_num)
  z20_12 = compare_result(microT_all[[1]][[1]], microT_all[[2]], thres_num)
  
  #count_result and rate_result are the final result
  count_result = rbind(z20_1[[1]], z20_1_con[[1]], z20_4[[1]], z20_6[[1]], z20_8[[1]], 
                       z20_9[[1]], z20_10[[1]], z20_11[[1]], z20_12[[1]])
  rownames(count_result) = c('miRACLe', 'miRACLe_cons', 'ProMISe_mrna_cons', 'TargetScan7_CWCS_cons', 'PearsonCorrelation', 
                             'LASSO', 'Elastic-net', 'Zscore', 'DIANA-microT-CDS')
  colnames(count_result) = seq(100, thres_num, 100)
  rate_result = rbind(z20_1[[2]], z20_1_con[[2]], z20_4[[2]], z20_6[[2]], z20_8[[2]], 
                      z20_9[[2]], z20_10[[2]], z20_11[[2]], z20_12[[2]])
  rownames(rate_result) = c('miRACLe', 'miRACLe_cons', 'ProMISe_mrna_cons', 'TargetScan7_CWCS_cons', 'PearsonCorrelation', 
                            'LASSO', 'Elastic-net', 'Zscore', 'DIANA-microT-CDS')
  colnames(rate_result) = seq(100, thres_num, 100)
  return(list(count_result, rate_result))
}

wil_test = function(valid1){
  valid = matrix(data = 0, nrow = nrow(valid1), ncol = ncol(valid1))
  for(i in 1:nrow(valid)){
    valid[i, 1] = valid1[i, 1]
    for(j in 2:ncol(valid)){
      valid[i, j] = valid1[i, j] * j - valid1[i, (j - 1)] * (j - 1)
    }
  }
  
  result = matrix(data = 0, ncol = 2, nrow = 7)
  for(i in 3:9){
    result[(i - 2), 1] = wilcox.test(as.vector(valid[1, ]), as.vector(valid[i, ]), alternative = "greater", paired = TRUE, correct = FALSE)$p.value
    result[(i - 2), 2] = wilcox.test(as.vector(valid[2, ]), as.vector(valid[i, ]), alternative = "greater", paired = TRUE, correct = FALSE)$p.value
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

Vset1 = read.table("Vset_all.txt", head = TRUE, sep = "\t")
Vset2 = read.table("Vset_hc.txt", head = TRUE, sep = "\t")

name_cancer = as.matrix(read.table("TCGA_HNSC_sampleMatch.txt", head = TRUE, sep = "\t"))
mirna_cancer = as.matrix(read.table("HNSC_miRNA_expression.txt", head = FALSE, sep = "\t"))
mrna_cancer = as.matrix(read.table("HNSC_mRNA_expression.txt", head = FALSE, sep = "\t"))

###################################################################################################
######################################### MAIN code ###############################################
###################################################################################################

thres_num = 5000
V1_result = validation_summary(thres_num, Vset1)
V2_result = validation_summary(thres_num, Vset2)
test_V1 = wil_test(V1_result[[1]])
test_V2 = wil_test(V2_result[[1]])

write.table(V1_result[[2]], file = "Precision result - Vall.txt", quote = FALSE, sep = "\t")
write.table(V2_result[[2]], file = "Precision result - Vhc.txt", quote = FALSE, sep = "\t")
write.table(test_V1, file = "p-value of Wilcoxon test - Vall.txt", quote = FALSE, sep = "\t")
write.table(test_V2, file = "p-value of Wilcoxon test - Vhc.txt", quote = FALSE, sep = "\t")
