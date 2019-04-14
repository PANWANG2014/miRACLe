######################################################################################################

   ### Part I: FUNCTIONS THAT ARE USED IN DATA PROCESSING AND MODELING ###

#the following three functions are used to match the names of samples, mirnas and mrnas
#By using these three functions, we can get the matched mrna expression data and mirna expression data
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
#express_pre is used to prepare the matched expression data without restrict of targetscan
express_pre = function(mirna_base, mrna_base, name_base, cutoff){
  mirna_sample_name = mirna_base[1, -1]
  mrna_sample_name = mrna_base[1, -1]
  x = match(name_base[, "miRNA"], mirna_sample_name)
  y = match(name_base[, "mRNA"], mrna_sample_name)
  mirna_name = mirna_base[-1, 1]
  mrna_fullname = mrna_base[-1, 1]
  mirna_use1 = mirna_base[-1, -1]
  mrna_use1 = mrna_base[-1, -1]
  mirna_use1[is.na(mirna_use1)] = 0
  mrna_use1[is.na(mrna_use1)] = 0
  mirna_use1[mirna_use1 < 0] = 0
  mrna_use1[mrna_use1 < 0] = 0
  mirna_use2 = mirna_use1[, x]
  mrna_use2 = mrna_use1[, y]
  mirna_use = matrix(as.numeric(mirna_use2), nrow = nrow(mirna_use2))
  mrna_use = matrix(as.numeric(mrna_use2), nrow = nrow(mrna_use2))
  mirna_sign = c()
  for(i in 1:nrow(mirna_use)){
    if(sum(mirna_use[i,] == 0) >= cutoff * ncol(mirna_use)){
      mirna_sign = c(mirna_sign, i)
    }
  }
  mrna_sign = c()
  for(i in 1:nrow(mrna_use)){
    if(sum(mrna_use[i,] == 0) >= cutoff * ncol(mrna_use)){
      mrna_sign = c(mrna_sign, i)
    }
  }
  if(length(mrna_sign) != 0){
    mrna_use = mrna_use[-mrna_sign,]
    mrna_fullname = mrna_fullname[-mrna_sign]
  }
  if(length(mirna_sign) != 0){
    mirna_use = mirna_use[-mirna_sign,]
    mirna_name = mirna_name[-mirna_sign]
  }
  for(i in 1:length(mirna_name)){
    mirna_name[i] = unlist(strsplit(mirna_name[i], "\\|"))[1]
  }
  rownames(mirna_use) = mirna_name
  rownames(mrna_use) = mrna_fullname
  mrna_name = rep(0, length(mrna_fullname))
  mrna_name_all = matrix(data = 0, ncol = 2, nrow = nrow(mrna_base), byrow = TRUE)
  for(i in 1 : length(mrna_fullname)){
    mrna_name[i] = unlist(strsplit(mrna_fullname[i], "\\|"))[1]
    mrna_name_all[i, 1] = unlist(strsplit(mrna_fullname[i], "\\|"))[1]
    mrna_name_all[i, 2] = unlist(strsplit(mrna_fullname[i], "\\|"))[2]
  }
  return(list(mirna_use, mrna_use, mirna_name, mrna_name, mrna_name_all, mrna_fullname))
}
#mirna_mrna_data is used to select the RNAs that are expressed in a certain percentage of samples
mirna_mrna_data = function(mirna_use, mrna_use, mirna_name, mrna_name, mrna_name_sp, mrna_fullname, cutoff){
  mirna_use[is.na(mirna_use)] = 0
  mrna_use[is.na(mrna_use)] = 0
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
#MMI_location is used to match the validation set with the expression to get the location
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
#mirna_mrna_loc is used to reshape the CWCS (or qMRE) matrix to fit the expression data
mirna_mrna_loc = function(mirna_name, mrna_name, mirna, mrna, CWCS, mrna_fullname){
  mirna_loc = match(mirna_name, mirna)
  mrna_loc = match(mrna_name, mrna)
  CWCS = t(CWCS)
  CWCS_use = CWCS[mirna_loc, mrna_loc]
  rownames(CWCS_use) = mirna_name
  colnames(CWCS_use) = mrna_fullname
  return(list(mirna_loc, mrna_loc, CWCS_use))
}
#sum_miracle is used to calculate the normalization coefficient of miracle algorithm
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

####the following codes are used to calculate the rank result of different methods
#miracle_pop is used to get the position of population-level result of miracle algorithm
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
#promise_pop is used to get the position of population-level result of promise algorithm
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
#pearson_cal is used to get the result of pearson correlation
pearson_cal = function(mirna_use1, mrna_use1, thres_num){
  mirna_use = matrix(as.numeric(mirna_use1),nrow = nrow(mirna_use1))
  mrna_use = matrix(as.numeric(mrna_use1),nrow = nrow(mrna_use1))
  corMat = cor(t(mirna_use), t(mrna_use))
  sort_one_rank = order(corMat, decreasing=FALSE)[1:thres_num]
  sort_one_prob = sort(corMat, decreasing=FALSE)[1:thres_num]
  return(list(sort_one_rank, sort_one_prob, corMat))
}
#LASSO_cal is used to get the result of LASSO
LASSO_cal = function(mirna_use1, mrna_use1, thres_num){
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
  sort_one_rank = order(res, decreasing = FALSE)[1 : thres_num]
  sort_one_prob = sort(res, decreasing = FALSE)[1 : thres_num]
  return(list(sort_one_rank, sort_one_prob, res))
}
#ELASTIC_cal is used to get the result of elastic_net
ELASTIC_cal = function(mirna_use1, mrna_use1, thres_num){ 
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
  sort_one_rank = order(res, decreasing = FALSE)[1 : thres_num]
  sort_one_prob = sort(res, decreasing = FALSE)[1 : thres_num]
  return(list(sort_one_rank, sort_one_prob, res))
}
#zscore_cal is use to get the result of zscore
Zscore_cal = function(mirna_use1, mrna_use1, thres_num){
  mirna_use = matrix(as.numeric(mirna_use1),nrow=nrow(mirna_use1))
  mrna_use = matrix(as.numeric(mrna_use1),nrow=nrow(mrna_use1))
  mirna_use = t(scale(t(mirna_use)))
  mrna_use = t(scale(t(mrna_use)))
  
  res = matrix(data=0,ncol=nrow(mrna_use),nrow=nrow(mirna_use),byrow=TRUE)
  #zAB=(xBminA-meanB)/sdB
  for (i in 1:nrow(mirna_use)){
    for (j in 1:nrow(mrna_use)){
      indexminA = which(mirna_use[i, ] == min(mirna_use[i, ]))
      xBminA = mrna_use[j, indexminA]
      res[i,j] = abs(median(xBminA))
    }
  }
  sort_one_rank = order(res, decreasing = TRUE)[1 : thres_num]
  sort_one_prob = sort(res, decreasing = TRUE)[1 : thres_num]
  return(list(sort_one_rank, sort_one_prob, res))
}
#targetscan_rank is used to get the result of targetscan
targetscan_rank = function(CWCS, thres_num){
  target_rank = order(CWCS, decreasing = TRUE)[1 : thres_num]
  target_prob = sort(CWCS, decreasing = TRUE)[1 : thres_num]
  return(list(target_rank, target_prob))
}
#compare_result is used to get the comparison result of different methods based on the validation set
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
#result_summary is used to get the summary table for different validation sets
result_summary = function(z5, z5_tar, z5_pre){
  z20_1 = compare_result(z9, z5, thres_num)
  z20_2 = compare_result(z10[[1]], z5, thres_num)
  z20_3 = compare_result(z10[[2]], z5, thres_num)
  z20_4 = compare_result(z10_con[[1]], z5, thres_num)
  z20_5 = compare_result(z10_con[[2]], z5, thres_num)
  z20_6 = compare_result(z15[[1]], z5_tar, thres_num)
  z20_7 = compare_result(z16[[1]], z5_tar, thres_num)
  z20_8 = compare_result(z11[[1]], z5_pre, thres_num)
  z20_9 = compare_result(z12[[1]], z5_pre, thres_num)
  z20_10 = compare_result(z13[[1]], z5_pre, thres_num)
  z20_11 = compare_result(z14[[1]], z5_pre, thres_num)
  
  #count_result and rate_result are the final result
  count_result = rbind(z20_1[[1]], z20_2[[1]], z20_3[[1]], z20_4[[1]], z20_5[[1]], z20_6[[1]],
                          z20_7[[1]], z20_8[[1]], z20_9[[1]], z20_10[[1]], z20_11[[1]])
  rownames(count_result) = c("miracle", "promise mrna all", "promise joint all", "promise mrna","promise joint", 
                                "targetscan qMRE", "targetscan CWCS", "pearson", "LASSO", "elastic net", "Zscore")
  colnames(count_result) = seq(100, 5000, 100)
  rate_result = rbind(z20_1[[2]], z20_2[[2]], z20_3[[2]], z20_4[[2]], z20_5[[2]], z20_6[[2]],
                           z20_7[[2]], z20_8[[2]], z20_9[[2]], z20_10[[2]], z20_11[[2]])
  rownames(rate_result) = c("miracle", "promise mrna all", "promise joint all", "promise mrna","promise joint", 
                                 "targetscan qMRE", "targetscan CWCS", "pearson", "LASSO", "elastic net", "Zscore")
  colnames(rate_result) = seq(100, 5000, 100)
  return(list(count_result, rate_result))
}

####################################################################################################

   ### Part II: DATA INPUT ###

#the following package is used for promise method and expression-based methods
library(Roleswitch)
library(glmnet)
#first read in the files we need to use, which contains:
#input for miracle algorithm: mirna expression, mrna expression, CWCS_matrix, mirna list, mrna list
#input for promise algorithm: qMRE_matrix
#input for validation: Vset
mrna = as.matrix(read.table("mrna_list.txt", head = TRUE, sep = "\t"))
mirna = as.matrix(read.table("mirna_list.txt", head = TRUE, sep = "\t"))
qMRE_matrix = as.matrix(read.table("qMRE_matrix.txt", head = TRUE, sep = "\t"))
qMRE_conserved_matrix = as.matrix(read.table("qMRE_con_matrix.txt", head = TRUE, sep = "\t"))
CWCS_matrix1 = as.matrix(read.table("CWCS_matrix.txt", head=TRUE, sep = "\t"))
CWCS_matrix = abs(CWCS_matrix1)
Vset1 = read.table("V1.txt", head = TRUE, sep = "\t")
Vset2 = read.table("V2.txt", head = TRUE, sep = "\t")
Vset3 = read.table("V3.txt", head = TRUE, sep = "\t")
Vset4 = read.table("V4.txt", head = TRUE, sep = "\t")
name_cancer = as.matrix(read.table("TCGA_miRNA_mRNA_paired.txt", head = TRUE, sep = "\t"))
mirna_cancer = as.matrix(read.table("TCGA_miRNA_expression.txt", head = FALSE, sep = "\t"))
mrna_cancer = as.matrix(read.table("TCGA_mRNA_expression.txt", head = FALSE, sep = "\t"))

####################################################################################################

   ### Part III: MAIN PROGRAM ###

#the following code is used for preparation of the analysis
thres_num = 5000
pre = express_pre(mirna_cancer, mrna_cancer, name_cancer, 0.8)
x = name_func(name_cancer, mirna_cancer, mrna_cancer)
z1 = mirna_matrix(name_cancer, mirna_cancer, x[[1]], mirna)
z2 = mrna_matrix(name_cancer, mrna_cancer, x[[2]], mrna)
z3 = mirna_mrna_data(z1[[1]], z2[[1]], z1[[2]], z2[[2]], z2[[3]], z2[[4]], 1)
z4 = mirna_mrna_data(z1[[1]], z2[[1]], z1[[2]], z2[[2]], z2[[3]], z2[[4]], 0.8)
z5_1 = MMI_location(Vset1, z3[[3]], z3[[5]]) #this is used to locate the validation set
z5_1_tar = MMI_location(Vset1, z4[[3]], z4[[5]]) #this is used to locate the validation set of selected targetscan
z5_1_pre = MMI_location(Vset1, pre[[3]], pre[[5]]) #this is used to locate the validation set of selected targetscan
z5_2 = MMI_location(Vset2, z3[[3]], z3[[5]])
z5_2_tar = MMI_location(Vset2, z4[[3]], z4[[5]])
z5_2_pre = MMI_location(Vset2, pre[[3]], pre[[5]])
z5_3 = MMI_location(Vset3, z3[[3]], z3[[5]])
z5_3_tar = MMI_location(Vset3, z4[[3]], z4[[5]])
z5_3_pre = MMI_location(Vset3, pre[[3]], pre[[5]])
z5_4 = MMI_location(Vset4, z3[[3]], z3[[5]])
z5_4_tar = MMI_location(Vset4, z4[[3]], z4[[5]])
z5_4_pre = MMI_location(Vset4, pre[[3]], pre[[5]])

z6 = mirna_mrna_loc(z3[[3]], z3[[4]], mirna, mrna, qMRE_matrix, z3[[6]])
z6_con = mirna_mrna_loc(z3[[3]], z3[[4]], mirna, mrna, qMRE_conserved_matrix, z3[[6]])
z7 = mirna_mrna_loc(z3[[3]], z3[[4]], mirna, mrna, CWCS_matrix, z3[[6]])
z6_tar = mirna_mrna_loc(z4[[3]], z4[[4]], mirna, mrna, qMRE_matrix, z4[[6]])
z7_tar = mirna_mrna_loc(z4[[3]], z4[[4]], mirna, mrna, CWCS_matrix, z4[[6]])

#the following code is used to calculate the network of 11 different methods
z8 = sum_miracle(z3[[1]], z3[[2]], z7[[3]])
z9 = miracle_pop(z3[[1]], z3[[2]], z7[[3]], z8,  thres_num) #get the result of miracle
z10 = promise_pop(z3[[1]], z3[[2]], z6[[3]], thres_num) #get the result of promise which uses qMRE
z10_con = promise_pop(z3[[1]], z3[[2]], z6_con[[3]], thres_num) #get the result of promise which uses conserved_qMRE
z11 = pearson_cal(pre[[1]], pre[[2]], thres_num) #get the result of Pearson correlation
z12 = LASSO_cal(pre[[1]], pre[[2]], thres_num) #get the result of LASSO
z13 = ELASTIC_cal(pre[[1]], pre[[2]], thres_num) #get the result of elastic net
z14 = Zscore_cal(pre[[1]], pre[[2]], thres_num) #get the result of Z-score
z15 = targetscan_rank(z6_tar[[3]], thres_num) #get the result of targetscan qMRE
z16 = targetscan_rank(z7_tar[[3]], thres_num) #get the result of targetscan CWCS

z20_1_1 = result_summary(z5_1, z5_1_tar, z5_1_pre)
z20_1_2 = result_summary(z5_2, z5_2_tar, z5_2_pre)
z20_1_3 = result_summary(z5_3, z5_3_tar, z5_3_pre)
z20_1_4 = result_summary(z5_4, z5_4_tar, z5_4_pre)

write.table(z20_1_1[[1]], file = "Validation count result-V1.txt", quote = FALSE, sep = "\t")
write.table(z20_1_1[[2]], file = "Validation rate result-V1.txt", quote = FALSE, sep = "\t")
write.table(z20_1_2[[1]], file = "Validation count result-V2.txt", quote = FALSE, sep = "\t")
write.table(z20_1_2[[2]], file = "Validation rate result-V2.txt", quote = FALSE, sep = "\t")
write.table(z20_1_3[[1]], file = "Validation count result-V3.txt", quote = FALSE, sep = "\t")
write.table(z20_1_3[[2]], file = "Validation rate result-V3.txt", quote = FALSE, sep = "\t")
write.table(z20_1_4[[1]], file = "Validation count result-V4.txt", quote = FALSE, sep = "\t")
write.table(z20_1_4[[2]], file = "Validation rate result-V4.txt", quote = FALSE, sep = "\t")

