###################################################################################################
## This Rscript generates benchmarking evaluation of MCC/NCI60 - precision analysis

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

mirna_matrix = function(mirna_base, mirna_name, mirna){  
  mirna_name1 = rep(0,nrow(mirna_base))
  for(i in 1:nrow(mirna_base)){
    mirna_name1[i] = unlist(strsplit(mirna_base[i,1],"\\|"))[1]
  }
  mirna_use = mirna_base[mirna_name1%in%mirna, mirna_name]
  mirna_name1 = mirna_name1[mirna_name1%in%mirna]
  return(list(mirna_use, mirna_name1))
}

mirna_matrix_IDA = function(mirna_base, mirna_name){
  mirna_name1 = rep(0,nrow(mirna_base))
  for(i in 1:nrow(mirna_base)){
    mirna_name1[i] = unlist(strsplit(mirna_base[i,1],"\\|"))[1]
  }
  mirna_use = mirna_base[, mirna_name]
  mirna_use = mirna_use[-1,]
  return(list(mirna_use, mirna_name1))
}

mrna_matrix = function(mrna_base, mrna_name, mrna){
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

mrna_matrix_IDA = function(mrna_base, mrna_name){
  mrna_name_all = matrix(data = 0, ncol = 2, nrow = nrow(mrna_base), byrow = TRUE)
  mrna_name1 = rep(0, nrow(mrna_base))
  for(i in 1 : nrow(mrna_base)){
    mrna_name1[i] = unlist(strsplit(mrna_base[i, 1], "\\|"))[1]
    mrna_name_all[i, 1] = unlist(strsplit(mrna_base[i, 1], "\\|"))[1]
    mrna_name_all[i, 2] = unlist(strsplit(mrna_base[i, 1], "\\|"))[2]
  }
  mrna_use = mrna_base[, mrna_name]
  mrna_use = mrna_use[-1,]
  mrna_fullname = mrna_base[, 1]
  return(list(mrna_use, mrna_name1, mrna_name_all, mrna_fullname))
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

RNA_name_sub = function(mirna_name){
  mirna = c()
  for(i in 1:length(mirna_name)){
    temp = unlist(strsplit(mirna_name[i],"\\|"))[1]
    mirna = c(mirna, temp)
  }
  return(mirna)
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
    temp2 = CWCS * temp
    res[i] = sum(temp2)
  }
  return(res)
}

miracle_pop = function(mirna_use1, mrna_use1, CWCS, sumup, thres_num){
  mirna_use = matrix(as.numeric(mirna_use1), nrow = nrow(mirna_use1))
  mrna_use = matrix(as.numeric(mrna_use1), nrow = nrow(mrna_use1))
  pro_matrix_temp = matrix(data = 0, ncol = nrow(mrna_use), nrow = nrow(mirna_use), byrow = TRUE)
  pro_matrix_agg = matrix(data = 0, ncol = nrow(mrna_use), nrow = nrow(mirna_use), byrow = TRUE)
  pos = which(sumup != 0)
  mirna_use = mirna_use[, pos]
  mrna_use = mrna_use[, pos]
  sumup = sumup[pos]
  
  for(m in 1:ncol(mirna_use)){
    pro_matrix_temp = (mirna_use[, m] %*% t(mrna_use[, m])) * CWCS / sumup[m]
    pro_matrix_agg = pro_matrix_agg + pro_matrix_temp
  }
  order_agg = order(pro_matrix_agg, decreasing = TRUE)[1 : thres_num]
  return(order_agg)
}

sIDA_cal = function(mirna_use1, mrna_use1, mirna_name, mrna_name, thres_num, score_matrix, alpha){
  mirna_use = matrix(as.numeric(mirna_use1),nrow=nrow(mirna_use1))
  mrna_use = matrix(as.numeric(mrna_use1),nrow=nrow(mrna_use1))
  
  data = cbind(t(mirna_use), t(mrna_use))
  cause = 1:nrow(mirna_use)
  effect = (1 + nrow(mirna_use)): (nrow(mrna_use) + nrow(mirna_use))
  allname = c(mirna_name, mrna_name[, 1])
  multiset = character(0)
  result = matrix(nrow = nrow(mrna_use), ncol = nrow(mirna_use))
  suffStat = list(C = cor(data), n = nrow(data))
  pcFit = pc(suffStat, indepTest = gaussCItest, p = ncol(data), alpha = alpha)
  for (l in cause){
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
  res = result * score_matrix
  sort_one_rank = order(res, decreasing=FALSE)[1:thres_num]
  sort_one_prob = sort(res, decreasing=FALSE)[1:thres_num]
  return(list(sort_one_rank, sort_one_prob, res))
}

cIDA_cal = function(mirna_use1, mrna_use1, mirna_name, mrna_name, thres_num, alpha){
  mirna_use = matrix(as.numeric(mirna_use1),nrow=nrow(mirna_use1))
  mrna_use = matrix(as.numeric(mrna_use1),nrow=nrow(mrna_use1))
  
  data = cbind(t(mirna_use), t(mrna_use))
  cause = 1:nrow(mirna_use)
  effect = (1 + nrow(mirna_use)): (nrow(mrna_use) + nrow(mirna_use))
  allname = c(mirna_name, mrna_name[, 1])
  multiset = character(0)
  result = matrix(nrow = nrow(mrna_use), ncol = nrow(mirna_use))
  suffStat = list(C = cor(data), n = nrow(data))
  pcFit = pc(suffStat, indepTest = gaussCItest, p = ncol(data), alpha = alpha)
  for (l in cause){
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
  sort_one_rank = order(result, decreasing=FALSE)[1:thres_num]
  sort_one_prob = sort(result, decreasing=FALSE)[1:thres_num]
  return(list(sort_one_rank, sort_one_prob, result))
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

result_summary = function(name_cancer, mirna_cancer, mrna_cancer, score_matrix, thres_num, Vset, diff_mirna, diff_mrna, alpha, method){
  score_list = matrix_transfer(score_matrix)
  x = name_func(name_cancer, mirna_cancer, mrna_cancer)
  z1 = mirna_matrix(mirna_cancer, x[[1]], score_list[[1]])
  z2 = mrna_matrix(mrna_cancer, x[[2]], score_list[[2]])
  if(method == "miracle"){
    z3 = mirna_mrna_data(z1[[1]], z2[[1]], z1[[2]], z2[[2]], z2[[3]], z2[[4]], 1)
    z5 = MMI_location(Vset, z3[[3]], z3[[5]])
    z7 = mirna_mrna_loc(z3[[3]], z3[[4]], score_list[[1]], score_list[[2]], score_list[[3]], z3[[6]])
    z8 = sum_miracle(z3[[1]], z3[[2]], z7[[3]])
    z9 = miracle_pop(z3[[1]], z3[[2]], z7[[3]], z8, thres_num) #get the result of miracle
    return(list(z9, z5))
  }
  else if(method == "miracleDE"){
    z3 = mirna_mrna_data(z1[[1]], z2[[1]], z1[[2]], z2[[2]], z2[[3]], z2[[4]], 1)
    temp1 = match(RNA_name_sub(diff_mirna), z3[[3]])
    diff_mirna_list = RNA_name_sub(diff_mirna)[which(is.na(temp1) != TRUE)]
    temp2 = match(RNA_name_sub(diff_mrna), z3[[4]])
    diff_mrna_list = RNA_name_sub(diff_mrna)[which(is.na(temp2) != TRUE)]
    temp3 = match(diff_mrna_list, z3[[4]])
    temp4 = match(diff_mirna_list, z3[[3]])
    
    mirna_modi = z3[[1]][temp4, ]
    mrna_modi = z3[[2]][temp3,]
    mirna_name_modi = z3[[3]][temp4]
    mrna_name_modi = z3[[4]][temp3]
    mrna_namefull_modi = z3[[6]][temp3]
    mrna_namepair_modi = z3[[5]][temp3, ]
    
    z5 = MMI_location(Vset, mirna_name_modi, mrna_namepair_modi)
    z7 = mirna_mrna_loc(mirna_name_modi, mrna_name_modi, score_list[[1]], score_list[[2]], score_list[[3]], mrna_namefull_modi)
    z8 = sum_miracle(mirna_modi, mrna_modi, z7[[3]])
    z9 = miracle_pop(mirna_modi, mrna_modi, z7[[3]], z8, thres_num) #get the result of miracle
    return(list(z9, z5))
  }
  else if(method == "sIDA"){
    z3 = mirna_mrna_data(z1[[1]], z2[[1]], z1[[2]], z2[[2]], z2[[3]], z2[[4]], 1)
    temp1 = match(RNA_name_sub(diff_mirna), z3[[3]])
    diff_mirna_list = RNA_name_sub(diff_mirna)[which(is.na(temp1) != TRUE)]
    temp2 = match(RNA_name_sub(diff_mrna), z3[[4]])
    diff_mrna_list = RNA_name_sub(diff_mrna)[which(is.na(temp2) != TRUE)]
    temp3 = match(diff_mrna_list, z3[[4]])
    temp4 = match(diff_mirna_list, z3[[3]])
    
    mirna_modi = z3[[1]][temp4, ]
    mrna_modi = z3[[2]][temp3,]
    mirna_name_modi = z3[[3]][temp4]
    mrna_name_modi = z3[[4]][temp3]
    mrna_namefull_modi = z3[[6]][temp3]
    mrna_namepair_modi = z3[[5]][temp3, ]
    
    z5 = MMI_location(Vset, mirna_name_modi, mrna_namepair_modi) #this is used to locate the validation set
    z7 = mirna_mrna_loc(mirna_name_modi, mrna_name_modi, score_list[[1]], score_list[[2]], score_list[[3]], mrna_namefull_modi)
    z9 = sIDA_cal(mirna_modi, mrna_modi, mirna_name_modi, mrna_namepair_modi, thres_num, z7[[3]], alpha) #get the result of sIDA
    return(list(z9, z5))
  }
  else if(method == "cIDA"){
    x = name_func(name_cancer, mirna_cancer, mrna_cancer)
    z1 = mirna_matrix_IDA(mirna_cancer, x[[1]])
    z2 = mrna_matrix_IDA(mrna_cancer, x[[2]])
    z3 = mirna_mrna_data(z1[[1]], z2[[1]], z1[[2]], z2[[2]], z2[[3]], z2[[4]], 1)
    temp1 = match(RNA_name_sub(diff_mirna), z3[[3]])
    diff_mirna_list = RNA_name_sub(diff_mirna)[which(is.na(temp1) != TRUE)]
    temp2 = match(RNA_name_sub(diff_mrna), z3[[4]])
    diff_mrna_list = RNA_name_sub(diff_mrna)[which(is.na(temp2) != TRUE)]
    temp3 = match(diff_mrna_list, z3[[4]])
    temp4 = match(diff_mirna_list, z3[[3]])
    
    mirna_modi = z3[[1]][temp4, ]
    mrna_modi = z3[[2]][temp3,]
    mirna_name_modi = z3[[3]][temp4]
    mrna_name_modi = z3[[4]][temp3]
    mrna_namefull_modi = z3[[6]][temp3]
    mrna_namepair_modi = z3[[5]][temp3, ]
    
    z5 = MMI_location(Vset, mirna_name_modi, mrna_namepair_modi) #this is used to locate the validation set
    z9 = cIDA_cal(mirna_modi, mrna_modi, mirna_name_modi, mrna_namepair_modi, thres_num, alpha) #get the result of classical IDA
    return(list(z9, z5))
  }
}

validation_summary = function(thres_num, alpha, Vset){
  miracle_all = result_summary(name_cancer, mirna_cancer, mrna_cancer, CWCS_all, thres_num, Vset, diff_mirna, diff_mrna, alpha, "miracle")
  miracle_conserved = result_summary(name_cancer, mirna_cancer, mrna_cancer, CWCS_conserved, thres_num, Vset, diff_mirna, diff_mrna, alpha, "miracle")
  miracleDE_all = result_summary(name_cancer, mirna_cancer, mrna_cancer, CWCS_all, thres_num, Vset, diff_mirna, diff_mrna, alpha, "miracleDE")
  miracleDE_conserved = result_summary(name_cancer, mirna_cancer, mrna_cancer, CWCS_conserved, thres_num, Vset, diff_mirna, diff_mrna, alpha, "miracleDE")
  cIDA = result_summary(name_cancer, mirna_cancer, mrna_cancer, combine, thres_num, Vset, diff_mirna, diff_mrna, alpha, "cIDA")
  sIDA = result_summary(name_cancer, mirna_cancer, mrna_cancer, combine, thres_num, Vset, diff_mirna, diff_mrna, alpha, "sIDA")

  z20_1 = compare_result(miracle_all[[1]], miracle_all[[2]], thres_num)
  z20_2 = compare_result(miracle_conserved[[1]], miracle_conserved[[2]], thres_num)
  z20_3 = compare_result(miracleDE_all[[1]], miracleDE_all[[2]], thres_num)
  z20_4 = compare_result(miracleDE_conserved[[1]], miracleDE_conserved[[2]], thres_num)
  z20_5 = compare_result(cIDA[[1]][[1]], cIDA[[2]], thres_num)
  z20_6 = compare_result(sIDA[[1]][[1]], sIDA[[2]], thres_num)
  
  count_result = rbind(z20_1[[1]], z20_2[[1]], z20_3[[1]], z20_4[[1]], z20_5[[1]], z20_6[[1]])
  rownames(count_result) = c('miRACLe', 'miRACle_cons', 'miRACLe_DE', 'miRACLe_cons_DE',
                             'cIDA', 'mIDA')
  colnames(count_result) = seq(100, thres_num, 100)
  rate_result = rbind(z20_1[[2]], z20_2[[2]], z20_3[[2]], z20_4[[2]], z20_5[[2]], z20_6[[2]])
  rownames(rate_result) = c('miRACLe', 'miRACle_cons', 'miRACLe_DE', 'miRACLe_cons_DE',
                            'cIDA', 'mIDA')
  colnames(rate_result) = seq(100, thres_num, 100)
  return(list(count_result, rate_result))
}

wil_test = function(valid){
  valid = matrix(data = 0, nrow = nrow(valid1), ncol = ncol(valid1))
  for(i in 1:nrow(valid)){
    valid[i, 1] = valid1[i, 1]
    for(j in 2:ncol(valid)){
      valid[i, j] = valid1[i, j] * j - valid1[i, (j - 1)] * (j - 1)
    }
  }
  a1 = wilcox.test(as.vector(valid[1, ]), as.vector(valid[3, ]), alternative = "greater", paired = TRUE, correct = FALSE)$p.value
  a2 = wilcox.test(as.vector(valid[2, ]), as.vector(valid[4, ]), alternative = "greater", paired = TRUE, correct = FALSE)$p.value
  a3 = wilcox.test(as.vector(valid[3, ]), as.vector(valid[5, ]), alternative = "greater", paired = TRUE, correct = FALSE)$p.value
  a4 = wilcox.test(as.vector(valid[3, ]), as.vector(valid[6, ]), alternative = "greater", paired = TRUE, correct = FALSE)$p.value
  a5 = wilcox.test(as.vector(valid[4, ]), as.vector(valid[5, ]), alternative = "greater", paired = TRUE, correct = FALSE)$p.value
  a6 = wilcox.test(as.vector(valid[4, ]), as.vector(valid[6, ]), alternative = "greater", paired = TRUE, correct = FALSE)$p.value
  a = c(a1, a2, a3, a4, a5, a6)
  result = matrix(a, ncol = 6)
  colnames(result) = c('miRACLe vs miRACLe_DE', 'miRACLe_cons vs miRACLe_cons_DE', 'miRACLe_DE vs cIDA', 'miRACLe_DE vs mIDA',
                       'miRACLe_cons_DE vs cIDA', 'miRACLe_cons_DE vs mIDA')
  return(result)
}

###################################################################################################
####################################### DATA INPUTS ###############################################
###################################################################################################

CWCS_all = as.matrix(read.table("TargetScan7_CWCS.txt", head=TRUE, sep = "\t"))
CWCS_conserved = as.matrix(read.table("TargetScan7_CWCS_cons.txt", head=TRUE, sep = "\t"))
combine = as.matrix(read.table("Combine_MMIs.txt", head = TRUE, sep = "\t"))
Vall = read.table("Vset_all.txt", head = TRUE, sep = "\t")
Vhc = read.table("Vset_hc.txt", head = TRUE, sep = "\t")

name_cancer = as.matrix(read.table("MCC_tumor_sampleMatch.txt", head = TRUE, sep = "\t"))
mirna_cancer = as.matrix(read.table("MCC_miRNA_expression.txt", head = FALSE, sep = "\t"))
mrna_cancer = as.matrix(read.table("MCC_mRNA_expression.txt", head = FALSE, sep = "\t"))
diff_mirna = as.matrix(read.table("MCC_diffmiRNAlist.txt", head = FALSE, sep = "\t"))
diff_mrna = as.matrix(read.table("MCC_diffmRNAlist.txt", head = FALSE, sep = "\t"))

# name_cancer = as.matrix(read.table("NCI60_mesenchymal_sampleMatch.txt", head = TRUE, sep = "\t"))
# mirna_cancer = as.matrix(read.table("NCI60_miRNA_expression.txt", head = FALSE, sep = "\t"))
# mrna_cancer = as.matrix(read.table("NCI60_mRNA_expression.txt", head = FALSE, sep = "\t"))
# diff_mirna = as.matrix(read.table("NCI60_diffmiRNAlist.txt", head = FALSE, sep = "\t"))
# diff_mrna = as.matrix(read.table("NCI60_diffmRNAlist.txt", head = FALSE, sep = "\t"))

###################################################################################################
######################################### MAIN code ###############################################
###################################################################################################

thres_num = 5000
alpha = 0.01

valid_Vall = validation_summary(thres_num, alpha, Vall)
valid_Vhc = validation_summary(thres_num, alpha, Vhc)
test_Vall = wil_test(valid_Vall[[1]])
test_Vhc = wil_test(valid_Vhc[[1]])

write.table(valid_Vall[[2]], file = "Precision result - Vall.txt", quote = FALSE, sep = "\t")
write.table(valid_Vhc[[2]], file = "Precision result - Vhc.txt", quote = FALSE, sep = "\t")
write.table(test_Vall, file = "Vall - p-value of Wilcoxon test.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(test_Vhc, file = "Vhc - p-value of Wilcoxon test.txt", quote = FALSE, sep = "\t", row.names = FALSE)

