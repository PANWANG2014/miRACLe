######################################################################################################

### Part I: FUNCTIONS THAT ARE USED IN DATA PROCESSING AND MODELING ###

#The following three functions are used to match the names of samples, mirnas and mrnas
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
#MMI_location is used to get the position of the validation set in the predicted matrix
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

####the following codes are used to calculate the rank result of different methods
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
#miracle_sam is used to get the position of individual-level result of miracle algorithm
miracle_sam = function(mirna_use1, mrna_use1, CWCS, sumup, thres_num, vset){
  num = thres_num / 100
  mirna_use = matrix(as.numeric(mirna_use1), nrow = nrow(mirna_use1))
  mrna_use = matrix(as.numeric(mrna_use1), nrow = nrow(mrna_use1))
  pro_matrix_samp = matrix(data = 0, ncol = ncol(mrna_use), nrow = thres_num, byrow = TRUE)
  match_matrix_samp = matrix(data = 0, ncol = ncol(mrna_use), nrow = thres_num, byrow = TRUE)
  match_samp = matrix(data = 0, ncol = ncol(mrna_use), nrow = num, byrow = TRUE)
  
  for(m in 1 : ncol(mirna_use)){
    pro_matrix_temp = (mirna_use[, m] %*% t(mrna_use[, m])) * CWCS/sumup[m]
    pro_matrix_samp[, m] = order(pro_matrix_temp, decreasing = TRUE)[1 : thres_num]
    match_matrix_samp[, m] = pro_matrix_samp[, m] %in% vset
  }
  for(m in 1: ncol(mirna_use)){
    for(k in 1 : num){
      match_samp[k, m] = sum(match_matrix_samp[(1 : (100 * k)), m]) / (100 * k)
    }
  }
  return(match_samp)
}
#promise_sam is used to get the position of individual-level result of promise algorithm
promise_sam = function(mirna_use1, mrna_use1, CWCS, thres_num, vset){
  num = thres_num / 100
  mirna_use = matrix(as.numeric(mirna_use1), nrow = nrow(mirna_use1))
  mrna_use = matrix(as.numeric(mrna_use1), nrow = nrow(mrna_use1))
  pro_samp_mrna = matrix(data = 0, ncol = ncol(mrna_use), nrow = thres_num, byrow = TRUE)
  pro_samp_full = matrix(data = 0, ncol = ncol(mrna_use), nrow = thres_num, byrow = TRUE)
  match_matrix_samp_mrna = matrix(data = 0, ncol = ncol(mrna_use), nrow = thres_num, byrow = TRUE)
  match_matrix_samp_full = matrix(data = 0, ncol = ncol(mrna_use), nrow = thres_num, byrow = TRUE)
  match_samp_mrna = matrix(data = 0, ncol = ncol(mrna_use), nrow = num, byrow = TRUE)
  match_samp_full = matrix(data = 0, ncol = ncol(mrna_use), nrow = num, byrow = TRUE)
  CWCS_new = t(CWCS)
  
  for(m in 1:ncol(mirna_use)){
    x = matrix(mrna_use[, m])
    rownames(x) = rownames(CWCS_new)
    z = matrix(mirna_use[, m])
    rownames(z) = colnames(CWCS_new)
    rs = roleswitch(x, z, CWCS_new)
    pro_matrix_mrna = t(rs$p.x)
    pro_matrix_full = t(rs$p.xz)
    pro_samp_mrna[, m] = order(pro_matrix_mrna, decreasing = TRUE)[1 : thres_num]
    pro_samp_full[, m] = order(pro_matrix_full, decreasing = TRUE)[1 : thres_num]
    match_matrix_samp_mrna[, m] = pro_samp_mrna[, m] %in% vset
    match_matrix_samp_full[, m] = pro_samp_full[, m] %in% vset
  }
  for(m in 1: ncol(mirna_use)){
    for(k in 1 : num){
      match_samp_mrna[k, m] = sum(match_matrix_samp_mrna[(1 : (100 * k)), m]) / (100 * k)
      match_samp_full[k, m] = sum(match_matrix_samp_full[(1 : (100 * k)), m]) / (100 * k)
    }
  }
  return(list(match_samp_mrna, match_samp_full))
}
#target_sam is used to get the position of individual-level result of targetscan
targetscan_sam = function(CWCS, thres_num, vset){
  num = thres_num / 100
  count_targetscan = rep(0, num)
  temp = order(CWCS, decreasing = TRUE)[1 : thres_num]
  match_targetscan = temp %in% vset
  for(k in 1 : num){
    count_targetscan[k] = sum(match_targetscan[1:(100 * k)]) / (100 * k)
  }
  return(count_targetscan)
}

####################################################################################################

### Part II: DATA INPUT ###

#the following package is used for promise method and expression-based methods
library(Roleswitch)
#first read in the files we need to use, which contains:
#input for miracle algorithm: mirna expression, mrna expression, CWCS_matrix, mirna list, mrna list
#input for promise algorithm: mirna expression, mrna expression, qMRE_matrix, mirna list, mrna list
#input for validation: Vset
mrna = as.matrix(read.table("mrna_list.txt", head = TRUE, sep = "\t"))
mirna = as.matrix(read.table("mirna_list.txt", head = TRUE, sep = "\t"))
qMRE_matrix = as.matrix(read.table("qMRE_matrix.txt", head = TRUE, sep = "\t"))
qMRE_conserved_matrix = as.matrix(read.table("qMRE_con_matrix.txt", head = TRUE, sep = "\t"))
CWCS_matrix1 = as.matrix(read.table("CWCS_matrix.txt", head=TRUE, sep = "\t"))
CWCS_matrix = abs(CWCS_matrix1)
Vset = read.table("V1.txt", head = TRUE, sep = "\t")
name_NCI = as.matrix(read.table("NCI60_miRNA_mRNA_paired.txt", head = TRUE, sep = "\t"))
mirna_NCI = as.matrix(read.table("NCI60_miRNA_expression.txt", head = FALSE, sep = "\t"))
mrna_NCI = as.matrix(read.table("NCI60_mRNA_expression.txt", head = FALSE, sep = "\t"))

####################################################################################################

### Part III: MAIN PROGRAM ###

#the following code is used for preparation of the analysis
thres_num = 5000
x = name_func(name_NCI, mirna_NCI, mrna_NCI)
z1 = mirna_matrix(name_NCI, mirna_NCI, x[[1]], mirna)
z2 = mrna_matrix(name_NCI, mrna_NCI, x[[2]], mrna)
z3 = mirna_mrna_data(z1[[1]], z2[[1]], z1[[2]], z2[[2]], z2[[3]], z2[[4]], 1)
z4 = mirna_mrna_data(z1[[1]], z2[[1]], z1[[2]], z2[[2]], z2[[3]], z2[[4]], 0.8)
z5 = MMI_location(Vset, z3[[3]], z3[[5]]) #this is used to locate the validation set
z6 = mirna_mrna_loc(z3[[3]], z3[[4]], mirna, mrna, qMRE_matrix, z3[[6]])
z6_con = mirna_mrna_loc(z3[[3]], z3[[4]], mirna, mrna, qMRE_conserved_matrix, z3[[6]])
z7 = mirna_mrna_loc(z3[[3]], z3[[4]], mirna, mrna, CWCS_matrix, z3[[6]])
z6_tar = mirna_mrna_loc(z4[[3]], z4[[4]], mirna, mrna, qMRE_matrix, z4[[6]])
z7_tar = mirna_mrna_loc(z4[[3]], z4[[4]], mirna, mrna, CWCS_matrix, z4[[6]])

z8 = sum_miracle(z3[[1]], z3[[2]], z7[[3]]) 
z9 = miracle_sam(z3[[1]], z3[[2]], z7[[3]], z8,  thres_num, z5) #result of miracle
z10 = promise_sam(z3[[1]], z3[[2]], z6[[3]], thres_num, z5) #result of promise all
z11 = promise_sam(z3[[1]], z3[[2]], z6_con[[3]], thres_num, z5) #result of promise
z12 = targetscan_sam(z6_tar[[3]], thres_num, z5) #get the rankers of targetscan_qMRE
z13 = targetscan_sam(z7_tar[[3]], thres_num, z5) #get the rankers of targetscan_CWCS

write.table(z9, file = "Validation rate for miracle.txt", quote = FALSE,sep = "\t")
write.table(z10[[1]], file = "Validation rate for promise mrna all.txt", quote = FALSE, sep = "\t")
write.table(z10[[2]], file = "Validation rate for promise joint all.txt", quote = FALSE, sep = "\t")
write.table(z11[[1]], file = "Validation rate for promise mrna.txt", quote = FALSE, sep = "\t")
write.table(z11[[2]], file = "Validation rate for promise joint.txt", quote = FALSE, sep = "\t")
write.table(z12, file = "Validation rate for targetscan qMRE.txt", quote = FALSE, sep = "\t")
write.table(z13, file = "Validation rate for targetscan CWCS.txt", quote = FALSE, sep = "\t")
