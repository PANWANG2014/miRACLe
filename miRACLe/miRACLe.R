######################################################################################################

   ### Part I: FUNCTIONS THAT ARE USED IN DATA PROCESSING AND MODELING ###

#matrix_transfer is used to get the sequence matching matrix
matrix_transfer = function(input_matrix){
  mirna_name = unique(input_matrix[, 2])
  mrna_name = unique(input_matrix[, 1])
  mod_matrix = matrix(data = 0, ncol = length(mirna_name), nrow = length(mrna_name))
  pos1 = match(input_matrix[, 2], mirna_name)
  pos2 = match(input_matrix[, 1], mrna_name)
  for(i in 1:nrow(input_matrix)){
    mod_matrix[pos2[i], pos1[i]] = as.numeric(input_matrix[i, 3])
  }
  return(list(mirna_name, mrna_name, mod_matrix))
}
#The following three functions are used to match samples, mirnas and mrnas
#By using these three functions, we can get the paired mrna-mirna expression data
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
mrna_matrix = function(mrna_base, mrna_name, mrna){
  mrna_exist = rep(0, nrow(mrna_base))
  mrna_name1 = rep(0, nrow(mrna_base))
  for(i in 1 : nrow(mrna_base)){
    mrna_name1[i] = unlist(strsplit(mrna_base[i, 1], "\\|"))[1]
    if(mrna_name1[i] %in% mrna == 1){
      mrna_exist[i] = i
    }
  }
  mrna_use = mrna_base[mrna_exist, mrna_name]
  mrna_name2 = mrna_name1[mrna_exist]
  return(list(mrna_use, mrna_name2))
}
#mirna_mrna_data will select the RNAs that are expressed in a certain percentage of samples
mirna_mrna_data = function(mirna_use, mrna_use, mirna_name, mrna_name, cutoff){
  if(is.vector(mirna_use) == TRUE){
    mirna_use = as.matrix(mirna_use)
  }
  if(is.vector(mrna_use) == TRUE){
    mrna_use = as.matrix(mrna_use)
  }
  mirna_use[is.na(mirna_use)] = 0
  mrna_use[is.na(mrna_use)] = 0
  mirna_use[mirna_use < 0] = 0
  mrna_use[mrna_use < 0] = 0
  mirna_sgn = seq(1, nrow(mirna_use), 1)
  mrna_sgn = seq(1, nrow(mrna_use), 1)
  for (i in 1:nrow(mirna_use)){
    if(sum(mirna_use[i,] == 0) >= cutoff * ncol(mirna_use)){
      mirna_sgn[i] = 0
    }
  }
  for (i in 1:nrow(mrna_use)){
    if(sum(mrna_use[i,] == 0) >= cutoff * ncol(mrna_use)){
      mrna_sgn[i] = 0
    }
  }  
  mirna_use1 = mirna_use[mirna_sgn, ]
  mrna_use1 = mrna_use[mrna_sgn, ]
  mirna_name1 = mirna_name[mirna_sgn]
  mrna_name1 = mrna_name[mrna_sgn]
  return(list(mirna_use1, mrna_use1, mirna_name1, mrna_name1))
}
#mirna_mrna_loc is used to get the sequence matrix that fits the expression data
mirna_mrna_loc = function(mirna_name, mrna_name, mirna, mrna, CWCS){
  mirna_loc = match(mirna_name, mirna)
  mrna_loc = match(mrna_name, mrna)
  CWCS = t(CWCS)
  CWCS_use = CWCS[mirna_loc, mrna_loc]
  rownames(CWCS_use) = mirna_name
  colnames(CWCS_use) = mrna_name
  return(list(mirna_loc, mrna_loc, CWCS_use))
}

#the following function sum_miracle and miracle_detail are the core functions of the miracle algorithm.
#sum_miracle is used to calculate the normalization coefficient of miracle algorithm
sum_miracle = function(mirna, mrna, CWCS){
  if(is.vector(mirna) == TRUE){
    mirna = as.matrix(mirna)
  }
  if(is.vector(mrna) == TRUE){
    mrna = as.matrix(mrna)
  }
  
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
#miracle_basic is used to get the detailed information of the top rankers by miracle score, both population-level and individual-level
#Here, we offer two sets of results.
#One is default, which is the first percent pairs, the other is the number of top rankers which is set by the users
miracle_detail = function(mirna_use1, mrna_use1, CWCS, sumup, ID, mirna_name, mrna_name, topList){
  if(is.vector(mirna_use1) == TRUE){
    mirna_use1 = as.matrix(mirna_use1)
  }
  if(is.vector(mrna_use1) == TRUE){
    mrna_use1 = as.matrix(mrna_use1)
  }
  numm = floor(length(CWCS[CWCS != 0])/100)
  mirna_use = matrix(as.numeric(mirna_use1), nrow = nrow(mirna_use1))
  mrna_use = matrix(as.numeric(mrna_use1), nrow = nrow(mrna_use1))
  #pro_matrix_temp saves the probability matrix of each single sample, while pro_matrix_agg is the integrated result
  pro_matrix_temp = matrix(data = 0, ncol = nrow(mrna_use), nrow = nrow(mirna_use), byrow = TRUE)
  pro_matrix_agg = matrix(data = 0, ncol = nrow(mrna_use), nrow = nrow(mirna_use), byrow = TRUE)
  
  detail_output1 = matrix(data = 0, ncol = (4 * ncol(mirna_use)), nrow = topList, byrow = TRUE)
  detail_output2 = matrix(data = 0, ncol = (4 * ncol(mirna_use)), nrow = numm, byrow = TRUE)
  temp_name1 = c()
  temp_name2 = c()
  num_max = 0
  
  for(m in 1:ncol(mirna_use)){
    pro_matrix_temp = (mirna_use[, m] %*% t(mrna_use[, m])) * CWCS / sumup[m]
    pro_matrix_agg = pro_matrix_agg + pro_matrix_temp
    sort_one_rank1 = order(pro_matrix_temp, decreasing = TRUE)[1 : topList]
    sort_one_prob1 = sort(pro_matrix_temp, decreasing = TRUE)[1 : topList]
    temp1 = ceiling(sort_one_rank1/length(mirna_name))
    temp2 = sort_one_rank1 - (temp1 - 1) * length(mirna_name)
    for(i in 1 : topList){
      detail_output1[i, (4 * m - 3)] = i
      detail_output1[i, (4 * m - 2)] = mrna_name[temp1[i]]
      detail_output1[i, (4 * m - 1)] = mirna_name[temp2[i]]
      detail_output1[i, (4 * m)] = sort_one_prob1[i]
    }
    temp_name1 = c(temp_name1, ID[m], "mrna name", "mirna name", "miracle score")
    
    num1 = floor(length(pro_matrix_temp[pro_matrix_temp != 0]))/100
    sort_one_rank2 = order(pro_matrix_temp, decreasing = TRUE)[1 : num1]
    sort_one_prob2 = sort(pro_matrix_temp, decreasing = TRUE)[1 : num1]
    temp3 = ceiling(sort_one_rank2/length(mirna_name))
    temp4 = sort_one_rank2 - (temp3 - 1) * length(mirna_name)
    for(i in 1 : num1){
      detail_output2[i, (4 * m - 3)] = i
      detail_output2[i, (4 * m - 2)] = mrna_name[temp3[i]]
      detail_output2[i, (4 * m - 1)] = mirna_name[temp4[i]]
      detail_output2[i, (4 * m)] = sort_one_prob2[i]
    }
    temp_name2 = c(temp_name2, ID[m], "mrna name", "mirna name", "miracle score")
    num_max = max(num_max, num1)
  }
  colnames(detail_output1) = temp_name1
  colnames(detail_output2) = temp_name2
  detail_output3 = detail_output2[1 : num_max, ]
  
  num2 = floor(length(pro_matrix_agg[pro_matrix_agg != 0])/100)
  out_matrix1 = matrix(data = 0, ncol = 4, nrow = topList, byrow = TRUE)
  out_matrix2 = matrix(data = 0, ncol = 4, nrow = num2, byrow = TRUE)
  
  temp_rank1 = order(pro_matrix_agg, decreasing = TRUE)[1 : topList]
  temp_prob1 = sort(pro_matrix_agg, decreasing = TRUE)[1 : topList]
  temp_rank2 = order(pro_matrix_agg, decreasing = TRUE)[1 : num2]
  temp_prob2 = sort(pro_matrix_agg, decreasing = TRUE)[1 : num2]

  mrna1 = ceiling(temp_rank1/length(mirna_name))
  mirna1 = temp_rank1 - (mrna1 - 1) * length(mirna_name)
  mrna2 = ceiling(temp_rank2/length(mirna_name))
  mirna2 = temp_rank2 - (mrna2 - 1) * length(mirna_name)
  
  for(i in 1:nrow(out_matrix1)){
    out_matrix1[i, 1] = i
    out_matrix1[i, 2] = mrna_name[mrna1[i]]
    out_matrix1[i, 3] = mirna_name[mirna1[i]]
    out_matrix1[i, 4] = temp_prob1[i]
  }
  for(i in 1:nrow(out_matrix2)){
    out_matrix2[i, 1] = i
    out_matrix2[i, 2] = mrna_name[mrna2[i]]
    out_matrix2[i, 3] = mirna_name[mirna2[i]]
    out_matrix2[i, 4] = temp_prob2[i]
  }
  colnames(out_matrix1) = c("rank", "mrna name", "mirna name", "miracle score")
  colnames(out_matrix2) = c("rank", "mrna name", "mirna name", "miracle score")
  
  return(list(out_matrix1, out_matrix2, detail_output1, detail_output3))
}
#miracle_score is the main function
miracle_score = function(seqScore, mirtaList, mirExpr, taExpr, topList){
  input_list = matrix_transfer(seqScore)
  x = name_func(mirtaList, mirExpr, taExpr)
  z1 = mirna_matrix(mirExpr, x[[1]], input_list[[1]])
  z2 = mrna_matrix(taExpr, x[[2]], input_list[[2]])
  z3 = mirna_mrna_data(z1[[1]], z2[[1]], z1[[2]], z2[[2]], 1)
  z7 = mirna_mrna_loc(z3[[3]], z3[[4]], input_list[[1]], input_list[[2]], input_list[[3]])
  z8 = sum_miracle(z3[[1]], z3[[2]], z7[[3]])
  z9 = miracle_detail(z3[[1]], z3[[2]], z7[[3]], z8, mirtaList[,1], z3[[3]], z3[[4]], topList)
  return(z9)
}
#####################################################################################################

   ### Part II: DATA INPUT ###

#The input data contains two parts.
#The first part is the sequence-based information obtained from Targetscan, which is CWCS
seqScore = as.matrix(read.table("CWCS.txt", head = TRUE, sep = "\t"))

#The other part is expression data
#mirExpr and taExpr are the mirna and mrna expression data, whose format should be strictly satisfied
#mirtaList provides the relationship of sample names in mirna expression and sample names in mrna expression
mirtaList = as.matrix(read.table("Sample_miRNA_mRNA_paired.txt", head = TRUE, sep = "\t"))
mirExpr = as.matrix(read.table("Sample_miRNA_expression.txt", head = FALSE, sep = "\t"))
taExpr = as.matrix(read.table("Sample_mRNA_expression.txt", head = FALSE, sep = "\t"))
#topList is the number of top rankers that users are interested in.
topList = 5000

#####################################################################################################

   ### Part III: MAIN PROGRAM ###

#Our algorithm provides both population-level result and sample-level result
final_output = miracle_score(seqScore, mirtaList, mirExpr, taExpr, topList)

write.table(final_output[[1]],file = paste("Population-level top", topList, "predictions.txt", sep = " "), quote = FALSE, sep = "\t", row.names = FALSE)
write.table(final_output[[2]],file = "Population-level top 1 percent predictions.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(final_output[[3]],file = paste("Individual-level top", topList, "predictions.txt", sep = " "), quote = FALSE, sep = "\t", row.names = FALSE)
write.table(final_output[[4]],file = "Individual-level top 1 percent predictions.txt", quote = FALSE, sep = "\t", row.names = FALSE)
