###################################################################################################
## This Rscript generates benchmarking evaluation of Hela additional - median fold change

###################################################################################################
## Loading required packages
library(Roleswitch)

###################################################################################################
######################################## FUNCTIONS ################################################
###################################################################################################

pair_set_mod = function(pair_set, symbol_to_ID){
  pos = match(pair_set[, 1], symbol_to_ID[, 1])
  pos1 = na.omit(pos)
  pair_set_filter = pair_set[which(is.na(pos) == FALSE), ]
  pair_set_filter[, 1] = symbol_to_ID[pos1, 2]
  return(pair_set_filter)
}

matrix_transfer = function(input_matrix){
  mirna_name = unique(input_matrix[, 2])
  mrna_name = unique(as.numeric(input_matrix[, 1]))
  mod_matrix = matrix(data = 0, ncol = length(mirna_name), nrow = length(mrna_name))
  pos1 = match(input_matrix[, 2], mirna_name)
  pos2 = match(as.numeric(input_matrix[, 1]), mrna_name)
  for(i in 1:nrow(input_matrix)){
    mod_matrix[pos2[i], pos1[i]] = mod_matrix[pos2[i], pos1[i]] + abs(as.numeric(input_matrix[i, 3]))
  }
  return(list(mirna_name, mrna_name, mod_matrix))
}

filter_exp = function(seq_matrix, miRNA_matrix, mRNA_matrix){
  temp1 = which(seq_matrix[[1]] %in% miRNA_matrix[, 1])
  temp2 = which(as.numeric(seq_matrix[[2]]) %in% as.numeric(mRNA_matrix[, 1]))
  filter_miRNA = seq_matrix[[1]][temp1]
  filter_mRNA = seq_matrix[[2]][temp2]
  filter_matrix = seq_matrix[[3]][temp2, temp1]
  return(list(filter_miRNA, filter_mRNA, filter_matrix))
}

miracle_score = function(miRNA_matrix, mRNA_matrix, matrix_set){
  miRNA_list = intersect(miRNA_matrix[, 1], matrix_set[[1]])
  mRNA_list = intersect(as.numeric(mRNA_matrix[, 1]), as.numeric(matrix_set[[2]]))
  miRNA_exp = as.numeric(miRNA_matrix[match(miRNA_list, miRNA_matrix[, 1]), 2])
  mRNA_exp = as.numeric(mRNA_matrix[match(mRNA_list, as.numeric(mRNA_matrix[, 1])), 2])
  pos1 = match(mRNA_list, as.numeric(matrix_set[[2]]))
  pos2 = match(miRNA_list, matrix_set[[1]])
  score_matrix = matrix_set[[3]][pos1, pos2]
  miracle_matrix = (mRNA_exp %*% t(miRNA_exp)) * score_matrix
  sumup = sum(miracle_matrix)
  miracle_matrix = miracle_matrix / sumup
  return(list(miRNA_list, mRNA_list, miracle_matrix))
}

promise_score = function(miRNA_matrix, mRNA_matrix, matrix_set){
  miRNA_list = intersect(miRNA_matrix[, 1], matrix_set[[1]])
  mRNA_list = intersect(as.numeric(mRNA_matrix[, 1]), as.numeric(matrix_set[[2]]))
  miRNA_exp = as.numeric(miRNA_matrix[match(miRNA_list, miRNA_matrix[, 1]), 2])
  mRNA_exp = as.numeric(mRNA_matrix[match(mRNA_list, as.numeric(mRNA_matrix[, 1])), 2])
  pos1 = match(mRNA_list, as.numeric(matrix_set[[2]]))
  pos2 = match(miRNA_list, matrix_set[[1]])
  score_matrix = matrix_set[[3]][pos1, pos2]
  z = matrix(miRNA_exp)
  x = matrix(mRNA_exp)
  rownames(x) = mRNA_list
  rownames(z) = miRNA_list
  colnames(score_matrix) = miRNA_list
  rownames(score_matrix) = mRNA_list
  rs = roleswitch(x, z, score_matrix) 
  score_mrna = rs$p.x
  score_full = rs$p.xz
  return(list(miRNA_list, mRNA_list, score_mrna, score_full))
}

method_medfc = function(score, vset, num){
  method_score = score[[3]]
  mrna = score[[2]]
  mirna = score[[1]]
  method_score = t(method_score)
  mirna_name1 = colnames(vset)
  mrna_name = vset[, 1]
  mirna_name = gsub("\\.", "-", mirna_name1)
  mirna_poss = match(mirna_name, mirna)
  mirna_pos = mirna_poss[which(is.na(mirna_poss) != TRUE)]
  if(length(mirna_pos) == 1){
    pro_matrix_sel = method_score[mirna_pos, ]
    sel_mrna_score = sort(pro_matrix_sel, decreasing = TRUE)[1:num]
    temp = order(pro_matrix_sel, decreasing = TRUE)[1:num]
    sel_mrna_name = as.numeric(mrna[temp])
    tar_num = length(which(pro_matrix_sel != 0))
    return(list(sel_mrna_name, sel_mrna_score, tar_num))
  }
  else if(length(mirna_pos) > 1){
    pro_matrix_sel = method_score[mirna_pos, ]
    sel_mrna_name = matrix(data = 0, ncol = length(mirna_pos), nrow = num, byrow = TRUE)
    sel_mrna_score = matrix(data = 0, ncol = length(mirna_pos), nrow = num, byrow = TRUE)
    tar_num = rep(0, length(mirna_pos))
    for(i in 1: length(mirna_pos)){
      sel_mrna_score[,i] = sort(pro_matrix_sel[i, ], decreasing = TRUE)[1:num]
      temp = order(pro_matrix_sel[i, ], decreasing = TRUE)[1:num]
      sel_mrna_name[, i] = as.numeric(mrna[temp])
      tar_num[i] = length(which(pro_matrix_sel[i, ] != 0))
    }
    return(list(sel_mrna_name, sel_mrna_score, tar_num))
  }
  else{
    sel_mrna_name = c()
    sel_mrna_score = c()
    tar_num = c()
    return(list(sel_mrna_name, sel_mrna_score, tar_num))
  }
}

medfc_test = function(vset, mirna_name, sel_mrna_name1, tar_num, num){
  if(length(tar_num) == 0){
    mrna_count = c()
    mrnafc_matrix = c()
    return(list(mrna_count, mrnafc_matrix))
  }
  else{
    tar_num1 = tar_num[which(tar_num != 0)]
    thres = pmin(tar_num1, num)
    mirna_vali1 = as.vector(colnames(vset)[-1])
    mirna_vali = gsub("\\.", "-", mirna_vali1)
    vset_name = as.vector(as.numeric(vset[-1,1]))
    vset_score = vset[-1,-1]
    mirna_pos = match(mirna_vali, mirna_name)
    mirna_pos_use = mirna_pos[which(is.na(mirna_pos) != TRUE)]
    if(length(which(is.na(mirna_pos) != TRUE)) == 1){
      vset_score_use = vset_score[, which(is.na(mirna_pos) != TRUE)]
      sel_mrna_name = sel_mrna_name1
      mrnafc_matrix = matrix(data = 0, ncol = 1, nrow = length(sel_mrna_name), byrow = TRUE)
      temp = match(as.numeric(sel_mrna_name[1:thres]), as.numeric(vset_name))
      temp2 = vset_score_use[temp[which(is.na(temp) != TRUE)]]
      mrna_count = length(which(is.na(temp2) != TRUE))
      mrnafc_matrix[1:mrna_count] = temp2[which(is.na(temp2) != TRUE)]
    }
    else{
      vset_score_use2 = vset_score[, which(is.na(mirna_pos) != TRUE)]
      vset_score_use = vset_score_use2[,which(tar_num != 0)]
      sel_mrna_name = sel_mrna_name1[,which(tar_num != 0)]
      mrna_count = rep(0, ncol(vset_score_use))
      mrnafc_matrix = matrix(data = 0, ncol = ncol(vset_score_use), nrow = nrow(sel_mrna_name), byrow = TRUE)
      for(i in 1:ncol(vset_score_use)){
        temp = match(as.numeric(sel_mrna_name[1:thres[i],i]), as.numeric(vset_name))
        #mrna_count[i] = length(which(is.na(temp) != TRUE))
        temp2 = vset_score_use[temp[which(is.na(temp) != TRUE)], i]
        mrna_count[i] = length(which(is.na(temp2) != TRUE))
        mrnafc_matrix[1:mrna_count[i], i] = temp2[which(is.na(temp2) != TRUE)]
      }
    }
    return(list(mrna_count, mrnafc_matrix))
  }
}

medfc_vector = function(mrna_count, mrnafc_matrix){
  if(length(mrna_count) == 0){
    fc_vector = c()
    return(fc_vector)
  }
  else{
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
}

medfc_summary = function(vset){
  num = 800
  z10 = method_medfc(miracle_miRDB_result, vset, num) #miracle + miRDB
  z11 = method_medfc(miracle_microT_result, vset, num) #miracle + microT
  z12 = method_medfc(miracle_miRanda_result, vset, num) #miracle + miRanda
  
  z20 = method_medfc(miRDB_matrix, vset, num) #for miRDB
  z21 = method_medfc(microT_matrix, vset, num) #for microT
  z22 = method_medfc(miRanda_matrix, vset, num) #for miRanda

  z27 = medfc_test(vset, miracle_miRDB_result[[1]], z10[[1]], z10[[3]], num)
  z28 = medfc_vector(z27[[1]], z27[[2]]) #for miracle + miRDB
  z29 = medfc_test(vset, miracle_microT_result[[1]], z11[[1]], z11[[3]], num)
  z30 = medfc_vector(z29[[1]], z29[[2]]) #for miracle + microT
  z31 = medfc_test(vset, miracle_miRanda_result[[1]], z12[[1]], z12[[3]], num)
  z32 = medfc_vector(z31[[1]], z31[[2]]) #for miracle + miRanda
  
  z47 = medfc_test(vset, miRDB_matrix[[1]], z20[[1]], z20[[3]], num)
  z48 = medfc_vector(z47[[1]], z47[[2]]) #for miRDB 
  z49 = medfc_test(vset, microT_matrix[[1]], z21[[1]], z21[[3]], num)
  z50 = medfc_vector(z49[[1]], z49[[2]]) #for microT
  z51 = medfc_test(vset, miRanda_matrix[[1]], z22[[1]], z22[[3]], num)
  z52 = medfc_vector(z51[[1]], z51[[2]]) #for miranda 

  final_median = rbind(z28[1:500], z30[1:500], z32[1:500], z48[1:500], z50[1:500], z52[1:500])
  rownames(final_median) = c('miRACLe_MirTarget4', 'miRACLe_DIANA-microT-CDS', 'miRACLe_miRanda-mirSVR', 
                             'MirTarget4', 'DIANA-microT-CDS', 'miRanda-mirSVR')
  colnames(final_median) = seq(1, 500, 1)
  return(final_median)
}

wil_test = function(valid){
  a1 = wilcox.test(as.vector(valid[1, ]), as.vector(valid[4, ]), alternative = "less", paired = TRUE, correct = FALSE)$p.value
  a2 = wilcox.test(as.vector(valid[2, ]), as.vector(valid[5, ]), alternative = "less", paired = TRUE, correct = FALSE)$p.value
  a3 = wilcox.test(as.vector(valid[3, ]), as.vector(valid[6, ]), alternative = "less", paired = TRUE, correct = FALSE)$p.value
  a = c(a1, a2, a3)
  result = matrix(a, ncol = 3)
  colnames(result) = c('miRACLe_MirTarget4 vs MirTarget4', 'miRACLe_DIANA-microT-CDS vs DIANA-microT-CDS', 
                       'miRACLe_miRanda-mirSVR vs miRanda-mirSVR')
  return(result)
}

###################################################################################################
####################################### DATA INPUTS ###############################################
###################################################################################################

symbol_to_ID = as.matrix(read.table("Symbol_to_ID.txt", head = TRUE, sep = "\t"))
miRDB = as.matrix(read.table("MirTarget4.txt", head = TRUE, sep = "\t"))
microT = as.matrix(read.table("DIANA_microT_CDS.txt", head = TRUE, sep = "\t"))
miRanda = as.matrix(read.table("miRanda_mirSVR.txt", head = TRUE, sep = "\t"))

Liu = as.matrix(read.table("Transet_HeLa_Seq.txt", head = TRUE, sep = "\t"))
Selbach = as.matrix(read.table("Transet_HeLa_Array.txt", head = TRUE, sep = "\t"))

###################################################################################################
######################################### MAIN code ###############################################
###################################################################################################

miRDB_mod = pair_set_mod(miRDB, symbol_to_ID)
microT_mod = pair_set_mod(microT, symbol_to_ID)
miRanda_mod = pair_set_mod(miRanda, symbol_to_ID)

miRDB_matrix1 = matrix_transfer(miRDB_mod)
microT_matrix1 = matrix_transfer(microT_mod)
miRanda_matrix1 = matrix_transfer(miRanda_mod)

########################## microarray ##############################################################

mRNA_matrix = as.matrix(read.table("HeLa_mRNA_array.txt", head = TRUE, sep = "\t"))
miRNA_matrix = as.matrix(read.table("HeLa_miRNA_array.txt", head = TRUE, sep = "\t"))

miRDB_matrix = filter_exp(miRDB_matrix1, miRNA_matrix, mRNA_matrix)
microT_matrix = filter_exp(microT_matrix1, miRNA_matrix, mRNA_matrix)
miRanda_matrix = filter_exp(miRanda_matrix1, miRNA_matrix, mRNA_matrix)

miracle_miRDB_result = miracle_score(miRNA_matrix, mRNA_matrix, miRDB_matrix)
miracle_microT_result = miracle_score(miRNA_matrix, mRNA_matrix, microT_matrix)
miracle_miRanda_result = miracle_score(miRNA_matrix, mRNA_matrix, miRanda_matrix)

final_Liu = medfc_summary(Liu)
final_Selbach = medfc_summary(Selbach)

test_Liu = wil_test(final_Liu)
test_Selbach = wil_test(final_Selbach)

write.table(final_Liu, file = "Microarray - mfc based on seq transet.txt", quote = FALSE, sep = "\t")
write.table(final_Selbach, file = "Microarray - mfc based on array transet.txt", quote = FALSE, sep = "\t")
write.table(test_Liu, file = "Microarray - p-value of Wilcoxon test for seq transet.txt", quote = FALSE, sep = "\t")
write.table(test_Selbach, file = "Microarray - p-value of Wilcoxon test for array transet.txt", quote = FALSE, sep = "\t")

########################## RNA-Seq ##############################################################

mRNA_matrix = as.matrix(read.table("HeLa_mRNA_Seq.txt", head = TRUE, sep = "\t"))
miRNA_matrix = as.matrix(read.table("HeLa_miRNA_Seq.txt", head = TRUE, sep = "\t"))

miRDB_matrix = filter_exp(miRDB_matrix1, miRNA_matrix, mRNA_matrix)
microT_matrix = filter_exp(microT_matrix1, miRNA_matrix, mRNA_matrix)
miRanda_matrix = filter_exp(miRanda_matrix1, miRNA_matrix, mRNA_matrix)

miracle_miRDB_result = miracle_score(miRNA_matrix, mRNA_matrix, miRDB_matrix)
miracle_microT_result = miracle_score(miRNA_matrix, mRNA_matrix, microT_matrix)
miracle_miRanda_result = miracle_score(miRNA_matrix, mRNA_matrix, miRanda_matrix)

final_Liu = medfc_summary(Liu)
final_Selbach = medfc_summary(Selbach)

test_Liu = wil_test(final_Liu)
test_Selbach = wil_test(final_Selbach)

write.table(final_Liu, file = "RNA-Seq - mfc based on seq transet.txt", quote = FALSE, sep = "\t")
write.table(final_Selbach, file = "RNA-Seq - mfc based on array transet.txt", quote = FALSE, sep = "\t")
write.table(test_Liu, file = "RNA-Seq - p-value of Wilcoxon test for seq transet.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(test_Selbach, file = "RNA-Seq - p-value of Wilcoxon test for array transet.txt", quote = FALSE, sep = "\t", row.names = FALSE)
