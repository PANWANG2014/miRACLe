###################################################################################################
## This Rscript generates benchmarking evaluation of Hela additional - precision analysis

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

promise_score_mRNA = function(miRNA_matrix, mRNA_matrix, matrix_set){
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
  return(list(miRNA_list, mRNA_list, score_mrna))
}

promise_score_joint = function(miRNA_matrix, mRNA_matrix, matrix_set){
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
  score_full = rs$p.xz
  return(list(miRNA_list, mRNA_list, score_full))
}

top_ranker = function(matrix, thres_num){
  miRNA_list = matrix[[1]]
  mRNA_list = matrix[[2]]
  score_matrix = matrix[[3]]
  order_agg = order(score_matrix, decreasing = TRUE)[1 : thres_num]
  sort_agg = sort(score_matrix, decreasing = TRUE)[1 : thres_num]
  miRNA_pos = ceiling(order_agg / length(mRNA_list))
  mRNA_pos = order_agg - length(mRNA_list) * (miRNA_pos - 1)
  top_infor = matrix(data = 0, nrow = thres_num, ncol = 4)
  colnames(top_infor) = c("rank", "miRNA name", "Gene ID", "score")
  for(i in 1:thres_num){
    top_infor[i, 1] = i
    top_infor[i, 2] = miRNA_list[miRNA_pos[i]]
    top_infor[i, 3] = mRNA_list[mRNA_pos[i]]
    top_infor[i, 4] = sort_agg[i]
  }
  return(top_infor)
}

accu_pred = function(matrix, pos_list, thres_num){
  miRNA_list = matrix[[1]]
  mRNA_list = matrix[[2]]
  score_matrix = matrix[[3]]
  temp1_pos = match(as.numeric(pos_list[, 1]), mRNA_list)
  temp2_pos = match(pos_list[, 2], miRNA_list)
  temp3_pos = intersect(which(is.na(temp1_pos) == FALSE), which(is.na(temp2_pos) == FALSE))
  pos_list_mod = pos_list[temp3_pos, ]
  temp4_pos = match(as.numeric(pos_list_mod[, 1]), mRNA_list)
  temp5_pos = match(pos_list_mod[, 2], miRNA_list)
  positive_position = (temp5_pos - 1) * length(mRNA_list) + temp4_pos
  
  num = thres_num/100
  pred_matrix = rep(0, num)
  pred_matrix_percent = rep(0, num)
  
  for(i in 1: num){
    ord = order(score_matrix, decreasing = TRUE)[1 : (100 * i)]
    pred_matrix[i] = sum(ord %in% positive_position)
    pred_matrix_percent[i] = pred_matrix[i]/(100 * i)
  }
  return(list(pred_matrix, pred_matrix_percent, length(positive_position)))
}

filter_exp = function(seq_matrix, miRNA_matrix, mRNA_matrix){
  temp1 = which(seq_matrix[[1]] %in% miRNA_matrix[, 1])
  temp2 = which(as.numeric(seq_matrix[[2]]) %in% as.numeric(mRNA_matrix[, 1]))
  filter_miRNA = seq_matrix[[1]][temp1]
  filter_mRNA = seq_matrix[[2]][temp2]
  filter_matrix = seq_matrix[[3]][temp2, temp1]
  return(list(filter_miRNA, filter_mRNA, filter_matrix))
}

wil_test = function(valid1){
  valid = matrix(data = 0, nrow = nrow(valid1), ncol = ncol(valid1))
  for(i in 1:ncol(valid)){
    valid[1, i] = valid1[1, i]
    for(j in 2:nrow(valid)){
      valid[j, i] = valid1[j, i] * j - valid1[(j - 1), i] * (j - 1)
    }
  }
  
  a1 = wilcox.test(as.vector(valid[, 1]), as.vector(valid[, 4]), alternative = "greater", paired = TRUE, correct = FALSE)$p.value
  a2 = wilcox.test(as.vector(valid[, 2]), as.vector(valid[, 5]), alternative = "greater", paired = TRUE, correct = FALSE)$p.value
  a3 = wilcox.test(as.vector(valid[, 3]), as.vector(valid[, 6]), alternative = "greater", paired = TRUE, correct = FALSE)$p.value
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
pos_list = as.matrix(read.table("Vset_HeLa.txt", head = TRUE, sep = "\t"))

miRDB = as.matrix(read.table("MirTarget4.txt", head = TRUE, sep = "\t"))
microT = as.matrix(read.table("DIANA_microT_CDS.txt", head = TRUE, sep = "\t"))
miRanda = as.matrix(read.table("miRanda_mirSVR.txt", head = TRUE, sep = "\t"))

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

num = 5000

miracle_miRDB_accu = accu_pred(miracle_miRDB_result, pos_list, num)
miracle_microT_accu = accu_pred(miracle_microT_result, pos_list, num)
miracle_miRanda_accu = accu_pred(miracle_miRanda_result, pos_list, num)

miRDB_accu = accu_pred(miRDB_matrix, pos_list, num)
microT_accu = accu_pred(microT_matrix, pos_list, num)
miRanda_accu = accu_pred(miRanda_matrix, pos_list, num)

percent_accu = cbind(miracle_miRDB_accu[[2]], miracle_microT_accu[[2]], miracle_miRanda_accu[[2]], miRDB_accu[[2]], 
                     microT_accu[[2]], miRanda_accu[[2]])
rownames(percent_accu) = seq(100, num, 100)
colnames(percent_accu) = c('miRACLe_MirTarget4', 'miRACLe_DIANA-microT-CDS', 'miRACLe_miRanda-mirSVR', 
                           'MirTarget4', 'DIANA-microT-CDS', 'miRanda-mirSVR')

test_result = wil_test(percent_accu)

write.table(percent_accu, file = "microarray - precision.txt", quote = FALSE, sep = "\t")
write.table(test_result, file = "microarray - p-value of Wilcoxon test.txt", quote = FALSE, sep = "\t", row.names = FALSE)

########################## RNA-seq ##############################################################

mRNA_matrix = as.matrix(read.table("HeLa_mRNA_Seq.txt", head = TRUE, sep = "\t"))
miRNA_matrix = as.matrix(read.table("HeLa_miRNA_Seq.txt", head = TRUE, sep = "\t"))


miRDB_matrix = filter_exp(miRDB_matrix1, miRNA_matrix, mRNA_matrix)
microT_matrix = filter_exp(microT_matrix1, miRNA_matrix, mRNA_matrix)
miRanda_matrix = filter_exp(miRanda_matrix1, miRNA_matrix, mRNA_matrix)

miracle_miRDB_result = miracle_score(miRNA_matrix, mRNA_matrix, miRDB_matrix)
miracle_microT_result = miracle_score(miRNA_matrix, mRNA_matrix, microT_matrix)
miracle_miRanda_result = miracle_score(miRNA_matrix, mRNA_matrix, miRanda_matrix)

num = 5000

miracle_miRDB_accu = accu_pred(miracle_miRDB_result, pos_list, num)
miracle_microT_accu = accu_pred(miracle_microT_result, pos_list, num)
miracle_miRanda_accu = accu_pred(miracle_miRanda_result, pos_list, num)

miRDB_accu = accu_pred(miRDB_matrix, pos_list, num)
microT_accu = accu_pred(microT_matrix, pos_list, num)
miRanda_accu = accu_pred(miRanda_matrix, pos_list, num)

percent_accu = cbind(miracle_miRDB_accu[[2]], miracle_microT_accu[[2]], miracle_miRanda_accu[[2]], miRDB_accu[[2]], 
                     microT_accu[[2]], miRanda_accu[[2]])

rownames(percent_accu) = seq(100, num, 100)
colnames(percent_accu) = c('miRACLe_MirTarget4', 'miRACLe_DIANA-microT-CDS', 'miRACLe_miRanda-mirSVR', 
                           'MirTarget4', 'DIANA-microT-CDS', 'miRanda-mirSVR')

test_result = wil_test(percent_accu)

write.table(percent_accu, file = "RNA-seq - precision.txt", quote = FALSE, sep = "\t")
write.table(test_result, file = "RNA-seq - p-value of Wilcoxon test.txt", quote = FALSE, sep = "\t", row.names = FALSE)

