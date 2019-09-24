###################################################################################################
## This Rscript generates benchmarking evaluation of HeLa - precision analysis

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
  result = matrix(data = 0, ncol = 2, nrow = 13)
  for(i in 3:15){
    result[(i - 2), 1] = wilcox.test(as.vector(valid[, 1]), as.vector(valid[, i]), alternative = "greater", paired = TRUE, correct = FALSE)$p.value
    result[(i - 2), 2] = wilcox.test(as.vector(valid[, 2]), as.vector(valid[, i]), alternative = "greater", paired = TRUE, correct = FALSE)$p.value
  }
  
  colnames(result) = c('miRACLe', 'miRACLe_cons')
  rownames(result) = c('ProMISe_mrna', 'ProMISe_joint', 'ProMISe_mrna_cons', 'ProMISe_joint_cons', 'TargetScan7_CWCS', 
                       'TargetScan7_CWCS_cons', 'PITA', 'miRTar2GO_HeLa', 'MirTarget4', 'DIANA-microT-CDS','miRanda-mirSVR', 
                       'miRmap', 'miRWalk3')
  return(result)
}

###################################################################################################
####################################### DATA INPUTS ###############################################
###################################################################################################

symbol_to_ID = as.matrix(read.table("Symbol_to_ID.txt", head = TRUE, sep = "\t"))
pos_list = as.matrix(read.table("Vset_HeLa.txt", head = TRUE, sep = "\t"))

CWCS_all = as.matrix(read.table("TargetScan7_CWCS.txt", head = TRUE, sep = "\t"))
CWCS_cons = as.matrix(read.table("TargetScan7_CWCS_cons.txt", head = TRUE, sep = "\t"))
qMRE_all = as.matrix(read.table("TargetScan7_qMRE.txt", head = TRUE, sep = "\t"))
qMRE_cons = as.matrix(read.table("TargetScan7_qMRE_cons.txt", head = TRUE, sep = "\t"))
PITA = as.matrix(read.table("PITA.txt", head = TRUE, sep = "\t"))
mirtar = as.matrix(read.table("miRTar2GO_HeLa.txt", head = TRUE, sep = "\t"))
miRDB = as.matrix(read.table("MirTarget4.txt", head = TRUE, sep = "\t"))
microT = as.matrix(read.table("DIANA_microT_CDS.txt", head = TRUE, sep = "\t"))
miRanda = as.matrix(read.table("miRanda_mirSVR.txt", head = TRUE, sep = "\t"))
miRmap = as.matrix(read.table("miRmap.txt", head = TRUE, sep = "\t"))
miRWalk = as.matrix(read.table("miRWalk3.txt", head = TRUE, sep = "\t"))

###################################################################################################
######################################### MAIN code ###############################################
###################################################################################################

CWCS_all_mod = pair_set_mod(CWCS_all, symbol_to_ID)
CWCS_cons_mod = pair_set_mod(CWCS_cons, symbol_to_ID)
qMRE_all_mod = pair_set_mod(qMRE_all, symbol_to_ID)
qMRE_cons_mod = pair_set_mod(qMRE_cons, symbol_to_ID)
PITA_mod = pair_set_mod(PITA, symbol_to_ID)
mirtar_mod = pair_set_mod(mirtar, symbol_to_ID)
miRDB_mod = pair_set_mod(miRDB, symbol_to_ID)
microT_mod = pair_set_mod(microT, symbol_to_ID)
miRanda_mod = pair_set_mod(miRanda, symbol_to_ID)
miRmap_mod = pair_set_mod(miRmap, symbol_to_ID)
miRWalk_mod = pair_set_mod(miRWalk, symbol_to_ID)

CWCS_all_matrix1 = matrix_transfer(CWCS_all_mod)
CWCS_cons_matrix1 = matrix_transfer(CWCS_cons_mod)
qMRE_all_matrix1 = matrix_transfer(qMRE_all_mod)
qMRE_cons_matrix1 = matrix_transfer(qMRE_cons_mod)
PITA_matrix1 = matrix_transfer(PITA_mod)
mirtar_matrix1 = matrix_transfer(mirtar_mod)
miRDB_matrix1 = matrix_transfer(miRDB_mod)
microT_matrix1 = matrix_transfer(microT_mod)
miRanda_matrix1 = matrix_transfer(miRanda_mod)
miRmap_matrix1 = matrix_transfer(miRmap_mod)
miRWalk_matrix1 = matrix_transfer(miRWalk_mod)

########################## microarray ##############################################################

mRNA_matrix = as.matrix(read.table("HeLa_mRNA_array.txt", head = TRUE, sep = "\t"))
miRNA_matrix = as.matrix(read.table("HeLa_miRNA_array.txt", head = TRUE, sep = "\t"))

CWCS_all_matrix = filter_exp(CWCS_all_matrix1, miRNA_matrix, mRNA_matrix)
CWCS_cons_matrix = filter_exp(CWCS_cons_matrix1, miRNA_matrix, mRNA_matrix)
qMRE_all_matrix = filter_exp(qMRE_all_matrix1, miRNA_matrix, mRNA_matrix)
qMRE_cons_matrix = filter_exp(qMRE_cons_matrix1, miRNA_matrix, mRNA_matrix)
PITA_matrix = filter_exp(PITA_matrix1, miRNA_matrix, mRNA_matrix)
mirtar_matrix = filter_exp(mirtar_matrix1, miRNA_matrix, mRNA_matrix)
miRDB_matrix = filter_exp(miRDB_matrix1, miRNA_matrix, mRNA_matrix)
microT_matrix = filter_exp(microT_matrix1, miRNA_matrix, mRNA_matrix)
miRanda_matrix = filter_exp(miRanda_matrix1, miRNA_matrix, mRNA_matrix)
miRmap_matrix = filter_exp(miRmap_matrix1, miRNA_matrix, mRNA_matrix)
miRWalk_matrix = filter_exp(miRWalk_matrix1, miRNA_matrix, mRNA_matrix)

miracle_all_result = miracle_score(miRNA_matrix, mRNA_matrix, CWCS_all_matrix)
miracle_cons_result = miracle_score(miRNA_matrix, mRNA_matrix, CWCS_cons_matrix)
promise_all_mRNA = promise_score_mRNA(miRNA_matrix, mRNA_matrix, qMRE_all_matrix)
promise_all_joint = promise_score_joint(miRNA_matrix, mRNA_matrix, qMRE_all_matrix)
promise_cons_mRNA = promise_score_mRNA(miRNA_matrix, mRNA_matrix, qMRE_cons_matrix)
promise_cons_joint = promise_score_joint(miRNA_matrix, mRNA_matrix, qMRE_cons_matrix)

num = 5000

miracle_all_accu = accu_pred(miracle_all_result, pos_list, num)
miracle_cons_accu = accu_pred(miracle_cons_result, pos_list, num)
promise_all_x_accu = accu_pred(promise_all_mRNA, pos_list, num)
promise_all_xz_accu = accu_pred(promise_all_joint, pos_list, num)
promise_cons_x_accu = accu_pred(promise_cons_mRNA, pos_list, num)
promise_cons_xz_accu = accu_pred(promise_cons_joint, pos_list, num)
CWCS_all_accu = accu_pred(CWCS_all_matrix, pos_list, num)
CWCS_cons_accu = accu_pred(CWCS_cons_matrix, pos_list, num)
PITA_accu = accu_pred(PITA_matrix, pos_list, num)
mirtar_accu = accu_pred(mirtar_matrix, pos_list, num)
miRDB_accu = accu_pred(miRDB_matrix, pos_list, num)
microT_accu = accu_pred(microT_matrix, pos_list, num)
miRanda_accu = accu_pred(miRanda_matrix, pos_list, num)
miRmap_accu = accu_pred(miRmap_matrix, pos_list, num)
miRWalk_accu = accu_pred(miRWalk_matrix, pos_list, num)

percent_accu = cbind(miracle_all_accu[[2]], miracle_cons_accu[[2]], promise_all_x_accu[[2]], promise_all_xz_accu[[2]], 
                     promise_cons_x_accu[[2]], promise_cons_xz_accu[[2]], CWCS_all_accu[[2]], CWCS_cons_accu[[2]],
                     PITA_accu[[2]], mirtar_accu[[2]], miRDB_accu[[2]], 
                     microT_accu[[2]], miRanda_accu[[2]], miRmap_accu[[2]], miRWalk_accu[[2]])

colnames(percent_accu) = c('miRACLe', 'miRACLe_cons', 'ProMISe_mrna', 'ProMISe_joint',
                          'ProMISe_mrna_cons', 'ProMISe_joint_cons', 'TargetScan7_CWCS', 'TargetScan7_CWCS_cons', 'PITA',
                          'miRTar2GO_HeLa', 'MirTarget4', 'DIANA-microT-CDS', 'miRanda-mirSVR', 'miRmap', 'miRWalk3')
rownames(percent_accu) = seq(100, num, 100)

test_result = wil_test(percent_accu)

write.table(percent_accu, file = "microarray - precision.txt", quote = FALSE, sep = "\t")
write.table(test_result, file = "microarray - p-value of Wilcoxon test.txt", quote = FALSE, sep = "\t")

########################## RNA-seq ##############################################################

mRNA_matrix = as.matrix(read.table("HeLa_mRNA_Seq.txt", head = TRUE, sep = "\t"))
miRNA_matrix = as.matrix(read.table("HeLa_miRNA_Seq.txt", head = TRUE, sep = "\t"))

CWCS_all_matrix = filter_exp(CWCS_all_matrix1, miRNA_matrix, mRNA_matrix)
CWCS_cons_matrix = filter_exp(CWCS_cons_matrix1, miRNA_matrix, mRNA_matrix)
qMRE_all_matrix = filter_exp(qMRE_all_matrix1, miRNA_matrix, mRNA_matrix)
qMRE_cons_matrix = filter_exp(qMRE_cons_matrix1, miRNA_matrix, mRNA_matrix)
PITA_matrix = filter_exp(PITA_matrix1, miRNA_matrix, mRNA_matrix)
mirtar_matrix = filter_exp(mirtar_matrix1, miRNA_matrix, mRNA_matrix)
miRDB_matrix = filter_exp(miRDB_matrix1, miRNA_matrix, mRNA_matrix)
microT_matrix = filter_exp(microT_matrix1, miRNA_matrix, mRNA_matrix)
miRanda_matrix = filter_exp(miRanda_matrix1, miRNA_matrix, mRNA_matrix)
miRmap_matrix = filter_exp(miRmap_matrix1, miRNA_matrix, mRNA_matrix)
miRWalk_matrix = filter_exp(miRWalk_matrix1, miRNA_matrix, mRNA_matrix)

miracle_all_result = miracle_score(miRNA_matrix, mRNA_matrix, CWCS_all_matrix)
miracle_cons_result = miracle_score(miRNA_matrix, mRNA_matrix, CWCS_cons_matrix)
promise_all_mRNA = promise_score_mRNA(miRNA_matrix, mRNA_matrix, qMRE_all_matrix)
promise_all_joint = promise_score_joint(miRNA_matrix, mRNA_matrix, qMRE_all_matrix)
promise_cons_mRNA = promise_score_mRNA(miRNA_matrix, mRNA_matrix, qMRE_cons_matrix)
promise_cons_joint = promise_score_joint(miRNA_matrix, mRNA_matrix, qMRE_cons_matrix)

num = 5000

miracle_all_accu = accu_pred(miracle_all_result, pos_list, num)
miracle_cons_accu = accu_pred(miracle_cons_result, pos_list, num)
promise_all_x_accu = accu_pred(promise_all_mRNA, pos_list, num)
promise_all_xz_accu = accu_pred(promise_all_joint, pos_list, num)
promise_cons_x_accu = accu_pred(promise_cons_mRNA, pos_list, num)
promise_cons_xz_accu = accu_pred(promise_cons_joint, pos_list, num)
CWCS_all_accu = accu_pred(CWCS_all_matrix, pos_list, num)
CWCS_cons_accu = accu_pred(CWCS_cons_matrix, pos_list, num)
PITA_accu = accu_pred(PITA_matrix, pos_list, num)
mirtar_accu = accu_pred(mirtar_matrix, pos_list, num)
miRDB_accu = accu_pred(miRDB_matrix, pos_list, num)
microT_accu = accu_pred(microT_matrix, pos_list, num)
miRanda_accu = accu_pred(miRanda_matrix, pos_list, num)
miRmap_accu = accu_pred(miRmap_matrix, pos_list, num)
miRWalk_accu = accu_pred(miRWalk_matrix, pos_list, num)


percent_accu = cbind(miracle_all_accu[[2]], miracle_cons_accu[[2]], promise_all_x_accu[[2]], promise_all_xz_accu[[2]], 
                     promise_cons_x_accu[[2]], promise_cons_xz_accu[[2]], CWCS_all_accu[[2]], CWCS_cons_accu[[2]],
                     PITA_accu[[2]], mirtar_accu[[2]], miRDB_accu[[2]], 
                     microT_accu[[2]], miRanda_accu[[2]], miRmap_accu[[2]], miRWalk_accu[[2]])

colnames(percent_accu) = c('miRACLe', 'miRACLe_cons', 'ProMISe_mrna', 'ProMISe_joint',
                          'ProMISe_mrna_cons', 'ProMISe_joint_cons', 'TargetScan7_CWCS', 'TargetScan7_CWCS_cons', 'PITA',
                          'miRTar2GO_HeLa', 'MirTarget4', 'DIANA-microT-CDS', 'miRanda-mirSVR', 'miRmap', 'miRWalk3')
rownames(percent_accu) = seq(100, num, 100)

test_result = wil_test(percent_accu)

write.table(percent_accu, file = "RNA-seq - precision.txt", quote = FALSE, sep = "\t")
write.table(test_result, file = "RNA-seq - p-value of Wilcoxon test.txt", quote = FALSE, sep = "\t")
