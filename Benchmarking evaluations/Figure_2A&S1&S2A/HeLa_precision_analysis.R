###################################################################################################
## This Rscript generates benchmarking evaluation of Hela - precision analysis

###################################################################################################
## Installing and loading required packages, also loading required data
library(Roleswitch)

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
## FUNCTIONS

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

###################################################################################################
## MAIN code

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
miracle_all_tr = top_ranker(miracle_all_result, num)
miracle_cons_tr = top_ranker(miracle_cons_result, num)
promise_all_x_tr = top_ranker(promise_all_mRNA, num)
promise_all_xz_tr = top_ranker(promise_all_joint, num)
promise_cons_x_tr = top_ranker(promise_cons_mRNA, num)
promise_cons_xz_tr = top_ranker(promise_cons_joint, num)
CWCS_all_tr = top_ranker(CWCS_all_matrix, num)
CWCS_cons_tr = top_ranker(CWCS_cons_matrix, num)
PITA_tr = top_ranker(PITA_matrix, num)
mirtar_hela_tr = top_ranker(mirtar_matrix, num)
miRDB_tr = top_ranker(miRDB_matrix, num)
microT_tr = top_ranker(microT_matrix, num)
miRanda_tr = top_ranker(miRanda_matrix, num)
miRmap_tr = top_ranker(miRmap_matrix, num)
miRWalk_tr = top_ranker(miRWalk_matrix, num)

write.table(miracle_all_tr, file = "microarray_Toprankers_miRACLe.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(miracle_cons_tr, file = "microarray_Toprankers_miRACLe_cons.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(promise_all_x_tr, file = "microarray_Toprankers_ProMISe_mrna.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(promise_all_xz_tr, file = "microarray_Toprankers_ProMISe_joint.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(promise_cons_x_tr, file = "microarray_Toprankers_ProMISe_mrna_cons.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(promise_cons_xz_tr, file = "microarray_Toprankers_ProMISe_joint_cons.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(CWCS_all_tr, file = "microarray_Toprankers_TargetScan_CWCS.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(CWCS_cons_tr, file = "microarray_Toprankers_TargetScan_CWCS_cons.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(PITA_tr, file = "microarray_Toprankers_PITA.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(mirtar_hela_tr, file = "microarray_Toprankers_miRTar2GO_HeLa.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(miRDB_tr, file = "microarray_Toprankers_MirTarget4.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(microT_tr, file = "microarray_Toprankers_DIANA_microT_CDS.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(miRanda_tr, file = "microarray_Toprankers_miRanda_mirSVR.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(miRmap_tr, file = "microarray_Toprankers_miRmap.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(miRWalk_tr, file = "microarray_Toprankers_miRWalk3.txt", quote = FALSE, sep = "\t", row.names = FALSE)

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

number_accu = cbind(miracle_all_accu[[1]], miracle_cons_accu[[1]], promise_all_x_accu[[1]], promise_all_xz_accu[[1]], 
                    promise_cons_x_accu[[1]], promise_cons_xz_accu[[1]], CWCS_all_accu[[1]], CWCS_cons_accu[[1]], 
                    PITA_accu[[1]], mirtar_accu[[1]], miRDB_accu[[1]], 
                    microT_accu[[1]], miRanda_accu[[1]], miRmap_accu[[1]], miRWalk_accu[[1]])

colnames(number_accu) = c('miRACLe', 'miRACLe_cons', 'ProMISe_mrna', 'ProMISe_joint',
                          'ProMISe_mrna_cons', 'ProMISe_joint_cons', 'TargetScan_CWCS', 'TargetScan_CWCS_cons', 'PITA',
                          'miRTar2GO_HeLa', 'MirTarget4', 'DIANA_microT_CDS', 'miRanda_mirSVR', 'miRmap', 'miRWalk3')

percent_accu = cbind(miracle_all_accu[[2]], miracle_cons_accu[[2]], promise_all_x_accu[[2]], promise_all_xz_accu[[2]], 
                     promise_cons_x_accu[[2]], promise_cons_xz_accu[[2]], CWCS_all_accu[[2]], CWCS_cons_accu[[2]],
                     PITA_accu[[2]], mirtar_accu[[2]], miRDB_accu[[2]], 
                     microT_accu[[2]], miRanda_accu[[2]], miRmap_accu[[2]], miRWalk_accu[[2]])

colnames(percent_accu) = c('miRACLe', 'miRACLe_cons', 'ProMISe_mrna', 'ProMISe_joint',
                          'ProMISe_mrna_cons', 'ProMISe_joint_cons', 'TargetScan_CWCS', 'TargetScan_CWCS_cons', 'PITA',
                          'miRTar2GO_HeLa', 'MirTarget4', 'DIANA_microT_CDS', 'miRanda_mirSVR', 'miRmap', 'miRWalk3')

write.table(number_accu, file = "microarray_validation_count.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(percent_accu, file = "microarray_precision.txt", quote = FALSE, sep = "\t", row.names = FALSE)

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
miracle_all_tr = top_ranker(miracle_all_result, num)
miracle_cons_tr = top_ranker(miracle_cons_result, num)
promise_all_x_tr = top_ranker(promise_all_mRNA, num)
promise_all_xz_tr = top_ranker(promise_all_joint, num)
promise_cons_x_tr = top_ranker(promise_cons_mRNA, num)
promise_cons_xz_tr = top_ranker(promise_cons_joint, num)
CWCS_all_tr = top_ranker(CWCS_all_matrix, num)
CWCS_cons_tr = top_ranker(CWCS_cons_matrix, num)
PITA_tr = top_ranker(PITA_matrix, num)
mirtar_hela_tr = top_ranker(mirtar_matrix, num)
miRDB_tr = top_ranker(miRDB_matrix, num)
microT_tr = top_ranker(microT_matrix, num)
miRanda_tr = top_ranker(miRanda_matrix, num)
miRmap_tr = top_ranker(miRmap_matrix, num)
miRWalk_tr = top_ranker(miRWalk_matrix, num)

write.table(miracle_all_tr, file = "Seq_Toprankers_miRACLe.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(miracle_cons_tr, file = "Seq_Toprankers_miRACLe_cons.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(promise_all_x_tr, file = "Seq_Toprankers_ProMISe_mrna.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(promise_all_xz_tr, file = "Seq_Toprankers_ProMISe_joint.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(promise_cons_x_tr, file = "Seq_Toprankers_ProMISe_mrna_cons.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(promise_cons_xz_tr, file = "Seq_Toprankers_ProMISe_joint_cons.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(CWCS_all_tr, file = "Seq_Toprankers_TargetScan_CWCS.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(CWCS_cons_tr, file = "Seq_Toprankers_TargetScan_CWCS_cons.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(PITA_tr, file = "Seq_Toprankers_PITA.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(mirtar_hela_tr, file = "Seq_Toprankers_miRTar2GO_HeLa.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(miRDB_tr, file = "Seq_Toprankers_MirTarget4.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(microT_tr, file = "Seq_Toprankers_DIANA_microT_CDS.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(miRanda_tr, file = "Seq_Toprankers_miRanda_mirSVR.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(miRmap_tr, file = "Seq_Toprankers_miRmap.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(miRWalk_tr, file = "Seq_Toprankers_miRWalk3.txt", quote = FALSE, sep = "\t", row.names = FALSE)

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

number_accu = cbind(miracle_all_accu[[1]], miracle_cons_accu[[1]], promise_all_x_accu[[1]], promise_all_xz_accu[[1]], 
                    promise_cons_x_accu[[1]], promise_cons_xz_accu[[1]], CWCS_all_accu[[1]], CWCS_cons_accu[[1]], 
                    PITA_accu[[1]], mirtar_accu[[1]], miRDB_accu[[1]], 
                    microT_accu[[1]], miRanda_accu[[1]], miRmap_accu[[1]], miRWalk_accu[[1]])

colnames(number_accu) = c('miRACLe', 'miRACLe_cons', 'ProMISe_mrna', 'ProMISe_joint',
                          'ProMISe_mrna_cons', 'ProMISe_joint_cons', 'TargetScan_CWCS', 'TargetScan_CWCS_cons', 'PITA',
                          'miRTar2GO_HeLa', 'MirTarget4', 'DIANA_microT_CDS', 'miRanda_mirSVR', 'miRmap', 'miRWalk3')

percent_accu = cbind(miracle_all_accu[[2]], miracle_cons_accu[[2]], promise_all_x_accu[[2]], promise_all_xz_accu[[2]], 
                     promise_cons_x_accu[[2]], promise_cons_xz_accu[[2]], CWCS_all_accu[[2]], CWCS_cons_accu[[2]],
                     PITA_accu[[2]], mirtar_accu[[2]], miRDB_accu[[2]], 
                     microT_accu[[2]], miRanda_accu[[2]], miRmap_accu[[2]], miRWalk_accu[[2]])

colnames(percent_accu) = c('miRACLe', 'miRACLe_cons', 'ProMISe_mrna', 'ProMISe_joint',
                          'ProMISe_mrna_cons', 'ProMISe_joint_cons', 'TargetScan_CWCS', 'TargetScan_CWCS_cons', 'PITA',
                          'miRTar2GO_HeLa', 'MirTarget4', 'DIANA_microT_CDS', 'miRanda_mirSVR', 'miRmap', 'miRWalk3')

write.table(number_accu, file = "Seq_validation_count.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(percent_accu, file = "Seq_precision.txt", quote = FALSE, sep = "\t", row.names = FALSE)