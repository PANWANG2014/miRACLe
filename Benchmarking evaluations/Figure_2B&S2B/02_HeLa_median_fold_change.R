###################################################################################################
## This Rscript generates benchmarking evaluation of HeLa - median fold change

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

method_medfc = function(method_score, mrna, mirna, vset, num){
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
  z10 = method_medfc(miracle_all_result[[3]], miracle_all_result[[2]], miracle_all_result[[1]], vset, num) #miracle all model
  z11 = method_medfc(miracle_conserved_result[[3]], miracle_conserved_result[[2]], miracle_conserved_result[[1]], vset, num) #miracle conserved model
  z12 = method_medfc(promise_all_result[[3]], promise_all_result[[2]], promise_all_result[[1]], vset, num) #promise all_mrna competition model
  z13 = method_medfc(promise_all_result[[4]], promise_all_result[[2]], promise_all_result[[1]], vset, num) #promise all_full competition model
  z14 = method_medfc(promise_conserved_result[[3]], promise_conserved_result[[2]], promise_conserved_result[[1]], vset, num) #promise conserved_mrna competition model
  z15 = method_medfc(promise_conserved_result[[4]], promise_conserved_result[[2]], promise_conserved_result[[1]], vset, num) #promise conserved_full competition model
  z16 = method_medfc(CWCS_all_matrix[[3]], CWCS_all_matrix[[2]], CWCS_all_matrix[[1]], vset, num) #CWCS all
  z17 = method_medfc(CWCS_conserved_matrix[[3]], CWCS_conserved_matrix[[2]], CWCS_conserved_matrix[[1]], vset, num) #CWCS conserved
  z18 = method_medfc(PITA_matrix[[3]], PITA_matrix[[2]], PITA_matrix[[1]], vset, num) #for PITA
  z19 = method_medfc(mirtar_matrix[[3]], mirtar_matrix[[2]], mirtar_matrix[[1]], vset, num) #for mirtar
  z20 = method_medfc(miRDB_matrix[[3]], miRDB_matrix[[2]], miRDB_matrix[[1]], vset, num) #for miRDB
  z21 = method_medfc(microT_matrix[[3]], microT_matrix[[2]], microT_matrix[[1]], vset, num) #for microT
  z22 = method_medfc(miRanda_matrix[[3]], miRanda_matrix[[2]], miRanda_matrix[[1]], vset, num) #for miRanda
  z23 = method_medfc(miRmap_matrix[[3]], miRmap_matrix[[2]], miRmap_matrix[[1]], vset, num) #for miRmap
  z24 = method_medfc(miRWalk_matrix[[3]], miRWalk_matrix[[2]], miRWalk_matrix[[1]], vset, num) #for miRWalk
  
  z27 = medfc_test(vset, miracle_all_result[[1]], z10[[1]], z10[[3]], num)
  z28 = medfc_vector(z27[[1]], z27[[2]]) #for miracle_all
  z29 = medfc_test(vset, miracle_conserved_result[[1]], z11[[1]], z11[[3]], num)
  z30 = medfc_vector(z29[[1]], z29[[2]]) #for miracle_conserved
  z31 = medfc_test(vset, promise_all_result[[1]], z12[[1]], z12[[3]], num)
  z32 = medfc_vector(z31[[1]], z31[[2]]) #for promise all_mrna model
  z33 = medfc_test(vset, promise_all_result[[1]], z13[[1]], z13[[3]], num)
  z34 = medfc_vector(z33[[1]], z33[[2]]) #for promise all_joint model
  z35 = medfc_test(vset, promise_conserved_result[[1]], z14[[1]], z14[[3]], num)
  z36 = medfc_vector(z35[[1]], z35[[2]]) #for promise conserved_mrna model
  z37 = medfc_test(vset, promise_conserved_result[[1]], z15[[1]], z15[[3]], num)
  z38 = medfc_vector(z37[[1]], z37[[2]]) #for promise conserved_joint model
  z39 = medfc_test(vset, CWCS_all_matrix[[1]], z16[[1]], z16[[3]], num)
  z40 = medfc_vector(z39[[1]], z39[[2]]) #for CWCS all
  z41 = medfc_test(vset, CWCS_conserved_matrix[[1]], z17[[1]], z17[[3]], num)
  z42 = medfc_vector(z41[[1]], z41[[2]]) #for CWCS conserved
  z43 = medfc_test(vset, PITA_matrix[[1]], z18[[1]], z18[[3]], num)
  z44 = medfc_vector(z43[[1]], z43[[2]]) #for PITA
  z45 = medfc_test(vset, mirtar_matrix[[1]], z19[[1]], z19[[3]], num)
  z46 = medfc_vector(z45[[1]], z45[[2]]) #for mirtar
  z47 = medfc_test(vset, miRDB_matrix[[1]], z20[[1]], z20[[3]], num)
  z48 = medfc_vector(z47[[1]], z47[[2]]) #for miRDB 
  z49 = medfc_test(vset, microT_matrix[[1]], z21[[1]], z21[[3]], num)
  z50 = medfc_vector(z49[[1]], z49[[2]]) #for microT
  z51 = medfc_test(vset, miRanda_matrix[[1]], z22[[1]], z22[[3]], num)
  z52 = medfc_vector(z51[[1]], z51[[2]]) #for miranda 
  z53 = medfc_test(vset, miRmap_matrix[[1]], z23[[1]], z23[[3]], num)
  z54 = medfc_vector(z53[[1]], z53[[2]]) #for miRmap
  z55 = medfc_test(vset, miRWalk_matrix[[1]], z24[[1]], z24[[3]], num)
  z56 = medfc_vector(z55[[1]], z55[[2]]) #for miRWalk
  
  if(length(z44) == 0){
    final_median = rbind(z28[1:500], z30[1:500], z32[1:500], z34[1:500], z36[1:500], z38[1:500], z40[1:500], z42[1:500],  
                         z46[1:500], z48[1:500], z50[1:500], z52[1:500], z54[1:500], z56[1:500])
    colnames(final_median) = seq(1, 500, 1)
    rownames(final_median) = c('miRACLe', 'miRACLe_cons', 'ProMISe_mrna', 'ProMISe_joint', 
                               'ProMISe_mrna_cons', 'ProMISe_joint_cons', 'TargetScan7_CWCS', 'TargetScan7_CWCS_cons',
                               'miRTar2GO_HeLa', 'MirTarget4', 'DIANA-microT-CDS', 'miRanda-mirSVR', 'miRmap', 'miRWalk3')
  }
  else{
    final_median = rbind(z28[1:500], z30[1:500], z32[1:500], z34[1:500], z36[1:500], z38[1:500], z40[1:500], z42[1:500],  
                         z44[1:500], z46[1:500], z48[1:500], z50[1:500], z52[1:500], z54[1:500], z56[1:500])
    colnames(final_median) = seq(1, 500, 1)
    rownames(final_median) = c('miRACLe', 'miRACLe_cons', 'ProMISe_mrna', 'ProMISe_joint',
                          'ProMISe_mrna_cons', 'ProMISe_joint_cons', 'TargetScan7_CWCS', 'TargetScan7_CWCS_cons', 'PITA',
                          'miRTar2GO_HeLa', 'MirTarget4', 'DIANA-microT-CDS', 'miRanda-mirSVR', 'miRmap', 'miRWalk3')
  }

  return(final_median)
}

wil_test = function(valid){
  if(nrow(valid) == 15){
    result = matrix(data = 0, ncol = 2, nrow = 13)
    for(i in 3:15){
      result[(i - 2), 1] = wilcox.test(as.vector(valid[1, ]), as.vector(valid[i, ]), alternative = "less", paired = TRUE, correct = FALSE)$p.value
      result[(i - 2), 2] = wilcox.test(as.vector(valid[2, ]), as.vector(valid[i, ]), alternative = "less", paired = TRUE, correct = FALSE)$p.value
    }
    colnames(result) = c('miRACLe', 'miRACLe_cons')
    rownames(result) = c('ProMISe_mrna', 'ProMISe_joint', 'ProMISe_mrna_cons', 'ProMISe_joint_cons', 'TargetScan7_CWCS', 
                       'TargetScan7_CWCS_cons', 'PITA', 'miRTar2GO_HeLa', 'MirTarget4', 'DIANA-microT-CDS','miRanda-mirSVR', 
                       'miRmap', 'miRWalk3')
  }
  else if(nrow(valid) == 14){
    result = matrix(data = 0, ncol = 2, nrow = 12)
    for(i in 3:14){
      result[(i - 2), 1] = wilcox.test(as.vector(valid[1, ]), as.vector(valid[i, ]), alternative = "less", paired = TRUE, correct = FALSE)$p.value
      result[(i - 2), 2] = wilcox.test(as.vector(valid[2, ]), as.vector(valid[i, ]), alternative = "less", paired = TRUE, correct = FALSE)$p.value
    }
    colnames(result) = c('miRACLe', 'miRACLe_cons')
    rownames(result) = c('ProMISe_mrna', 'ProMISe_joint', 'ProMISe_mrna_cons', 'ProMISe_joint_cons', 'TargetScan7_CWCS', 
                         'TargetScan7_CWCS_cons', 'miRTar2GO_HeLa', 'MirTarget4', 'DIANA-microT-CDS','miRanda-mirSVR', 
                         'miRmap', 'miRWalk3')
  }
  return(result)
}

###################################################################################################
####################################### DATA INPUTS ###############################################
###################################################################################################

symbol_to_ID = as.matrix(read.table("Symbol_to_ID.txt", head = TRUE, sep = "\t"))
CWCS_all = as.matrix(read.table("TargetScan7_CWCS.txt", head = TRUE, sep = "\t"))
CWCS_conserved = as.matrix(read.table("TargetScan7_CWCS_cons.txt", head = TRUE, sep = "\t"))
qMRE_all = as.matrix(read.table("TargetScan7_qMRE.txt", head = TRUE, sep = "\t"))
qMRE_conserved = as.matrix(read.table("TargetScan7_qMRE_cons.txt", head = TRUE, sep = "\t"))
PITA = as.matrix(read.table("PITA.txt", head = TRUE, sep = "\t"))
mirtar = as.matrix(read.table("miRTar2GO_HeLa.txt", head = TRUE, sep = "\t"))
miRDB = as.matrix(read.table("MirTarget4.txt", head = TRUE, sep = "\t"))
microT = as.matrix(read.table("DIANA_microT_CDS.txt", head = TRUE, sep = "\t"))
miRanda = as.matrix(read.table("miRanda_mirSVR.txt", head = TRUE, sep = "\t"))
miRmap = as.matrix(read.table("miRmap.txt", head = TRUE, sep = "\t"))
miRWalk = as.matrix(read.table("miRWalk3.txt", head = TRUE, sep = "\t"))

Liu = as.matrix(read.table("Transet_HeLa_Seq.txt", head = TRUE, sep = "\t"))
Selbach = as.matrix(read.table("Transet_HeLa_Array.txt", head = TRUE, sep = "\t"))


###################################################################################################
######################################### MAIN code ###############################################
###################################################################################################

CWCS_all_mod = pair_set_mod(CWCS_all, symbol_to_ID)
CWCS_conserved_mod = pair_set_mod(CWCS_conserved, symbol_to_ID)
qMRE_all_mod = pair_set_mod(qMRE_all, symbol_to_ID)
qMRE_conserved_mod = pair_set_mod(qMRE_conserved, symbol_to_ID)
PITA_mod = pair_set_mod(PITA, symbol_to_ID)
mirtar_mod = pair_set_mod(mirtar, symbol_to_ID)
miRDB_mod = pair_set_mod(miRDB, symbol_to_ID)
microT_mod = pair_set_mod(microT, symbol_to_ID)
miRanda_mod = pair_set_mod(miRanda, symbol_to_ID)
miRmap_mod = pair_set_mod(miRmap, symbol_to_ID)
miRWalk_mod = pair_set_mod(miRWalk, symbol_to_ID)

CWCS_all_matrix1 = matrix_transfer(CWCS_all_mod)
CWCS_conserved_matrix1 = matrix_transfer(CWCS_conserved_mod)
qMRE_all_matrix1 = matrix_transfer(qMRE_all_mod)
qMRE_conserved_matrix1 = matrix_transfer(qMRE_conserved_mod)
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
CWCS_conserved_matrix = filter_exp(CWCS_conserved_matrix1, miRNA_matrix, mRNA_matrix)
qMRE_all_matrix = filter_exp(qMRE_all_matrix1, miRNA_matrix, mRNA_matrix)
qMRE_conserved_matrix = filter_exp(qMRE_conserved_matrix1, miRNA_matrix, mRNA_matrix)
PITA_matrix = filter_exp(PITA_matrix1, miRNA_matrix, mRNA_matrix)
mirtar_matrix = filter_exp(mirtar_matrix1, miRNA_matrix, mRNA_matrix)
miRDB_matrix = filter_exp(miRDB_matrix1, miRNA_matrix, mRNA_matrix)
microT_matrix = filter_exp(microT_matrix1, miRNA_matrix, mRNA_matrix)
miRanda_matrix = filter_exp(miRanda_matrix1, miRNA_matrix, mRNA_matrix)
miRmap_matrix = filter_exp(miRmap_matrix1, miRNA_matrix, mRNA_matrix)
miRWalk_matrix = filter_exp(miRWalk_matrix1, miRNA_matrix, mRNA_matrix)

miracle_all_result = miracle_score(miRNA_matrix, mRNA_matrix, CWCS_all_matrix)
miracle_conserved_result = miracle_score(miRNA_matrix, mRNA_matrix, CWCS_conserved_matrix)
promise_all_result = promise_score(miRNA_matrix, mRNA_matrix, qMRE_all_matrix)
promise_conserved_result = promise_score(miRNA_matrix, mRNA_matrix, qMRE_conserved_matrix)

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

CWCS_all_matrix = filter_exp(CWCS_all_matrix1, miRNA_matrix, mRNA_matrix)
CWCS_conserved_matrix = filter_exp(CWCS_conserved_matrix1, miRNA_matrix, mRNA_matrix)
qMRE_all_matrix = filter_exp(qMRE_all_matrix1, miRNA_matrix, mRNA_matrix)
qMRE_conserved_matrix = filter_exp(qMRE_conserved_matrix1, miRNA_matrix, mRNA_matrix)
PITA_matrix = filter_exp(PITA_matrix1, miRNA_matrix, mRNA_matrix)
mirtar_matrix = filter_exp(mirtar_matrix1, miRNA_matrix, mRNA_matrix)
miRDB_matrix = filter_exp(miRDB_matrix1, miRNA_matrix, mRNA_matrix)
microT_matrix = filter_exp(microT_matrix1, miRNA_matrix, mRNA_matrix)
miRanda_matrix = filter_exp(miRanda_matrix1, miRNA_matrix, mRNA_matrix)
miRmap_matrix = filter_exp(miRmap_matrix1, miRNA_matrix, mRNA_matrix)
miRWalk_matrix = filter_exp(miRWalk_matrix1, miRNA_matrix, mRNA_matrix)

miracle_all_result = miracle_score(miRNA_matrix, mRNA_matrix, CWCS_all_matrix)
miracle_conserved_result = miracle_score(miRNA_matrix, mRNA_matrix, CWCS_conserved_matrix)
promise_all_result = promise_score(miRNA_matrix, mRNA_matrix, qMRE_all_matrix) 
promise_conserved_result = promise_score(miRNA_matrix, mRNA_matrix, qMRE_conserved_matrix)

final_Liu = medfc_summary(Liu)
final_Selbach = medfc_summary(Selbach)
test_Liu = wil_test(final_Liu)
test_Selbach = wil_test(final_Selbach)

write.table(final_Liu, file = "RNA-Seq - mfc based on seq transet.txt", quote = FALSE, sep = "\t")
write.table(final_Selbach, file = "RNA-Seq - mfc based on array transet.txt", quote = FALSE, sep = "\t")
write.table(test_Liu, file = "RNA-Seq - p-value of Wilcoxon test for seq transet.txt", quote = FALSE, sep = "\t")
write.table(test_Selbach, file = "RNA-Seq - p-value of Wilcoxon test for array transet.txt", quote = FALSE, sep = "\t")
