###################################################################################################
## This Rscript generates benchmarking evaluation of NCI60 - precision analysis

###################################################################################################
## Loading required packages
library(Roleswitch)

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

result_summary = function(name_cancer, mirna_cancer, mrna_cancer, score_matrix, Vset, thres_num, method){
  score_list = matrix_transfer(score_matrix)
  x = name_func(name_cancer, mirna_cancer, mrna_cancer)
  z1 = mirna_matrix(mirna_cancer, x[[1]], score_list[[1]])
  z2 = mrna_matrix(mrna_cancer, x[[2]], score_list[[2]])
  if(method == "miracle"){
    z3 = mirna_mrna_data(z1[[1]], z2[[1]], z1[[2]], z2[[2]], z2[[3]], z2[[4]], 1)
    z5 = MMI_location(Vset, z3[[3]], z3[[5]]) 
    z7 = mirna_mrna_loc(z3[[3]], z3[[4]], score_list[[1]], score_list[[2]], score_list[[3]], z3[[6]])
    z8 = sum_miracle(z3[[1]], z3[[2]], z7[[3]]) 
    z9 = miracle_sam(z3[[1]], z3[[2]], z7[[3]], z8, thres_num, z5) #result of miracle
    return(z9)
  }
  else if(method == "promise"){
    z3 = mirna_mrna_data(z1[[1]], z2[[1]], z1[[2]], z2[[2]], z2[[3]], z2[[4]], 1)
    z5 = MMI_location(Vset, z3[[3]], z3[[5]]) 
    z6 = mirna_mrna_loc(z3[[3]], z3[[4]], score_list[[1]], score_list[[2]], score_list[[3]], z3[[6]])
    z10 = promise_sam(z3[[1]], z3[[2]], z6[[3]], thres_num, z5) #result of promise
    return(z10)
  }
  else if(method == "sequence"){
    z4 = mirna_mrna_data(z1[[1]], z2[[1]], z1[[2]], z2[[2]], z2[[3]], z2[[4]], 0.8)
    z13 = MMI_location(Vset, z4[[3]], z4[[5]])
    z14 = mirna_mrna_loc(z4[[3]], z4[[4]], score_list[[1]], score_list[[2]], score_list[[3]], z4[[6]])
    z12 = targetscan_sam(z14[[3]], thres_num, z13) #result of sequence method
    return(z12)
  }
  
}

diff_rate = function(valid){
  valid_new = rep(0, length(valid))
  valid_new[1] = valid[1]
  for(j in 2:length(valid_new)){
    valid_new[j] = valid[j] * j - valid[(j - 1)] * (j - 1)
  }
  return(valid_new)
}


wil_test_pre = function(str1, str2){
  str1_new = diff_rate(str1)
  str2_new = diff_rate(str2)
  value = wilcox.test(str1_new, str2_new, alternative = "greater", paired = TRUE, correct = FALSE)$p.value
  return(value)
}

wil_test = function(miracle_all, miracle_cons, promise_all, promise_cons, targetscan_CWCS, 
                    targetscan_CWCS_cons, PITA_all, mirtar2go_all, mirTarget4_all, microT_all,
                    miRanda_all, miRmap_all, miRWalk_all){
  result_all = matrix(data = 0, nrow = ncol(miracle_all), ncol = 13)
  result_con = matrix(data = 0, nrow = ncol(miracle_all), ncol = 13)
  for(i in 1:nrow(result_all)){
    result_all[i, 1] = wil_test_pre(miracle_all[, i], promise_all[[1]][, i])
    result_all[i, 2] = wil_test_pre(miracle_all[, i], promise_all[[2]][, i])
    result_all[i, 3] = wil_test_pre(miracle_all[, i], promise_cons[[1]][, i])
    result_all[i, 4] = wil_test_pre(miracle_all[, i], promise_cons[[2]][, i])
    result_all[i, 5] = wil_test_pre(miracle_all[, i], targetscan_CWCS)
    result_all[i, 6] = wil_test_pre(miracle_all[, i], targetscan_CWCS_cons)
    result_all[i, 7] = wil_test_pre(miracle_all[, i], PITA_all)
    result_all[i, 8] = wil_test_pre(miracle_all[, i], mirtar2go_all)
    result_all[i, 9] = wil_test_pre(miracle_all[, i], mirTarget4_all)
    result_all[i, 10] = wil_test_pre(miracle_all[, i], microT_all)
    result_all[i, 11] = wil_test_pre(miracle_all[, i], miRanda_all)
    result_all[i, 12] = wil_test_pre(miracle_all[, i], miRmap_all)
    result_all[i, 13] = wil_test_pre(miracle_all[, i], miRWalk_all)
    
    result_con[i, 1] = wil_test_pre(miracle_cons[, i], promise_all[[1]][, i])
    result_con[i, 2] = wil_test_pre(miracle_cons[, i], promise_all[[2]][, i])
    result_con[i, 3] = wil_test_pre(miracle_cons[, i], promise_cons[[1]][, i])
    result_con[i, 4] = wil_test_pre(miracle_cons[, i], promise_cons[[2]][, i])
    result_con[i, 5] = wil_test_pre(miracle_cons[, i], targetscan_CWCS)
    result_con[i, 6] = wil_test_pre(miracle_cons[, i], targetscan_CWCS_cons)
    result_con[i, 7] = wil_test_pre(miracle_cons[, i], PITA_all)
    result_con[i, 8] = wil_test_pre(miracle_cons[, i], mirtar2go_all)
    result_con[i, 9] = wil_test_pre(miracle_cons[, i], mirTarget4_all)
    result_con[i, 10] = wil_test_pre(miracle_cons[, i], microT_all)
    result_con[i, 11] = wil_test_pre(miracle_cons[, i], miRanda_all)
    result_con[i, 12] = wil_test_pre(miracle_cons[, i], miRmap_all)
    result_con[i, 13] = wil_test_pre(miracle_cons[, i], miRWalk_all)
  }
  
  colnames(result_all) = c('ProMISe_mrna', 'ProMISe_joint', 'ProMISe_mrna_cons', 'ProMISe_joint_cons', 
                           'TargetScan_CWCS', 'TargetScan_CWCS_cons', 'PITA', 'miRTar2GO', 'MirTarget4', 'DIANA-microT-CDS', 
                           'miRanda-mirSVR', 'miRmap', 'miRWalk3')
  colnames(result_con) = c('ProMISe_mrna', 'ProMISe_joint', 'ProMISe_mrna_cons', 'ProMISe_joint_cons', 
                           'TargetScan_CWCS', 'TargetScan_CWCS_cons', 'PITA', 'miRTar2GO', 'MirTarget4', 'DIANA-microT-CDS', 
                           'miRanda-mirSVR', 'miRmap', 'miRWalk3')
  return(list(result_all, result_con))
}

###################################################################################################
####################################### DATA INPUTS ###############################################
###################################################################################################

CWCS_all = as.matrix(read.table("TargetScan7_CWCS.txt", head = TRUE, sep = "\t"))
CWCS_conserved = as.matrix(read.table("TargetScan7_CWCS_cons.txt", head = TRUE, sep = "\t"))
qMRE_all = as.matrix(read.table("TargetScan7_qMRE.txt", head = TRUE, sep = "\t"))
qMRE_conserved = as.matrix(read.table("TargetScan7_qMRE_cons.txt", head = TRUE, sep = "\t"))
PITA = as.matrix(read.table("PITA.txt", head = TRUE, sep = "\t"))
mirtar = as.matrix(read.table("miRTar2GO.txt", head = TRUE, sep = "\t"))
miRDB = as.matrix(read.table("MirTarget4.txt", head = TRUE, sep = "\t"))
microT = as.matrix(read.table("DIANA_microT_CDS.txt", head = TRUE, sep = "\t"))
miRanda = as.matrix(read.table("miRanda_mirSVR.txt", head = TRUE, sep = "\t"))
miRmap = as.matrix(read.table("miRmap.txt", head = TRUE, sep = "\t"))
miRWalk = as.matrix(read.table("miRWalk3.txt", head = TRUE, sep = "\t"))

Vset = read.table("Vset_Celllines.txt", head = TRUE, sep = "\t")
name_NCI = as.matrix(read.table("NCI60_sampleMatch.txt", head = TRUE, sep = "\t"))
mirna_NCI = as.matrix(read.table("NCI60_miRNA_expression.txt", head = FALSE, sep = "\t"))
mrna_NCI = as.matrix(read.table("NCI60_mRNA_expression.txt",head = FALSE, sep = "\t"))

###################################################################################################
######################################### MAIN code ###############################################
###################################################################################################

num = 5000
miracle_all = result_summary(name_NCI, mirna_NCI, mrna_NCI, CWCS_all, Vset, num, "miracle")
miracle_cons = result_summary(name_NCI, mirna_NCI, mrna_NCI, CWCS_conserved, Vset, num, "miracle")
promise_all = result_summary(name_NCI, mirna_NCI, mrna_NCI, qMRE_all, Vset, num, "promise")
promise_cons = result_summary(name_NCI, mirna_NCI, mrna_NCI, qMRE_conserved, Vset, num, "promise")
targetscan_CWCS = result_summary(name_NCI, mirna_NCI, mrna_NCI, CWCS_all, Vset, num, "sequence")
targetscan_CWCS_cons = result_summary(name_NCI, mirna_NCI, mrna_NCI, CWCS_conserved, Vset, num, "sequence")
PITA_all = result_summary(name_NCI, mirna_NCI, mrna_NCI, PITA, Vset, num, "sequence")
mirtar2go_all = result_summary(name_NCI, mirna_NCI, mrna_NCI, mirtar, Vset, num, "sequence")
mirTarget4_all = result_summary(name_NCI, mirna_NCI, mrna_NCI, miRDB, Vset, num, "sequence")
microT_all = result_summary(name_NCI, mirna_NCI, mrna_NCI, microT, Vset, num, "sequence")
miRanda_all = result_summary(name_NCI, mirna_NCI, mrna_NCI, miRanda, Vset, num, "sequence")
miRmap_all = result_summary(name_NCI, mirna_NCI, mrna_NCI, miRmap, Vset, num, "sequence")
miRWalk_all = result_summary(name_NCI, mirna_NCI, mrna_NCI, miRWalk, Vset, num, "sequence")
colnames(miracle_all) = name_NCI[, 1]
colnames(miracle_cons) = name_NCI[, 1]
colnames(promise_all[[1]]) = name_NCI[, 1]
colnames(promise_all[[2]]) = name_NCI[, 1]
colnames(promise_cons[[1]]) = name_NCI[, 1]
colnames(promise_cons[[2]]) = name_NCI[, 1]

rownames(miracle_all) = seq(100, num, 100)
rownames(miracle_cons) = seq(100, num, 100)
rownames(promise_all[[1]]) = seq(100, num, 100)
rownames(promise_all[[2]]) = seq(100, num, 100)
rownames(promise_cons[[1]]) = seq(100, num, 100)
rownames(promise_cons[[2]]) = seq(100, num, 100)

wil_result = wil_test(miracle_all, miracle_cons, promise_all, promise_cons, targetscan_CWCS, targetscan_CWCS_cons, PITA_all, 
                      mirtar2go_all, mirTarget4_all, microT_all, miRanda_all, miRmap_all, miRWalk_all)

rownames(wil_result[[1]]) = name_NCI[, 1]
rownames(wil_result[[2]]) = name_NCI[, 1]

write.table(miracle_all, file = "Precision for miRACLe.txt", quote = FALSE,sep = "\t")
write.table(miracle_cons, file = "Precision for miRACLe_cons.txt", quote = FALSE,sep = "\t")
write.table(promise_all[[1]], file = "Precision for ProMISe_mrna.txt", quote = FALSE, sep = "\t")
write.table(promise_all[[2]], file = "Precision for ProMISe_joint.txt", quote = FALSE, sep = "\t")
write.table(promise_cons[[1]], file = "Precision for ProMISe_mrna_cons.txt", quote = FALSE, sep = "\t")
write.table(promise_cons[[2]], file = "Precision for ProMISe_joint_cons.txt", quote = FALSE, sep = "\t")
write.table(targetscan_CWCS, file = "Precision for TargetScan_CWCS.txt", quote = FALSE, sep = "\t")
write.table(targetscan_CWCS_cons, file = "Precision for TargetScan_CWCS_cons.txt", quote = FALSE, sep = "\t")
write.table(PITA_all, file = "Precision for PITA.txt", quote = FALSE, sep = "\t")
write.table(mirtar2go_all, file = "Precision for miRTar2GO.txt", quote = FALSE, sep = "\t")
write.table(mirTarget4_all, file = "Precision for MirTarget4.txt", quote = FALSE, sep = "\t")
write.table(microT_all, file = "Precision for DIANA_microT_CDS.txt", quote = FALSE, sep = "\t")
write.table(miRanda_all, file = "Precision for miRanda_mirSVR.txt", quote = FALSE, sep = "\t")
write.table(miRmap_all, file = "Precision for miRmap.txt", quote = FALSE, sep = "\t")
write.table(miRWalk_all, file = "Precision for miRWalk3.txt", quote = FALSE, sep = "\t")

write.table(wil_result[[1]], file = "p-value of Wilcoxon test for miRACLe.txt", quote = FALSE, sep = "\t")
write.table(wil_result[[2]], file = "p-value of Wilcoxon test for miRACLe_cons.txt", quote = FALSE, sep = "\t")
