###################################################################################################
## This Rscript generates benchmarking evaluation of TCGA - additional precision analysis

###################################################################################################
######################################## FUNCTIONS ################################################
###################################################################################################

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
  return(list(mrna_use, mrna_name1))
}

mirna_mrna_data = function(mirna_use, mrna_use, mirna_name, mrna_name, cutoff){
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
  rownames(mirna_use1) = mirna_name1
  rownames(mrna_use1) = mrna_name1
  return(list(mirna_use1, mrna_use1, mirna_name1, mrna_name1))
}

exp_fit_method = function(name_cancer, mirna_cancer, mrna_cancer, score_list){
  x = name_func(name_cancer, mirna_cancer, mrna_cancer)
  z1 = mirna_matrix(name_cancer, mirna_cancer, x[[1]], score_list[[1]])
  z2 = mrna_matrix(name_cancer, mrna_cancer, x[[2]], score_list[[2]])
  z3 = mirna_mrna_data(z1[[1]], z2[[1]], z1[[2]], z2[[2]], 1)
  return(z3)
}

pair_set_mod = function(pair_set, symbol_to_ID){
  pos = match(pair_set[, 1], symbol_to_ID[, 1])
  pos1 = na.omit(pos)
  pair_set_filter = pair_set[which(is.na(pos) == FALSE), ]
  pair_set_filter[, 1] = symbol_to_ID[pos1, 2]
  return(pair_set_filter)
}

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
  temp1_pos = match(pos_list[, 1], mRNA_list)
  temp2_pos = match(pos_list[, 2], miRNA_list)
  temp3_pos = intersect(which(is.na(temp1_pos) == FALSE), which(is.na(temp2_pos) == FALSE))
  
  pos_list_mod = pos_list[temp3_pos, ]
  temp4_pos = match(pos_list_mod[, 1], mRNA_list)
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

valid_trans = function(vset, trans_ID){
  vsett = vset[which(as.numeric(vset[, 1]) %in% as.numeric(trans_ID[, 2])),]
  vsett[, 1] = trans_ID[match(as.numeric(vsett[, 1]), as.numeric(trans_ID[, 2])), 1]
  return(vsett)
}

filter_exp = function(seq_matrix, exp_matrix){
  temp1 = which(seq_matrix[[1]] %in% exp_matrix[[3]])
  temp2 = which(seq_matrix[[2]] %in% exp_matrix[[4]])
  filter_miRNA1 = seq_matrix[[1]][temp1]
  filter_miRNA = filter_miRNA1[match(exp_matrix[[3]], filter_miRNA1)]
  filter_mRNA1 = seq_matrix[[2]][temp2]
  filter_mRNA = filter_mRNA1[match(exp_matrix[[4]], filter_mRNA1)]
  filter_matrix1 = seq_matrix[[3]][temp2, temp1]
  filter_matrix = filter_matrix1[match(exp_matrix[[4]], filter_mRNA1), match(exp_matrix[[3]], filter_miRNA1)]
  return(list(filter_miRNA, filter_mRNA, filter_matrix))
}

miracle_pop = function(exp_trans, CWCS, sumup){
  CWCS = t(CWCS)
  mirna_use = matrix(as.numeric(exp_trans[[1]]), nrow = nrow(exp_trans[[1]]))
  mrna_use = matrix(as.numeric(exp_trans[[2]]), nrow = nrow(exp_trans[[2]]))
  pro_matrix_temp = matrix(data = 0, ncol = nrow(mrna_use), nrow = nrow(mirna_use), byrow = TRUE)
  pro_matrix_agg = matrix(data = 0, ncol = nrow(mrna_use), nrow = nrow(mirna_use), byrow = TRUE)
  for(m in 1:ncol(mirna_use)){
    pro_matrix_temp = (mirna_use[, m] %*% t(mrna_use[, m])) * CWCS / sumup[m]
    pro_matrix_agg = pro_matrix_agg + pro_matrix_temp
  }
  pro_matrix = t(pro_matrix_agg)
  #order_agg = order(pro_matrix_agg, decreasing = TRUE)[1 : thres_num]
  return(list(exp_trans[[3]], exp_trans[[4]], pro_matrix))
}

sum_miracle = function(mirna, mrna, CWCS){
  CWCS = t(CWCS)
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

###################################################################################################
####################################### DATA INPUTS ###############################################
###################################################################################################

symbol_to_ID = as.matrix(read.table("Symbol to ID.txt", head = TRUE, sep = "\t"))
Vset = as.matrix(read.table("Vset_all.txt", head = TRUE, sep = "\t"))

miRDB = as.matrix(read.table("MirTarget4.txt", head = TRUE, sep = "\t"))
microT = as.matrix(read.table("DIANA_microT_CDS.txt", head = TRUE, sep = "\t"))
miRanda = as.matrix(read.table("miRanda_mirSVR.txt", head = TRUE, sep = "\t"))

name_cancer = as.matrix(read.table("TCGA_BRCA_sampleMatch.txt", head = TRUE, sep = "\t"))
mirna_cancer = as.matrix(read.table("BRCA_miRNA_expression.txt", head = FALSE, sep = "\t"))
mrna_cancer = as.matrix(read.table("BRCA_mRNA_expression.txt", head = FALSE, sep = "\t"))

###################################################################################################
######################################### MAIN code ###############################################
###################################################################################################

V1_trans = valid_trans(Vset, symbol_to_ID)

miRDB_matrix1 = matrix_transfer(miRDB)
microT_matrix1 = matrix_transfer(microT)
miRanda_matrix1 = matrix_transfer(miRanda)

exp_trans_miRDB = exp_fit_method(name_cancer, mirna_cancer, mrna_cancer, miRDB_matrix1)
exp_trans_microT = exp_fit_method(name_cancer, mirna_cancer, mrna_cancer, microT_matrix1)
exp_trans_miRanda = exp_fit_method(name_cancer, mirna_cancer, mrna_cancer, miRanda_matrix1)

miRDB_matrix = filter_exp(miRDB_matrix1, exp_trans_miRDB)
microT_matrix = filter_exp(microT_matrix1, exp_trans_microT)
miRanda_matrix = filter_exp(miRanda_matrix1, exp_trans_miRanda)

num = 5000
sumup_miRDB = sum_miracle(exp_trans_miRDB[[1]], exp_trans_miRDB[[2]], miRDB_matrix[[3]])
miracle_miRDB_TCGA = miracle_pop(exp_trans_miRDB, miRDB_matrix[[3]], sumup_miRDB)
miracle_miRDB_V1 = accu_pred(miracle_miRDB_TCGA, V1_trans, num)

sumup_microT = sum_miracle(exp_trans_microT[[1]], exp_trans_microT[[2]], microT_matrix[[3]])
miracle_microT_TCGA = miracle_pop(exp_trans_microT, microT_matrix[[3]], sumup_microT)
miracle_microT_V1 = accu_pred(miracle_microT_TCGA, V1_trans, num)

sumup_miRanda = sum_miracle(exp_trans_miRanda[[1]], exp_trans_miRanda[[2]], miRanda_matrix[[3]])
miracle_miRanda_TCGA = miracle_pop(exp_trans_miRanda, miRanda_matrix[[3]], sumup_miRanda)
miracle_miRanda_V1 = accu_pred(miracle_miRanda_TCGA, V1_trans, num)

miRDB_V1 = accu_pred(miRDB_matrix, V1_trans, num)
microT_V1 = accu_pred(microT_matrix, V1_trans, num)
miRanda_V1 = accu_pred(miRanda_matrix, V1_trans, num)
percent_V1 = cbind(miracle_miRDB_V1[[2]], miracle_microT_V1[[2]], miracle_miRanda_V1[[2]], 
                   miRDB_V1[[2]], microT_V1[[2]], miRanda_V1[[2]])

rownames(percent_V1) = seq(100, num, 100)
colnames(percent_V1) = c('miracle_mirTarget4', 'miracle_DIANA-microT-CDS', 'miracle_miRanda-mirSVR', 
                         'mirTarget4', 'DIANA-microT-CDS', 'miRanda-mirSVR')

write.table(percent_V1, file = "BRCA_additional_precision_Vall.txt", quote = FALSE, sep = "\t")

