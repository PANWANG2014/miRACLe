######################################################################################################

   ### Part I: FUNCTIONS THAT ARE USED IN DATA PROCESSING AND MODELING ###

#the following three functions are used to match the names of samples, mirnas and mrnas
#By using these three functions, we can get the matched mrna expression data and mirna expression data
name_func = function(name_base, mirna_base, mrna_base){
  x = match(name_base[, 2], mirna_base[1, ])
  y = match(name_base[, 3], mrna_base[1, ])
  return(list(x,y))
}
mirna_matrix = function(name_base, mirna_base, mirna_name, mirna){  
  mirna_name1 = rep(0, nrow(mirna_base))
  for(i in 1:nrow(mirna_base)){
    mirna_name1[i] = unlist(strsplit(mirna_base[i, 1], "\\|"))[1]
  }
  mirna_use = mirna_base[mirna_name1 %in% mirna, mirna_name]
  mirna_name1 = mirna_name1[mirna_name1 %in% mirna]
  return(list(mirna_use, mirna_name1))
}
mrna_matrix = function(name_base, mrna_base, mrna_name, mrna){
  #mirna_name1 is used to find the location of identified pairs in the final matrix
  mrna_exist = rep(0, nrow(mrna_base))
  mrna_name_all = matrix(data = 0, ncol = 2, nrow = nrow(mrna_base), byrow = TRUE)
  mrna_name1 = rep(0, nrow(mrna_base))
  
  for(i in 1:nrow(mrna_base)){
    mrna_name1[i] = unlist(strsplit(mrna_base[i, 1], "\\|"))[1]
    mrna_name_all[i, 1] = unlist(strsplit(mrna_base[i, 1], "\\|"))[1]
    mrna_name_all[i, 2] = unlist(strsplit(mrna_base[i, 1], "\\|"))[2]
    if(mrna_name1[i] %in% mrna == 1){
      mrna_exist[i] = i
    }
  }
  mrna_use = mrna_base[mrna_exist,mrna_name]
  mrna_name1 = mrna_name1[mrna_exist]
  mrna_name_sp = mrna_name_all[mrna_exist,]
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
  mrna_name_sp1 = mrna_name_sp[mrna_sgn, ]
  mrna_fullname1 = mrna_fullname[mrna_sgn]
  rownames(mirna_use1) = mirna_name1
  rownames(mrna_use1) = mrna_fullname1
  return(list(mirna_use1, mrna_use1, mirna_name1, mrna_name1, mrna_name_sp1, mrna_fullname1))
}
#mirna_mrna_loc is used to get the CWCS matrix that fits the expression data
mirna_mrna_loc = function(mirna_name, mrna_name, mirna, mrna, CWCS){
  mirna_loc = match(mirna_name, mirna)
  mrna_loc = match(mrna_name, mrna)
  CWCS = t(CWCS)
  CWCS_use = CWCS[mirna_loc, mrna_loc]
  rownames(CWCS_use) = mirna_name
  colnames(CWCS_use) = mrna_name
  return(list(mirna_loc, mrna_loc, CWCS_use))
}
#cate_part is used to do stratified sampling using clinical information, num is the number to partition the data
cate_part = function(basic, num){
  whole = seq(1, length(basic[, 3]), 1)
  cate = as.numeric(basic[, 3]) * 10 + as.numeric(basic[, 4])
  cate[which(cate == 21)] = 11
  training_pos = list()
  test_pos = list()
  uni_cate = unique(cate)
  for(i in 1:num){
    pos_temp = c()
    for(j in 1:length(uni_cate)){
      pos = which(cate == uni_cate[j])
      age = as.numeric(basic[pos, 2])
      training_age = createDataPartition(age, 1, p = 0.8, groups = 2)
      pos_temp = c(pos_temp, pos[training_age[[1]]])
    }
    training_pos[[i]] = pos_temp
    test_pos[[i]] = whole[-pos_temp]
  }
  return(list(training_pos, test_pos))
}
#part_method is used to randomly partition the samples into subsets to do 10-fold cross validation
part_method = function(label, times){
  t = seq(1, length(label), 1)
  part_strategy = matrix(data = 0, ncol = length(label), nrow = times)
  if(length(label) %% 10 == 0){
    part_num = rep(length(label)/10, 10)
  }
  else{
    k = length(label) %% 10
    part_num = c(rep(ceiling(length(label)/10), k), rep(floor(length(label)/10), (10-k)))
  }
  for(time in 1:times){
    s = sample(t, size = part_num[1], replace = FALSE)
    part_strategy[time, s] = 1
    c = s
    for(i in 2:9){
      s = sample(t[-c], size = part_num[i], replace = FALSE)
      c = c(c,s)
      part_strategy[time, s] = i
    }
    s = t[-c]
    part_strategy[time, s] = 10
  }
  return(part_strategy)
}
#match_part is used to match the expression data with the order of clinical information
match_part = function(name, basic){
  name_sub = rep(0, length(name[, 1]))
  for(i in 1:length(name_sub)){
    name_sub[i] = substr(name[i, 1], 1, nchar(name[i, 1]) - 1)
  }
  basic_sub = rep(0, length(basic[, 1]))
  for(i in 1:length(basic_sub)){
    basic_sub[i] = substr(basic[i, 1], 3, nchar(basic[i, 1]))
  }
  pair = match(basic_sub, name_sub)
  return(pair)
}
#feature_select is used to get the features for mrna expression
feature_select = function(mrna_E, mrna_T, temp){
  train_E_mrna = matrix(as.numeric(mrna_E[,temp]), nrow = nrow(mrna_E[,temp]))
  train_T_mrna = matrix(as.numeric(mrna_T[,temp]), nrow = nrow(mrna_T[,temp]))
  pv_mrna = rep(0, nrow(train_E_mrna))
  for(i in 1:nrow(train_E_mrna)){
    pv_mrna[i] = t.test(as.numeric(train_E_mrna[i,]), as.numeric(train_T_mrna[i,]), paired = TRUE)$p.value
  }
  mean_E = rowMeans(train_E_mrna)
  mean_T = rowMeans(train_T_mrna)
  pos1 = which(pv_mrna < 0.05/length(pv_mrna))
  pos2 = which(abs(mean_T - mean_E) > 1)
  pos = intersect(pos1, pos2)
  mean_dif = abs(mean_T[pos] - mean_E[pos])
  pos_new = pos[order(mean_dif, decreasing = TRUE)]
  return(pos_new)
}
#para_sep is used to get the parameters of seperate version of naive bayes
para_sep = function(mrna_E, mrna_T, temp, fc_order){
  exp_E = matrix(as.numeric(mrna_E[,temp]), nrow = nrow(mrna_E[,temp]))
  exp_T = matrix(as.numeric(mrna_T[,temp]), nrow = nrow(mrna_T[,temp]))
  
  para_E_sep = matrix(data = 0, ncol = 4, nrow = length(fc_order))
  para_T_sep = matrix(data = 0, ncol = 4, nrow = length(fc_order))
  cut_point = matrix(data = 0, ncol = 3, nrow = length(fc_order))
  exp = cbind(exp_E, exp_T)
  for(i in 1:length(fc_order)){
    cut_point[i,] = as.vector(quantile(exp[fc_order[i],]))[2:4]
    para_E_sep[i,1] = (length(which(exp_E[fc_order[i], ] <= cut_point[i, 1])) + 1)/ (ncol(exp_E) + 4)
    para_E_sep[i,2] = (length(which(exp_E[fc_order[i], ] <= cut_point[i, 2])) - 
                         length(which(exp_E[fc_order[i], ] <= cut_point[i, 1])) + 1)/(ncol(exp_E) + 4)
    para_E_sep[i,3] = (length(which(exp_E[fc_order[i], ] <= cut_point[i, 3])) - 
                         length(which(exp_E[fc_order[i], ] <= cut_point[i, 2])) + 1)/(ncol(exp_E) + 4)
    para_E_sep[i,4] = (length(which(exp_E[fc_order[i], ] > cut_point[i, 3])) + 1)/ (ncol(exp_E) + 4)
    para_T_sep[i,1] = (length(which(exp_T[fc_order[i], ] <= cut_point[i, 1])) + 1)/ (ncol(exp_T) + 4)
    para_T_sep[i,2] = (length(which(exp_T[fc_order[i], ] <= cut_point[i, 2])) - 
                         length(which(exp_T[fc_order[i], ] <= cut_point[i, 1])) + 1)/(ncol(exp_T) + 4)
    para_T_sep[i,3] = (length(which(exp_T[fc_order[i], ] <= cut_point[i, 3])) - 
                         length(which(exp_T[fc_order[i], ] <= cut_point[i, 2])) + 1)/(ncol(exp_T) + 4)
    para_T_sep[i,4] = (length(which(exp_T[fc_order[i], ] > cut_point[i, 3])) + 1)/ (ncol(exp_T) + 4)
  }
  return(list(para_E_sep, para_T_sep, cut_point))
}
#test_sep is used to calculate the log probability of the test data when using seperate version of naive bayes
test_sep = function(para_sep_mrna, mrna, temp, pos_mrna){
  test_mrna = matrix(as.numeric(mrna[,temp]), nrow = nrow(mrna[,temp]))
  prob_test = matrix(data = 0, nrow = 2, ncol = length(temp))
  rownames(prob_test) = c("prob of E", "prob of T")
  for(j in 1:ncol(test_mrna)){
    for(i in 1:length(pos_mrna)){
      temp = test_mrna[pos_mrna[i], j]
      pos_fea = (temp > para_sep_mrna[[3]][i, 1])+(temp > para_sep_mrna[[3]][i, 2])+(temp > para_sep_mrna[[3]][i, 3])+1
      prob_test[1, j] = prob_test[1, j] + log(para_sep_mrna[[1]][i, pos_fea])
      prob_test[2, j] = prob_test[2, j] + log(para_sep_mrna[[2]][i, pos_fea])
    }
  }
  return(prob_test)
}
#index_cal is used to calculate related indexes of predicted results
index_cal = function(top_feature, feature_detail_mrna, accu_sep_E, accu_sep_T){
  summ_result = matrix(data = 0, ncol = 9, nrow = length(top_feature))
  colnames(summ_result) = c("chosen feature number", "intersect feature number", "union feature number", "stability score",
                            "ACC", "PPV", "TPR", "TNR", "F score")
  #the following code is used to calculate the stablility score
  temp_inter = list()
  temp_uni = list()
  for(j in 1:length(top_feature)){
    temp_inter[[j]] =  feature_detail_mrna[[1]][1:top_feature[j],1]
    temp = feature_detail_mrna[[1]][1:top_feature[j], 1]
    for(k in 1:10){
      temp = intersect(temp, feature_detail_mrna[[time]][1:top_feature[j], k])
    }
    temp_uni[[j]] = temp
    for(time in 1:times){
      temp = feature_detail_mrna[[time]][1:top_feature[j], 1]
      for(k in 1:10){
        temp = intersect(temp, feature_detail_mrna[[time]][1:top_feature[j], k])
      }
      temp_inter[[j]] = intersect(temp_inter[[j]], temp)
      temp_uni[[j]] = union(temp_uni[[j]], temp)
    }
  }
  #the following code is used to calculate different indexes
  for(i in 1:length(top_feature)){
    summ_result[i, 2] = length(temp_inter[[i]])
    summ_result[i, 3] = length(temp_uni[[i]])
    summ_result[i, 4] = length(temp_inter[[i]])/length(temp_uni[[i]])
  }
  summ_result[, 1] = top_feature
  summ_result[, 5] = (rowSums(accu_sep_E) + rowSums(accu_sep_T))/(rowSums(total_val_vali) + rowSums(total_val_vali))
  summ_result[, 6] = rowSums(accu_sep_T)/(rowSums(total_val_vali) - rowSums(accu_sep_E) + rowSums(accu_sep_T))
  summ_result[, 7] = rowSums(accu_sep_T)/rowSums(total_val_vali)
  summ_result[, 8] = rowSums(accu_sep_E)/rowSums(total_val_vali)
  summ_result[, 9] = 2 * summ_result[, 6] * summ_result[, 7] / (summ_result[, 6] + summ_result[, 7])
  return(summ_result)
}

######################################################################################################

   ### Part II: DATA INPUT ###

#the following package is used to do stratified sampling
library(caret)
#first read in the files we need to use, which contains:
#basic_information: information about samples
#expression data: mirna_CICAMS_T, mrna_CICAMS_T, mirna_CICAMS_E, mrna_CICAMS_E
#RNA information: mirna list, mrna list
basic_information = as.matrix(read.table("clinical information.txt", head = TRUE, sep = "\t"))
mirna_CICAMS_T = as.matrix(read.table("ESCC_Tumor_miRNA_expression.txt", head = FALSE, sep = "\t"))
mrna_CICAMS_T = as.matrix(read.table("ESCC_Tumor_mRNA_expression.txt", head = FALSE, sep = "\t"))
name_CICAMS_T = as.matrix(read.table("ESCC_Tumor_miRNA_mRNA_paired.txt", head = TRUE, sep = "\t"))
mirna_CICAMS_E = as.matrix(read.table("ESCC_Normal_miRNA_expression.txt", head = FALSE, sep = "\t"))
mrna_CICAMS_E = as.matrix(read.table("ESCC_Normal_mRNA_expression.txt", head = FALSE, sep = "\t"))
name_CICAMS_E = as.matrix(read.table("ESCC_Normal_miRNA_mRNA_paired.txt", head = TRUE, sep = "\t"))
mrna = as.matrix(read.table("mrna_list.txt", head = TRUE, sep = "\t"))
mirna = as.matrix(read.table("mirna_list.txt", head = TRUE, sep = "\t"))

######################################################################################################

   ### Part III: MAIN PROGRAM ###

#the following code is used for preparation of the analysis
x_E = name_func(name_CICAMS_E, mirna_CICAMS_E, mrna_CICAMS_E)
z1_E = mirna_matrix(name_CICAMS_E, mirna_CICAMS_E, x_E[[1]], mirna)
z2_E = mrna_matrix(name_CICAMS_E, mrna_CICAMS_E, x_E[[2]], mrna)
z3_E = mirna_mrna_data(z1_E[[1]], z2_E[[1]], z1_E[[2]], z2_E[[2]], z2_E[[3]], z2_E[[4]], 1)
x_T = name_func(name_CICAMS_T, mirna_CICAMS_T, mrna_CICAMS_T)
z1_T = mirna_matrix(name_CICAMS_T, mirna_CICAMS_T, x_T[[1]], mirna)
z2_T = mrna_matrix(name_CICAMS_T, mrna_CICAMS_T, x_T[[2]], mrna)
z3_T = mirna_mrna_data(z1_T[[1]], z2_T[[1]], z1_T[[2]], z2_T[[2]], z2_T[[3]], z2_T[[4]], 1)

#match the order of expression data with clinical information, and partition the data into training data and test data
times = 100
part = cate_part(basic_information, times) 
compare = match_part(name_CICAMS_E, basic_information)
test_record = matrix(unlist(part[[2]]), nrow = length(part[[2]][[1]]))
t = seq(1, length(part[[1]][[1]]), 1)
part_strategy = part_method(part[[1]][[1]], times)
top_feature = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100) #these are the numbers of features we select
feature_number_mrna = matrix(data = 0, ncol = 10, nrow = times) #used to store the features selected in each CV round
feature_detail_mrna = list()
for(i in 1:times){
  feature_detail_mrna[[i]] = matrix(data = 0, ncol = 10, nrow = 1000)
}
feature_vali_number = matrix(data = 0, ncol = times, nrow = length(top_feature)) #used to record number of features in validation step
accu_sep_T = matrix(data = 0, ncol = times, nrow = length(top_feature))
accu_sep_E = matrix(data = 0, ncol = times, nrow = length(top_feature))
total_val_vali = matrix(data = 0, ncol = times, nrow = length(top_feature))
for(time in 1:times){
  #partition the data into training data and validation data
  mrna_E = z3_E[[2]][,compare[part[[1]][[time]]]]
  mrna_T = z3_T[[2]][,compare[part[[1]][[time]]]]
  vali_mrna_E = z3_E[[2]][,compare[part[[2]][[time]]]]
  vali_mrna_T = z3_T[[2]][,compare[part[[2]][[time]]]]
  #the following code is used to do the cross validation and select the features from mrna
  for(k in 1:10){
    temp2 = which(part_strategy[time,] == k) #for test
    temp1 = t[-temp2] #for training
    pos_mrna = feature_select(mrna_E, mrna_T, temp1) #pos_mrna is the ranked features
    feature_number_mrna[time, k] = length(pos_mrna)
    feature_detail_mrna[[time]][, k] = pos_mrna[1:1000]
  }
  #the following part is used to do the internal independent validation
  temp3 = seq(1, ncol(mrna_E), 1)
  temp4 = seq(1, ncol(vali_mrna_E), 1)
  for(j in 1:length(top_feature)){
    fea_mrna = feature_detail_mrna[[time]][1:top_feature[j], 1]
    for(k in 1:10){
      fea_mrna = intersect(fea_mrna, feature_detail_mrna[[time]][1:top_feature[j], k])
    }
    vali_sep_mrna = para_sep(mrna_E, mrna_T, temp3, fea_mrna)
    prob_E_sep_vali = test_sep(vali_sep_mrna, vali_mrna_E, temp4, fea_mrna)
    prob_T_sep_vali = test_sep(vali_sep_mrna, vali_mrna_T, temp4, fea_mrna)
    accu_sep_E[j, time] = accu_sep_E[j, time] + sum(prob_E_sep_vali[1,] >= prob_E_sep_vali[2,])
    accu_sep_T[j, time] = accu_sep_T[j, time] + sum(prob_T_sep_vali[1,] < prob_T_sep_vali[2,])
    total_val_vali[j, time] = total_val_vali[j, time] + ncol(vali_mrna_E)
  }
  print(time)
}
summ_result = index_cal(top_feature, feature_detail_mrna, accu_sep_E, accu_sep_T)
#the following code is used to get the information of the features ranked
feature_detail = feature_detail_mrna[[1]]
for(time in 2:times){
  feature_detail = cbind(feature_detail, feature_detail_mrna[[time]])
}
feature_name = matrix(data = 0, ncol = ncol(feature_detail), nrow = nrow(feature_detail))
feature_ID = matrix(data = 0, ncol = ncol(feature_detail), nrow = nrow(feature_detail))
for(i in 1:ncol(feature_name)){
  feature_name[, i] = z3_E[[5]][feature_detail[, i], 1]
  feature_ID[, i] = z3_E[[5]][feature_detail[, i], 2]
}

write.table(summ_result, file = "summary result of NB_mrna.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(feature_name, file = "selected mrna symbol.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(feature_ID, file = "selected mrna ID.txt", quote = FALSE, sep = "\t", row.names = FALSE)

