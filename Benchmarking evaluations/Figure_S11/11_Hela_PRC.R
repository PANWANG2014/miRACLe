###################################################################################################
## This Rscript generates benchmarking evaluation of HeLa - precision-recall curve

###################################################################################################
## Loading required packages
library(Roleswitch)
library(ROCR)

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

vali_pos = function(matrix, pos_list){
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
  
  pos = which(score_matrix > 0)
  poss = positive_position[which(positive_position %in% pos)]
  
  labels = rep(0, length(pos))
  labels[which(pos %in% poss)] = 1
  predictions = score_matrix[pos]
  simple = list(predictions, labels)
  names(simple) = c("predictions", "labels")
  df = data.frame(simple)
  pred = prediction(df$predictions, df$labels)
  
  return(pred)
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
####################################### DATA INPUTS ###############################################
###################################################################################################

symbol_to_ID = as.matrix(read.table("Symbol to ID.txt", head = TRUE, sep = "\t"))
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

########################## RNA-seq #################################################################
#we can choose to load the RNA-seq data or the microarray data, the following analysis are the same
#mRNA_matrix = as.matrix(read.table("HeLa_mRNA_Seq.txt", head = TRUE, sep = "\t"))
#miRNA_matrix = as.matrix(read.table("HeLa_miRNA_Seq.txt", head = TRUE, sep = "\t"))

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

pred_miracle_all = vali_pos(miracle_all_result, pos_list)
pred_miracle_cons = vali_pos(miracle_cons_result, pos_list)
pred_promise_all_mRNA = vali_pos(promise_all_mRNA, pos_list)
pred_promise_all_joint = vali_pos(promise_all_joint, pos_list)
pred_promise_cons_mRNA = vali_pos(promise_cons_mRNA, pos_list)
pred_promise_cons_joint = vali_pos(promise_cons_joint, pos_list)

pred_CWCS_all = vali_pos(CWCS_all_matrix, pos_list)
pred_CWCS_cons = vali_pos(CWCS_cons_matrix, pos_list)
pred_PITA = vali_pos(PITA_matrix, pos_list)
pred_mirtar = vali_pos(mirtar_matrix, pos_list)
pred_miRDB = vali_pos(miRDB_matrix, pos_list)
pred_microT = vali_pos(microT_matrix, pos_list)
pred_miRanda = vali_pos(miRanda_matrix, pos_list)
pred_miRmap = vali_pos(miRmap_matrix, pos_list)
pred_miRWalk = vali_pos(miRWalk_matrix, pos_list)

perf_miracle_all = performance(pred_miracle_all, "prec", "rec")
perf_miracle_cons = performance(pred_miracle_cons, "prec", "rec")
perf_promise_all_mRNA = performance(pred_promise_all_mRNA, "prec", "rec")
perf_promise_all_joint = performance(pred_promise_all_joint, "prec", "rec")
perf_promise_cons_mRNA = performance(pred_promise_cons_mRNA, "prec", "rec")
perf_promise_cons_joint = performance(pred_promise_cons_joint, "prec", "rec")

perf_CWCS_all = performance(pred_CWCS_all, "prec", "rec")
perf_CWCS_cons = performance(pred_CWCS_cons, "prec", "rec")
perf_PITA = performance(pred_PITA, "prec", "rec")
perf_mirtar = performance(pred_mirtar, "prec", "rec")
perf_miRDB = performance(pred_miRDB, "prec", "rec")
perf_microT = performance(pred_microT, "prec", "rec")
perf_miRanda = performance(pred_miRanda, "prec", "rec")
perf_miRmap = performance(pred_miRmap, "prec", "rec")
perf_miRWalk = performance(pred_miRWalk, "prec", "rec")


pdf(file = "PRC_HeLa_microarray.pdf", height = 5, width = 5)
#pdf(file = "PRC_HeLa_RNA-seq.pdf", height = 5, width = 5) #generate the corresponding plot regarding input data

plot(perf_miracle_cons, col = rgb(r=228, g=26, b=28, maxColorValue = 255), cex.lab = 1, 
     lwd = 1, main = "Precision-Recall curve for HeLa microarray", cex.main = 1, cex.axis = 1)
plot(perf_miracle_all, add = TRUE, col = rgb(r=55, g=126, b=184, maxColorValue = 255), lwd = 1)
plot(perf_promise_all_mRNA, add = TRUE, col = rgb(r=102, g=194, b=165, maxColorValue = 255), lwd = 1)
plot(perf_promise_all_joint, add = TRUE, col = rgb(r=179, g=222, b=105, maxColorValue = 255), lwd = 1)
plot(perf_promise_cons_mRNA, add = TRUE, col = rgb(r=0, g=68, b=27, maxColorValue = 255), lwd = 1)
plot(perf_promise_cons_joint, add = TRUE, col = rgb(r=77, g=175, b=74, maxColorValue = 255), lwd = 1)
plot(perf_CWCS_all, add = TRUE, col = rgb(r=222, g=139, b=249, maxColorValue = 255), lwd = 1)
plot(perf_CWCS_cons, add = TRUE, col = rgb(r=152, g=78, b=163, maxColorValue = 255), lwd = 1)

plot(perf_PITA, add = TRUE, col = rgb(r=255, g=217, b=47, maxColorValue = 255), lwd = 1)
plot(perf_mirtar, add = TRUE, col = rgb(r=223, g=194, b=125, maxColorValue = 255), lwd = 1)
plot(perf_miRDB, add = TRUE, col = rgb(r=166, g=206, b=227, maxColorValue = 255), lwd = 1)
plot(perf_microT, add = TRUE, col = rgb(r=53, g=151, b=143, maxColorValue = 255), lwd = 1)
plot(perf_miRanda, add = TRUE, col = rgb(r=166, g=86, b=40, maxColorValue = 255), lwd = 1)
plot(perf_miRmap, add = TRUE, col = rgb(r=247, g=129, b=191, maxColorValue = 255), lwd = 1)
plot(perf_miRWalk, add = TRUE, col = rgb(r=253, g=192, b=134, maxColorValue = 255), lwd = 1)
dev.off()

