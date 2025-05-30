
# -------------------------------------------------------------------------

library(ISLR2)
data(NCI60)
dim(NCI60$data)
table(NCI60$labs)

sel <- c("BREAST", "COLON", "MELANOMA", "NSCLC", "RENAL")
nci.data <- NCI60$data[(NCI60$labs %in% sel), ]
nci.labs <- NCI60$labs[(NCI60$labs %in% sel)]

dim(nci.data)
length(nci.labs)
table(nci.labs)

data_scaled = scale(nci.data)

rownames(data_scaled) = as.character(nci.labs)
colnames(data_scaled) = colnames(nci.data)

data_scaled[1:5, 1:5]


# 1 -----------------------------------------------------------------------

ans_1 = list()

for (subtype in sel){
  
  sub_data = data_scaled[rownames(data_scaled) == subtype, ]
  
  dist_mat = dist(sub_data, method = "euclidean")
  
  min_dist = min(dist_mat)
  median_dist = median(dist_mat)
  max_dist = max(dist_mat)
  
  ans_1[[subtype]] = data.frame(
    Min_dist = min_dist,
    Median_dist = median_dist,
    Max_dist = max_dist
  )
}

ans_1


# 2 -----------------------------------------------------------------------

dist_data = dist(data_scaled, method = "euclidean")

linkage_types = c("single", "complete", "average", "centroid") 

ans_2 = list()

par(mfrow = c(2,2))

for (type in linkage_types){
  
  hclust_result = hclust(dist_data, method = type)
  
  plot(hclust_result, 
       labels = rownames(data_scaled), 
       main = paste(type, "linkage"), ylab = "Dist", sub = NA, cex = 0.6)
  
  clusters = cutree(hclust_result, k = 5)
  
  comp_table = table(clusters, nci.labs)
  ans_2[[type]] = comp_table
}

ans_2

# Explain the reason of your choice
#'4가지 linkage 방법(single, complete, average, centroid) 중에서 
#'complete linkage 방법을 선택하는 것이 가장 타당하다.
#'single, average, centroid 방법 모두 1번 Cluster에 거의 모든 subtype이 몰려 있다.
#'single linkage의 1번 Vluster에는
#'BREAST 5개, COLON 7개, MELANOMA 8개, NSCLC 8개, RENAL 8개로
#'모든 subtype이 한 Cluster에 섞여 있는 나쁜 결과를 보인다. 
#'하지만 상대적으로 complete linkage는 비교적 subtype별로 분리되는 경향을 보인다
#'Cluster 3: MELANOMA 3개 -> MELANOMA 중심
#'Cluster 5: NSCLC 6개 -> NSCLC 중
#'Cluster 1: RENAL 5개 -> RENAL 중심
#'나머지 Cluster들도 특정 subtype 위주로 분포되어 있어 상대적으로 해석이 가능해진다.
#'정리하면 다른 방법들은 한 Cluster에 모든 subtype이 뒤섞이는 반면,
#'complete linkage는 subtype 간 분리가 상대적으로 잘 되어 있다.


# 3 -----------------------------------------------------------------------

sel = c("BREAST", "COLON", "MELANOMA", "NSCLC", "RENAL")

ans_3 = matrix(1, nrow=5, ncol=5, dimnames = list(sel, sel))

for (i in 1:5){
  for (j in 1:5){

    A = data_scaled[rownames(data_scaled) == sel[i], ]
    B = data_scaled[rownames(data_scaled) == sel[j], ]
    
    #' cor-> 열 별
    #' 우리는 샘플 cor니깐 transpose rr
    pcc_mat = cor(t(A), t(B))
    
    max_abs_corr = pcc_mat[which.max(abs(pcc_mat))]
    
    ans_3[i,j] = max_abs_corr
  }
}

ans_3


# 4 -----------------------------------------------------------------------

factor = as.factor(nci.labs)
pvals = NULL

#anova(lm(data_scaled[,1] ~ factor))$Pr[1]
pvals = apply(data_scaled, 2, function(data_) anova(lm(data_~factor))$Pr[1])

top_20_gene = as.integer(names(sort(pvals)[1:20]))

top_20_data = data_scaled[ ,top_20_gene]

set.seed(1234)

k_means_ = kmeans(top_20_data, 5, nstart = 20)
W = k_means_$tot.withinss

W_j = NULL

for (i in 1:ncol(top_20_data)){
  
  W_j[i] = kmeans(top_20_data[, -i], 5, nstart = 20)$tot.withinss
  
}

ans_4 = data.frame(
  gene_Num = colnames(top_20_data),
  WW_j = abs(W-W_j)
)

ans_4_sorted = ans_4[order(-ans_4$WW_j), ]

ans_4 = list(
  List_of_WW_j = ans_4_sorted,
  Most_sign_gene = ans_4_sorted$gene_Num[1]
)

ans_4


# 5 -----------------------------------------------------------------------

set.seed(1234)
boot <- apply(matrix(1:40, 40, 1000), 2, sample, replace=TRUE)

pve = numeric(1000)
for (i in 1:1000){
  
  index = boot[, i]
  
  X_b = data_scaled[index, ]
  
  C = cor(t(X_b))                     
  
  eig = eigen(C)
  
  values = eig$values
  
  pve[i] = sum(values[1:2]) / sum(values)   
  
}

CI = quantile(pve, c(0.025, 0.975))
CI


# 6-------------------------------------------------------------------------

pca = prcomp(data_scaled, scale=TRUE)
PC_scores = pca$x[, 1:20]  # 40x20

factor_ = as.factor(nci.labs)
                    
F_value = apply(PC_scores, 2, function(pc) anova(lm(pc ~ factor_))$F[1])

set.seed(1234)
perm <- apply(matrix(1:40, 40, 1000), 2, sample, replace=FALSE)

F_perm = matrix(NA, nrow=20, ncol=1000)


for (i in 1:1000) {
  
  factor_perm = factor_[perm[, i]]
  
  F_perm[, i] = apply(PC_scores, 2, function(pc) anova(lm(pc ~ factor_perm))$F[1])
}

#(sum(F_perm[1, ] >= F_value[1]) + 1) / (1000 + 1)
perm_pval = (rowSums(F_perm >= F_value) + 1) / (1000 + 1)


sig_idx = which(perm_pval < 0.05)
ans_6 = data.frame(
  PC = paste0("PC", sig_idx),
  perm_pval = perm_pval[sig_idx]
)

ans_6








































