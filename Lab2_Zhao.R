## 1. Load the data set.
dream <- read.delim(file = "DREAM_data.txt", sep = "\t", header = TRUE)

## 2. Clean up the data set and divide it into gene expression data and drug sensitivity data.

### Set the row names
rownames(dream) <- dream$HGNC_ID

### Remove the first column that contains names of drugs and genes
dream <- dream[, -1]

### Generate gene expression data and drug sensitivity data
gene_expression <- dream[c(6:nrow(dream)), ]
drug_sensitivity <- dream[c(1:5), ]

### Clean up the cell lines with NA values for each drug

clean_up <- function(drug_number) {
  ### Extract the cell line names for the drug of interest
  drug_sensitivity_cell_line <- drug_sensitivity[drug_number, ]
  ### Eliminate all cell lines that have NA values
  cell_lines_no_na <- drug_sensitivity_cell_line[, !is.na(drug_sensitivity_cell_line)]
}

drug_1_cell_line <- clean_up(1)
drug_2_cell_line <- clean_up(2)
drug_3_cell_line <- clean_up(3)
drug_4_cell_line <- clean_up(4)
drug_5_cell_line <- clean_up(5)

### Each of the data.frame above contains all cell lines without any NAs. They will be used in further functions and analysis. 

###############################################################################

# A. k-NN implementation

## 1. Create a kNN classifier function for each drug. 

### Everolimus(mTOR): cell_line_no_na <- drug_1_cell_line
### Disulfiram(ALDH2: cell_line_no_na <- drug_2_cell_line
### Methylglyoxol(Pyruvate): cell_line_no_na <- drug_3_cell_line
### Mebendazole(Tubulin): cell_line_no_na <- drug_4_cell_line
### 4-HC(DNA alkylator):cell_line_no_na <- drug_5_cell_line
### k: # of nearest neighbors
### unknown_gene_expression: gene expression profile of unknown sample. This profile has the same genes as the trained gene expression data. 

kNN_classifier <- function(unknown_gene_expression, cell_line_no_na, k) {
  ### Extract gene expression data based on the names of cell lines that don't have NA value
  gene_expression_matrix <- gene_expression[, colnames(cell_line_no_na)]
  ### Calculate Pearson correlation coefficient between the gene expression data from an unknown sample and all the cell lines 
  pearson_coefficient <- t(cor(unknown_gene_expression, gene_expression_matrix, use = "complete.obs"))
  ### Pick up 5 most "close" neighbors to determine the score of unknown sample
  ordered <- pearson_coefficient[order(pearson_coefficient, decreasing = TRUE), , drop = FALSE]
  scores <- as.numeric(cell_line_no_na[, row.names(ordered)[1:k]])
  ### Calculate the fraction of the k-nearest neighbors that are sensitive to the corresponding drug (fraction of kNNs that are 1)
  final_score <- sum(scores == 1)/k
  return (final_score)
}

## 2. trial...
unknown_cell_line_1 <- gene_expression[, "X184B5"]
unknown_cell_line_2 <- gene_expression[, "BT549"]
kNN_classifier(unknown_cell_line_1, drug_1_cell_line, 5)
kNN_classifier(unknown_cell_line_2, drug_1_cell_line, 5)

###############################################################################

# B. kNN performance evaluation using Leve-One-Out Cross Validation:

### Generate a function to calculate the prediction scores in LOOCV process

### Everolimus(mTOR): cell_line_no_na <- drug_1_cell_line
### Disulfiram(ALDH2): cell_line_no_na <- drug_2_cell_line
### Methylglyoxol(Pyruvate): cell_line_no_na <- drug_3_cell_line
### Mebendazole(Tubulin): cell_line_no_na <- drug_4_cell_line
### 4-HC(DNA alkylator): cell_line_no_na <- drug_5_cell_line

kNN_LOOCV <- function(cell_line_no_na, k) {
  ### Extract gene expression data based on the names of cell lines that don't have NA value
  gene_expression_matrix <- gene_expression[, colnames(cell_line_no_na)]
  ### Generate a variable called predictions that store scores from kNN classifier (length of this variable is the total number of cell lines that don't have NA values)
  predictions <- numeric(ncol(gene_expression_matrix))
  ### Calculate the scores of cell line using Leave-One-Out Cross Validation
  for (i in c(1:ncol(gene_expression_matrix))) {
    ### Set train set and label, and test set
    train_set <- gene_expression_matrix[, -i]
    train_label <- cell_line_no_na[, -i]
    test_set <- (gene_expression_matrix[, i])
    ### Calculate Pearson correlation coefficient between the gene expression data from an unknown sample and all the cell lines
    pearson_coefficient <- t(cor(test_set, train_set, use = "complete.obs"))
    ### Pick up 5 most "close" neighbors to determine the score of unknown sample
    ordered <- pearson_coefficient[order(pearson_coefficient, decreasing = TRUE), , drop = FALSE]
    scores <- as.numeric(train_label[, row.names(ordered)[1:k]])
    predictions[i] = sum(scores == 1)
  }
  return(predictions)
}

prediction_drug1 <- kNN_LOOCV(drug_1_cell_line, 5)
prediction_drug2 <- kNN_LOOCV(drug_2_cell_line, 5)
prediction_drug3 <- kNN_LOOCV(drug_3_cell_line, 5)
prediction_drug4 <- kNN_LOOCV(drug_4_cell_line, 5)
prediction_drug5 <- kNN_LOOCV(drug_5_cell_line, 5)

### Generate ROC curves for 5 different drugs.
library(pROC)
pdf(file = "ROCs_of_Five_Drug.pdf")
par(pty = "s")
roc_curve <- roc(as.numeric(drug_1_cell_line), 
                 prediction_drug1, plot = TRUE, legacy.axes = TRUE, 
                 print.auc = TRUE, percent=TRUE, 
                 print.auc.x = 20, 
                 print.auc.y = 60,
                 xlab="False Positive Percentage (1-Specificity)", 
                 ylab="True Postive Percentage (Sensitivity)", 
                 col="#e3716e", lwd=4, 
                 main = "ROC Curve of Five Drugs")
plot.roc(as.numeric(drug_2_cell_line), prediction_drug2, percent=TRUE, col="#7ac7e2", lwd=4, print.auc = TRUE, add = TRUE, print.auc.x = 20, print.auc.y = 55)
plot.roc(as.numeric(drug_3_cell_line), prediction_drug3, percent=TRUE, col="#eca680", lwd=4, print.auc = TRUE, add = TRUE, print.auc.x = 20, print.auc.y = 50)
plot.roc(as.numeric(drug_4_cell_line), prediction_drug4, percent=TRUE, col="#54beaa", lwd=4, print.auc = TRUE, add = TRUE, print.auc.x = 20, print.auc.y = 45)
plot.roc(as.numeric(drug_5_cell_line), prediction_drug5, percent=TRUE, col="#f7df87", lwd=4, print.auc = TRUE, add = TRUE, print.auc.x = 20, print.auc.y = 40)
legend("bottomright", legend=c("Everolimus", "Disulfiram", "Methylglyoxol", "Mebendazole", "4-HC"), col=c("#e3716e", "#7ac7e2", "#eca680", "#54beaa", "#f7df87"), lwd=4, bty = "n")
dev.off()

## For Question 2b
ROC_Methylglyoxol <- roc(as.numeric(drug_3_cell_line), prediction_drug3)
ROC_Methylglyoxol$sensitivities
1-ROC_Methylglyoxol$specificities

###############################################################################

# C. Exploration of parameters affecting kNN performance:

## 1. Examine how k parameters affect kNN performance in each drug.

### Everolimus(mTOR): cell_line_no_na = drug_1_cell_line; drug_number = 1
### Disulfiram(ALDH2): cell_line_no_na = drug_2_cell_line; drug_number = 2
### Methylglyoxol(Pyruvate): cell_line_no_na = drug_3_cell_line; drug_number = 3
### Mebendazole(Tubulin): cell_line_no_na = drug_4_cell_line; drug_number = 4
### 4-HC(DNA alkylator): cell_line_no_na = drug_5_cell_line; drug_number = 5

### Generate figures 

figure_of_k <- function(cell_line_no_na, drug_number) {
  
  observed <- as.numeric(cell_line_no_na)
  predicted_1 <- kNN_LOOCV(cell_line_no_na, 3)
  predicted_2 <- kNN_LOOCV(cell_line_no_na, 5)
  predicted_3 <- kNN_LOOCV(cell_line_no_na, 7)
  
  pdf(file = paste("ROC Curve of ", rownames(drug_sensitivity)[drug_number], " with different k.pdf", sep = ""))
  par(pty = "s")
  roc_curve <- roc(observed, predicted_1, plot = TRUE, legacy.axes = TRUE, 
                   print.auc = TRUE, print.auc.y = 45, percent=TRUE, 
                   xlab="False Positive Percentage (1-Specificity)", 
                   ylab="True Postive Percentage (Sensitivity)", 
                   col="#e17771", lwd=4, 
                   main = paste(paste("ROC Curve of ", rownames(drug_sensitivity)[drug_number], "with different k")))
  plot.roc(observed, predicted_2, percent=TRUE, col="#4daf4a", lwd=4, print.auc=TRUE, add=TRUE, print.auc.y=40)
  plot.roc(observed, predicted_3, percent=TRUE, col="#ad7826", lwd=4, print.auc=TRUE, add=TRUE, print.auc.y=35)
  legend("bottomright", legend=c("k=3", "k=5", "k=7"), col=c("#e17771", "#4daf4a", "#ad7826"), lwd=4)
  dev.off()
}

figure_of_k(drug_1_cell_line, 1)
figure_of_k(drug_2_cell_line, 2)
figure_of_k(drug_3_cell_line, 3)
figure_of_k(drug_4_cell_line, 4)
figure_of_k(drug_5_cell_line, 5)

## 2. kNN classifier with weighted score.

### Everolimus(mTOR): cell_line_no_na = drug_1_cell_line; drug_number = 1
### Disulfiram(ALDH2): cell_line_no_na = drug_2_cell_line; drug_number = 2
### Methylglyoxol(Pyruvate): cell_line_no_na = drug_3_cell_line; drug_number = 3
### Mebendazole(Tubulin): cell_line_no_na = drug_4_cell_line; drug_number = 4
### 4-HC(DNA alkylator): cell_line_no_na = drug_5_cell_line; drug_number = 5

kNN_LOOCV_weighted_score <- function(cell_line_no_na, k) {
  ### Extract gene expression data based on the names of cell lines that don't have NA value
  gene_expression_matrix <- gene_expression[, colnames(cell_line_no_na)]
  ### Generate a variable called predictions that store scores from kNN classifier (length of this variable is the total number of cell lines that don't have NA values)
  predictions <- numeric(ncol(gene_expression_matrix))
  ### Calculate the scores of cell line using Leave-One-Out Cross Validation
  for (i in c(1:ncol(gene_expression_matrix))) {
    ### Set train set and label, and test set
    train_set <- gene_expression_matrix[, -i]
    train_label <- cell_line_no_na[, -i]
    test_set <- (gene_expression_matrix[, i])
    ### Calculate Pearson correlation coefficient between the gene expression data from an unknown sample and all the cell lines
    pearson_coefficient <- t(cor(test_set, train_set, use = "complete.obs"))
    ### Pick up 5 most "close" neighbors to determine the score of unknown sample
    ordered <- pearson_coefficient[order(pearson_coefficient, decreasing = TRUE), , drop = FALSE]
    score <- as.numeric(train_label[, row.names(ordered)[1:k]])
    weight <- ifelse(score==1, 1, -1)
    predictions[i] = sum(ordered[1:k] * weight)
  }
  return(predictions)
}

figure_of_weighted_score_by_pROC <- function(cell_line_no_na, drug_number, k) {
  observed <- as.numeric(cell_line_no_na)
  predicted <- kNN_LOOCV_weighted_score(cell_line_no_na, k)
  pdf(file = paste("ROC Curve of ", rownames(drug_sensitivity)[drug_number], " with Weighted Scores by pROC.pdf", sep = ""))
  par(pty = "s")
  roc_curve <- roc(observed, predicted, plot = TRUE, legacy.axes = TRUE, 
                   print.auc = TRUE, print.auc.y = 45, percent=TRUE, 
                   xlab="False Positive Percentage (1-Specificity)", 
                   ylab="True Postive Percentage (Sensitivity)", 
                   col="#e17771", lwd=4)
  legend("bottomright", legend=rownames(drug_sensitivity)[drug_number], col = "#e17771", lwd=4)
  dev.off()
}

figure_of_weighted_score_by_pROC(drug_1_cell_line, 1, 5)
figure_of_weighted_score_by_pROC(drug_2_cell_line, 2, 5)
figure_of_weighted_score_by_pROC(drug_3_cell_line, 3, 5)
figure_of_weighted_score_by_pROC(drug_4_cell_line, 4, 5)
figure_of_weighted_score_by_pROC(drug_5_cell_line, 5, 5)

figure_of_weighted_score_by_ROCit <- function(cell_line_no_na, drug_number, k) {
  observed <- as.numeric(cell_line_no_na)
  predicted <- kNN_LOOCV_weighted_score(cell_line_no_na, k)
  pdf(file = paste("ROC Curve of ", rownames(drug_sensitivity)[drug_number], " with Weighted Scores by ROCit.pdf", sep = ""))
  # install.packages("ROCit")
  library(ROCit)
  ROCit_obj <- rocit(score = predicted, class = observed)
  plot(ROCit_obj)
  legend("topleft", legend = rownames(drug_sensitivity)[drug_number], col = "black", lty = 4)
  summary_roc <- summary(ROCit_obj)
  
  dev.off()
}

figure_of_weighted_score_by_ROCit(drug_1_cell_line, 1, 5)
figure_of_weighted_score_by_ROCit(drug_2_cell_line, 2, 5)
figure_of_weighted_score_by_ROCit(drug_3_cell_line, 3, 5)
figure_of_weighted_score_by_ROCit(drug_4_cell_line, 4, 5)
figure_of_weighted_score_by_ROCit(drug_5_cell_line, 5, 5)